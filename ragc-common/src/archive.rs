// Archive I/O
// Custom binary format for storing compressed genome data in streams/parts

use crate::varint::{read_varint, write_varint};
use anyhow::{Context, Result};
use std::collections::HashMap;
use std::fs::File;
use std::io::{BufReader, BufWriter, Read, Seek, SeekFrom, Write};
use std::path::Path;

/// A part within a stream (offset and size in file)
#[derive(Debug, Clone)]
struct Part {
    offset: u64,
    size: u64,
}

impl Part {
    fn new(offset: u64, size: u64) -> Self {
        Part { offset, size }
    }
}

/// A stream containing multiple parts
#[derive(Debug)]
struct Stream {
    stream_name: String,
    cur_id: usize,
    raw_size: u64,
    packed_size: u64,
    packed_data_size: u64,
    parts: Vec<Part>,
}

impl Stream {
    fn new(stream_name: String) -> Self {
        Stream {
            stream_name,
            cur_id: 0,
            raw_size: 0,
            packed_size: 0,
            packed_data_size: 0,
            parts: Vec::new(),
        }
    }
}

/// Archive for reading/writing AGC data
pub struct Archive {
    input_mode: bool,
    file: Option<File>,
    reader: Option<BufReader<File>>,
    writer: Option<BufWriter<File>>,
    f_offset: u64,
    streams: Vec<Stream>,
    stream_map: HashMap<String, usize>,
}

impl Archive {
    /// Create a new archive in input (read) mode
    pub fn new_reader() -> Self {
        Archive {
            input_mode: true,
            file: None,
            reader: None,
            writer: None,
            f_offset: 0,
            streams: Vec::new(),
            stream_map: HashMap::new(),
        }
    }

    /// Create a new archive in output (write) mode
    pub fn new_writer() -> Self {
        Archive {
            input_mode: false,
            file: None,
            reader: None,
            writer: None,
            f_offset: 0,
            streams: Vec::new(),
            stream_map: HashMap::new(),
        }
    }

    /// Open an archive file
    pub fn open<P: AsRef<Path>>(&mut self, path: P) -> Result<()> {
        if self.input_mode {
            let file = File::open(path).context("Failed to open archive for reading")?;
            self.reader = Some(BufReader::new(file.try_clone()?));
            self.file = Some(file);
            self.deserialize()?;
        } else {
            let file = File::create(path).context("Failed to create archive for writing")?;
            self.writer = Some(BufWriter::new(file.try_clone()?));
            self.file = Some(file);
        }
        self.f_offset = 0;
        Ok(())
    }

    /// Close the archive (writes footer in write mode)
    pub fn close(&mut self) -> Result<()> {
        if !self.input_mode {
            if let Some(ref mut writer) = self.writer {
                writer.flush()?;
            }
            self.serialize()?;
        }

        self.reader = None;
        self.writer = None;
        self.file = None;
        Ok(())
    }

    /// Register a new stream and return its ID
    pub fn register_stream(&mut self, stream_name: &str) -> usize {
        // Check if already registered
        if let Some(&id) = self.stream_map.get(stream_name) {
            return id;
        }

        let id = self.streams.len();
        self.streams.push(Stream::new(stream_name.to_string()));
        self.stream_map.insert(stream_name.to_string(), id);
        id
    }

    /// Get stream ID by name (returns None if not found)
    pub fn get_stream_id(&self, stream_name: &str) -> Option<usize> {
        self.stream_map.get(stream_name).copied()
    }

    /// Get list of all stream names
    pub fn get_stream_names(&self) -> Vec<String> {
        self.streams.iter().map(|s| s.stream_name.clone()).collect()
    }

    /// Add a part to a stream
    pub fn add_part(&mut self, stream_id: usize, data: &[u8], metadata: u64) -> Result<()> {
        if stream_id >= self.streams.len() {
            anyhow::bail!("Invalid stream ID: {}", stream_id);
        }

        let writer = self
            .writer
            .as_mut()
            .context("Archive not open for writing")?;

        // Record part offset (before writing anything)
        let part_offset = self.f_offset;

        // Write metadata as varint
        let mut metadata_buf = Vec::new();
        write_varint(&mut metadata_buf, metadata)?;
        writer.write_all(&metadata_buf)?;
        self.f_offset += metadata_buf.len() as u64;

        // Write data
        writer.write_all(data)?;
        self.f_offset += data.len() as u64;

        // Record part (size is only the data size, not including metadata)
        self.streams[stream_id]
            .parts
            .push(Part::new(part_offset, data.len() as u64));

        // packed_size includes both metadata and data
        let total_size = (self.f_offset - part_offset) as u64;
        self.streams[stream_id].packed_size += total_size;
        self.streams[stream_id].packed_data_size += data.len() as u64;

        Ok(())
    }

    /// Set raw (uncompressed) size for a stream
    pub fn set_raw_size(&mut self, stream_id: usize, raw_size: u64) {
        if stream_id < self.streams.len() {
            self.streams[stream_id].raw_size = raw_size;
        }
    }

    /// Get raw size for a stream
    pub fn get_raw_size(&self, stream_id: usize) -> u64 {
        if stream_id < self.streams.len() {
            self.streams[stream_id].raw_size
        } else {
            0
        }
    }

    /// Get number of streams
    pub fn get_num_streams(&self) -> usize {
        self.streams.len()
    }

    /// Get number of parts in a stream
    pub fn get_num_parts(&self, stream_id: usize) -> usize {
        if stream_id < self.streams.len() {
            self.streams[stream_id].parts.len()
        } else {
            0
        }
    }

    /// Get the next part from a stream (sequential reading)
    pub fn get_part(&mut self, stream_id: usize) -> Result<Option<(Vec<u8>, u64)>> {
        if stream_id >= self.streams.len() {
            anyhow::bail!("Invalid stream ID: {}", stream_id);
        }

        let stream = &mut self.streams[stream_id];
        if stream.cur_id >= stream.parts.len() {
            return Ok(None); // No more parts
        }

        let part = stream.parts[stream.cur_id].clone();
        stream.cur_id += 1;

        self.read_part_data(&part)
    }

    /// Get a specific part by ID from a stream (random access)
    pub fn get_part_by_id(&mut self, stream_id: usize, part_id: usize) -> Result<(Vec<u8>, u64)> {
        if stream_id >= self.streams.len() {
            anyhow::bail!("Invalid stream ID: {}", stream_id);
        }

        let stream = &self.streams[stream_id];
        if part_id >= stream.parts.len() {
            anyhow::bail!("Invalid part ID: {}", part_id);
        }

        let part = stream.parts[part_id].clone();
        self.read_part_data(&part)
            .map(|opt| opt.expect("Part should exist"))
    }

    /// Read part data from file
    fn read_part_data(&mut self, part: &Part) -> Result<Option<(Vec<u8>, u64)>> {
        if part.size == 0 {
            return Ok(Some((Vec::new(), 0)));
        }

        let reader = self
            .reader
            .as_mut()
            .context("Archive not open for reading")?;

        // Seek to part offset
        reader.seek(SeekFrom::Start(part.offset))?;

        // Read metadata
        let (metadata, _) = read_varint(reader)?;

        // Read data (part.size is the data size, not including metadata)
        let mut data = vec![0u8; part.size as usize];
        reader.read_exact(&mut data)?;

        Ok(Some((data, metadata)))
    }

    /// Serialize footer to file (write mode)
    fn serialize(&mut self) -> Result<()> {
        let writer = self
            .writer
            .as_mut()
            .context("Archive not open for writing")?;

        let mut footer = Vec::new();

        // Write number of streams
        write_varint(&mut footer, self.streams.len() as u64)?;

        // Write each stream's metadata
        for stream in &mut self.streams {
            // Stream name (null-terminated string)
            footer.extend_from_slice(stream.stream_name.as_bytes());
            footer.push(0);

            // Number of parts
            write_varint(&mut footer, stream.parts.len() as u64)?;

            // Raw size
            write_varint(&mut footer, stream.raw_size)?;

            // Part offsets and sizes
            for part in &stream.parts {
                write_varint(&mut footer, part.offset)?;
                write_varint(&mut footer, part.size)?;
            }

            // Update packed size to include footer overhead
            stream.packed_size += footer.len() as u64;
        }

        // Write footer
        writer.write_all(&footer)?;

        // Write footer size as fixed 8-byte value
        let footer_size = footer.len() as u64;
        writer.write_all(&footer_size.to_le_bytes())?;

        writer.flush()?;
        Ok(())
    }

    /// Deserialize footer from file (read mode)
    fn deserialize(&mut self) -> Result<()> {
        let file = self.file.as_mut().context("Archive not open")?;

        // Get file size
        let file_size = file.metadata()?.len();

        // Read footer size (last 8 bytes)
        file.seek(SeekFrom::End(-8))?;
        let mut footer_size_bytes = [0u8; 8];
        file.read_exact(&mut footer_size_bytes)?;
        let footer_size = u64::from_le_bytes(footer_size_bytes);

        // Seek to start of footer
        file.seek(SeekFrom::Start(file_size - 8 - footer_size))?;

        // Read footer into buffer
        let mut footer = vec![0u8; footer_size as usize];
        file.read_exact(&mut footer)?;

        // Parse footer
        let mut cursor = std::io::Cursor::new(&footer);

        // Read number of streams
        let (num_streams, _) = read_varint(&mut cursor)?;

        self.streams.clear();
        self.stream_map.clear();

        for i in 0..num_streams {
            // Read stream name (null-terminated)
            let mut stream_name = String::new();
            loop {
                let mut byte = [0u8; 1];
                cursor.read_exact(&mut byte)?;
                if byte[0] == 0 {
                    break;
                }
                stream_name.push(byte[0] as char);
            }

            // Read number of parts
            let (num_parts, _) = read_varint(&mut cursor)?;

            // Read raw size
            let (raw_size, _) = read_varint(&mut cursor)?;

            // Create stream
            let mut stream = Stream::new(stream_name.clone());
            stream.raw_size = raw_size;

            // Read parts
            for _ in 0..num_parts {
                let (offset, _) = read_varint(&mut cursor)?;
                let (size, _) = read_varint(&mut cursor)?;
                stream.parts.push(Part::new(offset, size));
            }

            // Reset cur_id for reading
            stream.cur_id = 0;

            self.streams.push(stream);
            self.stream_map.insert(stream_name, i as usize);
        }

        // Seek back to beginning for reading parts
        file.seek(SeekFrom::Start(0))?;
        if let Some(ref mut reader) = self.reader {
            reader.seek(SeekFrom::Start(0))?;
        }

        Ok(())
    }
}

impl Drop for Archive {
    fn drop(&mut self) {
        let _ = self.close();
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use std::fs;

    #[test]
    fn test_archive_write_read() {
        let path = "test_archive.agc";

        // Write
        {
            let mut archive = Archive::new_writer();
            archive.open(path).unwrap();

            let stream_id = archive.register_stream("test_stream");
            archive.add_part(stream_id, b"Hello", 42).unwrap();
            archive.add_part(stream_id, b"World", 99).unwrap();
            archive.set_raw_size(stream_id, 100);

            archive.close().unwrap();
        }

        // Read
        {
            let mut archive = Archive::new_reader();
            archive.open(path).unwrap();

            let stream_id = archive.get_stream_id("test_stream").unwrap();
            assert_eq!(archive.get_num_parts(stream_id), 2);
            assert_eq!(archive.get_raw_size(stream_id), 100);

            let (data1, meta1) = archive.get_part(stream_id).unwrap().unwrap();
            assert_eq!(data1, b"Hello");
            assert_eq!(meta1, 42);

            let (data2, meta2) = archive.get_part(stream_id).unwrap().unwrap();
            assert_eq!(data2, b"World");
            assert_eq!(meta2, 99);

            assert!(archive.get_part(stream_id).unwrap().is_none());
        }

        fs::remove_file(path).unwrap();
    }

    #[test]
    fn test_multiple_streams() {
        let path = "test_multi_stream.agc";

        // Write
        {
            let mut archive = Archive::new_writer();
            archive.open(path).unwrap();

            let stream1 = archive.register_stream("stream1");
            let stream2 = archive.register_stream("stream2");

            archive.add_part(stream1, b"Data1", 1).unwrap();
            archive.add_part(stream2, b"Data2", 2).unwrap();
            archive.add_part(stream1, b"Data3", 3).unwrap();

            archive.close().unwrap();
        }

        // Read
        {
            let mut archive = Archive::new_reader();
            archive.open(path).unwrap();

            let stream1 = archive.get_stream_id("stream1").unwrap();
            let stream2 = archive.get_stream_id("stream2").unwrap();

            assert_eq!(archive.get_num_parts(stream1), 2);
            assert_eq!(archive.get_num_parts(stream2), 1);

            let (data, meta) = archive.get_part_by_id(stream1, 0).unwrap();
            assert_eq!(data, b"Data1");
            assert_eq!(meta, 1);

            let (data, meta) = archive.get_part_by_id(stream2, 0).unwrap();
            assert_eq!(data, b"Data2");
            assert_eq!(meta, 2);
        }

        fs::remove_file(path).unwrap();
    }
}
