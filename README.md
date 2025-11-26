# ragc

Rust reimplementation of [AGC](https://github.com/refresh-bio/agc) (Assembled Genomes Compressor). Produces bit-compatible archives readable by C++ AGC and vice versa.

## Installation

```bash
git clone https://github.com/ekg/ragc.git
cd ragc
cargo build --release
```

## Usage

```bash
# Create archive from FASTA
ragc create --output genomes.agc sample.fasta

# Create from multiple files
ragc create --output genomes.agc *.fasta

# Extract samples
ragc getset genomes.agc sample_name > output.fasta

# List samples/contigs
ragc listset genomes.agc
ragc listctg genomes.agc
```

## Library Usage

Add to `Cargo.toml`:
```toml
ragc-reader = { git = "https://github.com/ekg/ragc.git" }
```

```rust
use ragc_reader::{Decompressor, DecompressorConfig};

let mut dec = Decompressor::open("data.agc", DecompressorConfig::default())?;
for sample in dec.list_samples() {
    let contigs = dec.get_sample(&sample)?;
    for (name, sequence) in contigs {
        println!(">{}\n{}", name, String::from_utf8_lossy(&sequence));
    }
}
dec.close()?;
```

## Project Structure

```
ragc-core/     # Shared types: archive I/O, collection metadata, varint encoding
ragc-reader/   # Decompressor (for reading AGC files)
ragc-writer/   # Compressor, k-mer extraction, LZ diff (for creating AGC files)
ragc-cli/      # CLI: create, getset, listset, listctg commands
```

## Status

**Working:**
- Archive creation and reading (C++ AGC compatible)
- Multi-sample, multi-contig support
- Streaming compression with bounded memory
- PanSN format support

**Limitations:**
- Compression is suboptimal (archives larger than C++ AGC)
- Some CLI commands not implemented (getcol, etc.)

## Development

```bash
cargo build --release    # Build
cargo test               # Run tests
cargo clippy             # Lint
```

## License

MIT - same as the original C++ AGC project.

## Acknowledgments

Reimplementation of [AGC](https://github.com/refresh-bio/agc) by Sebastian Deorowicz and Adam Gudyś ([REFRESH Bioinformatics Group](https://github.com/refresh-bio)).
