#!/bin/bash
# Test C++ ↔ Rust compatibility

set -e

echo "=== Creating test FASTA ==="
cat > /tmp/test_compat.fasta << 'EOF'
>chr1
ACGTACGTACGTACGTACGTACGTACGTACGTACGTACGT
>chr2
TGCATGCATGCATGCATGCATGCATGCATGCATGCATGCA
EOF

echo ""
echo "=== Cleaning up old archives ==="
rm -f /tmp/cpp_compat.agc /tmp/rust_compat.agc

echo ""
echo "=== Creating archive with C++ ==="
/home/erik/agc/bin/agc create -o /tmp/cpp_compat.agc /tmp/test_compat.fasta

echo ""
echo "=== Creating archive with Rust ==="
cd /home/erik/agc/rust-agc
cargo run --release --bin agc -- create --output /tmp/rust_compat.agc /tmp/test_compat.fasta

echo ""
echo "=== Comparing archive sizes ==="
ls -lh /tmp/cpp_compat.agc /tmp/rust_compat.agc

echo ""
echo "=== Hex dump comparison (first 512 bytes) ==="
echo "--- C++ archive ---"
xxd -l 512 /tmp/cpp_compat.agc > /tmp/cpp_hex.txt
head -20 /tmp/cpp_hex.txt

echo ""
echo "--- Rust archive ---"
xxd -l 512 /tmp/rust_compat.agc > /tmp/rust_hex.txt
head -20 /tmp/rust_hex.txt

echo ""
echo "=== Testing Rust reading C++ archive ==="
cd /home/erik/agc/rust-agc
cargo run --release --bin agc -- listset /tmp/cpp_compat.agc || echo "FAILED: Rust cannot read C++ archive"

echo ""
echo "=== Testing C++ reading Rust archive ==="
/home/erik/agc/bin/agc listset /tmp/rust_compat.agc || echo "FAILED: C++ cannot read Rust archive"
