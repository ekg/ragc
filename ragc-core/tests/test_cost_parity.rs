#![cfg(feature = "ffi_cost")]

use ragc_core::{ragc_ffi, LZDiff};

fn cumulative_left(v: &[u32]) -> Vec<u32> {
    let mut out = Vec::with_capacity(v.len());
    let mut acc = 0u32;
    for &c in v { acc = acc.saturating_add(c); out.push(acc); }
    out
}

fn cumulative_right_suffix(v: &[u32]) -> Vec<u32> {
    let mut out = v.to_vec();
    let mut acc = 0u32;
    for c in out.iter_mut().rev() { acc = acc.saturating_add(*c); *c = acc; }
    out
}

#[test]
fn ffi_cost_all_literals() {
    // Text shorter than key_len -> all literals cost 1
    let mml = 20u32;
    let reference: Vec<u8> = (0..100).map(|i| (i % 4) as u8).collect();
    let text: Vec<u8> = vec![1u8; 10];

    let v_prefix = ragc_ffi::cost_vector(true, &reference, &text, mml);
    let v_suffix = ragc_ffi::cost_vector(false, &reference, &text, mml);

    assert_eq!(v_prefix.len(), text.len());
    assert_eq!(v_suffix.len(), text.len());
    assert!(v_prefix.iter().all(|&c| c == 1));
    assert!(v_suffix.iter().all(|&c| c == 1));
}

// Note: A stricter shape test will be added using golden vectors generated
// via scripts/cost_verifier.cpp to avoid brittle assumptions here.

#[test]
fn ffi_vs_rust_vectors_document_current_gap() {
    // This test documents the current gap: Rust cost model may differ from C++.
    // It ensures both paths run and produce vectors of the same length on non-trivial input.
    let mml = 20u32;
    let reference: Vec<u8> = (0..256).map(|i| (i % 4) as u8).collect();
    let text: Vec<u8> = (0..180).map(|i| ((i * 3 + 1) % 4) as u8).collect();

    // FFI (C++)
    let v_cpp_pref = ragc_ffi::cost_vector(true, &reference, &text, mml);
    let v_cpp_suf = ragc_ffi::cost_vector(false, &reference, &text, mml);

    // Pure Rust
    let mut lz = LZDiff::new(mml);
    lz.prepare(&reference);
    let v_rs_pref = lz.get_coding_cost_vector(&text, true);
    let v_rs_suf = lz.get_coding_cost_vector(&text, false);

    assert_eq!(v_cpp_pref.len(), text.len());
    assert_eq!(v_cpp_suf.len(), text.len());
    assert_eq!(v_rs_pref.len(), text.len());
    assert_eq!(v_rs_suf.len(), text.len());

    // Do not assert equality yet; this codifies the current state and unblocks using FFI as ground truth.
}
