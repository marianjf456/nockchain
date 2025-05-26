//! jets/stark_jets.rs
//!
//! Native Rust jet for forward NTT (`bp-ntt`).
use nockvm::interpreter::Context;
use nockvm::jets::util::slot;
use nockvm::jets::JetErr;
use nockvm::noun::{Noun, CellMemory, NounAllocator};
use crate::form::math::base::PRIME;

/// STARK base field modulus: 2^64 − 2^32 + 1
const MOD: u64 = PRIME;

/// Fast modular exponentiation.
#[inline(always)]
fn mod_exp(mut base: u64, mut exp: u64) -> u64 {
    let mut acc = 1u64;
    base %= MOD;
    while exp != 0 {
        if (exp & 1) != 0 {
            acc = ((acc as u128 * base as u128) % (MOD as u128)) as u64;
        }
        base = ((base as u128 * base as u128) % (MOD as u128)) as u64;
        exp >>= 1;
    }
    acc
}

/// In-place Cooley–Tuk radix-2 NTT.
#[inline(always)]
fn radix2_ntt_inplace(a: &mut [u64], root: u64) {
    let n = a.len();
    assert!(n.is_power_of_two(), "NTT size must be a power of two");

    // Bit-reversal permutation
    {
        let mut j = 0;
        for i in 1..n {
            let mut bit = n >> 1;
            while (j & bit) != 0 {
                j ^= bit;
                bit >>= 1;
            }
            j |= bit;
            if i < j {
                a.swap(i, j);
            }
        }
    }

    // Cooley–Tukey butterfly loops
    let mut len = 2;
    while len <= n {
        let wlen = mod_exp(root, (n / len) as u64);
        for i in (0..n).step_by(len) {
            let mut w = 1u64;
            for j in 0..(len / 2) {
                let u = a[i + j];
                let v = ((a[i + j + len / 2] as u128 * w as u128) % (MOD as u128)) as u64;
                a[i + j] = (u + v) % MOD;
                a[i + j + len / 2] = (u + MOD - v) % MOD;
                w = ((w as u128 * wlen as u128) % (MOD as u128)) as u64;
            }
        }
        len <<= 1;
    }
}

/// Convert a Hoon list of atoms into a `Vec<u64>`.
fn list_to_vec_atoms(mut list: Noun) -> Result<Vec<u64>, JetErr> {
    let mut out = Vec::new();
    while let Ok(cell) = list.as_cell() {
        let atom = cell.head().as_direct()?;
        out.push(atom.data());
        list = cell.tail();
    }
    Ok(out)
}

/// The hot-entry for `%stark-ntt` (forward NTT).
/// Expects a 2-cell `[ bpoly , root ]`.
pub fn stark_ntt_jet(context: &mut Context, subject: Noun) -> Result<Noun, JetErr> {
    // Extract the “door” and its arguments
    let door      = slot(subject, 7)?;
    let bp_noun   = slot(door,   6)?;
    let root_noun = slot(door,   7)?;

    // Deserialize bpoly into Vec<u64>
    let mut data = list_to_vec_atoms(bp_noun)?;
    let root = root_noun.as_direct()?.data();

    // Perform the in-place NTT
    radix2_ntt_inplace(&mut data, root);

    // Pack Vec<u64> back into a Hoon list: empty list = raw 0
    let mut list_n = unsafe { Noun::from_raw(0) };
    for &x in data.iter().rev() {
        let atom_n = unsafe { Noun::from_raw(x) };
        // allocate a new cons-cell
        let cell_ptr: *mut CellMemory = unsafe { context.stack.alloc_cell() };
        unsafe {
            (*cell_ptr).head = atom_n.into();
            (*cell_ptr).tail = list_n.into();
        }
        list_n = unsafe { Noun::from_raw(cell_ptr as u64) };
    }

    Ok(list_n)
}