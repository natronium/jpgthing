use std::{
    collections::HashMap,
    iter::{once, repeat},
    usize,
};

use bitvec::prelude::*;

use crate::jpeg::HuffmanTable;

//Implementation of the huffman decoder as specified by t.81
//Figure C.1
fn generate_size_table(huff_data: &HuffmanTable) -> (Vec<u8>, usize) {
    let huffsize: Vec<u8> = huff_data
        .ls
        .iter()
        .enumerate()
        .flat_map(|(idx, count)| repeat((idx + 1) as u8).take(*count as usize))
        .chain(once(0))
        .collect();
    let last_k = huffsize.len() - 1;
    (huffsize, last_k)
}

//figure C.2
//TODO: make this less imperative? more idiomatic? surely there's a mathier way to represent this
fn generate_code_table(huffsize: &[u8]) -> Vec<u16> {
    let mut k = 0;
    let mut code = 0;
    let mut si = huffsize[0];
    let mut huffcode = Vec::new();

    loop {
        loop {
            //huffcode[k] = code;
            huffcode.push(code);
            code += 1;
            k += 1;

            if huffsize[k] != si {
                break;
            }
        }

        if huffsize[k] == 0 {
            return huffcode;
        }

        while huffsize[k] != si {
            code <<= 1;
            si += 1;
        }
    }
}

//figure C.3
// huffval = symbol_length_assignment, but flat (e.g. sla.into_iter().flatten().collect)
pub fn order_codes(
    huffsize: &[u8],
    huffcode: &[u8],
    huffval: &[u8],
) -> (HashMap<u8, u8>, HashMap<u8, u8>) {
    let mut ehufsi = HashMap::new();
    let mut ehufco = HashMap::new();

    for (idx, value) in huffval.iter().enumerate() {
        ehufco.insert(*value, huffcode[idx]);
        ehufsi.insert(*value, huffsize[idx]);
    }
    (ehufsi, ehufco)
}

//figure F.15
fn decoder_tables(bits: &[u8; 16], huffcode: &[u16]) -> ([u16; 16], [i32; 16], [usize; 16]) {
    let mut j: usize = 0;
    let mut mincode = [0u16; 16];
    let mut maxcode = [0i32; 16];
    let mut valptr = [0usize; 16];
    for i in 0..16 {
        if bits[i] == 0 {
            maxcode[i] = -1;
            continue;
        }

        valptr[i] = j;
        mincode[i] = huffcode[j].into();
        j += bits[i] as usize - 1; //TODO: why can't i just .into here?
        maxcode[i] = huffcode[j].into();
        j += 1;
    }
    (mincode, maxcode, valptr)
}

pub struct BitReader {
    base: BitVec<u8, Msb0>,
    offset: usize,
}

impl BitReader {
    pub fn new(backing_store: BitVec<u8, Msb0>) -> Self {
        Self {
            base: backing_store,
            offset: 0,
        }
    }

    fn peek_nbits(&self, n: usize) -> &BitSlice<u8, Msb0> {
        &self.base[self.offset..(n + self.offset)]
    }

    fn read_bits(&mut self, num_bits: usize) -> &BitSlice<u8, Msb0> {
        let bits = &self.base[self.offset..(num_bits + self.offset)];
        self.offset += num_bits;
        bits
    }

    //Figure F.16
    fn decode(&mut self, huf_info: &HuffmanInfo) -> u8 {
        let mut i = 1;
        //TODO: i hate this. surely there is a better way to do this
        while self.peek_nbits(i).load_be::<u16>() as i32 > huf_info.maxcode[i - 1] {
            i += 1;
        }

        let code_slice = self.read_bits(i);

        let code: usize = code_slice.load_be();

        let j: usize = huf_info.valptr[i - 1] + code - huf_info.mincode[i - 1] as usize;

        huf_info.hufval[j]
    }

    pub fn is_empty(&self) -> bool {
        //if there are 7 or fewer bits left
        // AND all the remaining bits are 1
        // we're done.

        let remaining_bit_count = self.base.len() - self.offset;
        if remaining_bit_count < 8 {
            let remaining_bits = self.peek_nbits(remaining_bit_count);
            remaining_bits.count_ones() == remaining_bit_count
        } else {
            false
        }
    }
}

pub struct HuffmanInfo {
    hufval: Vec<u8>,
    mincode: [u16; 16],
    maxcode: [i32; 16],
    valptr: [usize; 16],
}

impl HuffmanInfo {
    pub fn new(huff_data: &HuffmanTable) -> Self {
        let (huffsize, _last_k) = generate_size_table(huff_data);
        let huffcode = generate_code_table(&huffsize);
        let (mincode, maxcode, valptr) = decoder_tables(&huff_data.ls, &huffcode);
        let hufval: Vec<u8> = huff_data
            .symbol_length_assignment
            .clone()
            .into_iter()
            .flatten()
            .collect();
        HuffmanInfo {
            hufval,
            mincode,
            maxcode,
            valptr,
        }
    }
}

//figure F.12
fn extend(v: u16, t: u8) -> i16 {
    let v_t = 1 << (t - 1);
    if v < v_t {
        (v as i16) + ((-1 << t) + 1)
    } else {
        v as i16
    }
}

/// returns all 63 ac coefficients in the scan in zigzag order. the slice starts at 0!
// Figure F.13
pub fn read_ac<'a, 'b>(stream: &mut BitReader, huff_info: &HuffmanInfo) -> [i16; 63] {
    let mut k: usize = 0;
    let mut ac_coeffs = [0; 63];
    loop {
        //rs = 0bRRRR_SSSS
        let rs = stream.decode(huff_info);

        let zero_prerun = rs >> 4; //RRRR aka R
        let amplitude = rs & 0b0000_1111; //SSSS

        if amplitude == 0 {
            if zero_prerun == 15 {
                //pre-run of 15 zeros + current zero = 16
                k += 16;
                continue;
            } else {
                //EOB (but also funkiness with nonzero, non-fifteen values for R?)
                break;
            }
        } else {
            k += usize::from(zero_prerun);
            ac_coeffs[k] = extend(stream.read_bits(amplitude.into()).load_be(), amplitude)
        }

        if k == 62 {
            break;
        } else {
            k += 1;
        }
    }

    ac_coeffs
}

// Section F.2.2.1
pub fn read_dc_diff<'a>(stream: &mut BitReader, huff_info: &HuffmanInfo) -> i16 {
    let amplitude = stream.decode(huff_info);
    if amplitude == 0 {
        //no bits to read. value is just 0
        return 0;
    }

    let bits: u16 = stream.read_bits(amplitude.into()).load_be(); //TODO: am i off-by-one?

    extend(bits, amplitude)
}
