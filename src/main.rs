use binrw::BinRead;
use jpeg::{HuffmanTable, JIF};
use std::{
    collections::HashMap,
    f64::consts::{FRAC_1_SQRT_2, PI},
    fs::File,
    iter::{once, repeat},
};

use crate::jpeg::MiscSegment;

mod helpers;
pub mod jpeg;

fn main() {
    //let filepath = "/tmp/testimages/peppers-50.jpg";
    let filepath = "/tmp/flowers-customdht.jpg";
    let mut file = File::open(filepath).unwrap();
    let jif_data = JIF::read_be(&mut file).unwrap();

    let dhts = jif_data
        .jpeg_frame
        .scan_1
        .misc_tables
        .iter()
        .filter(|e| matches!(e, MiscSegment::DefineHuffmanTable { len: _, tables: _ }));

    for dht in dhts {
        if let MiscSegment::DefineHuffmanTable { len: _, tables: ts } = dht {
            for t in ts {
                println!("ls: {:?}", t.ls);
                let (huffsize, _last_k) = also_generate_size_table(t);
                let huffcode = generate_code_table(&huffsize);
                let huffval = t
                    .symbol_length_assignment
                    .clone()
                    .into_iter()
                    .flatten()
                    .collect::<Vec<u8>>();
                let gen = order_codes(&huffsize, &huffcode, &huffval);
                //let also_gen = also_generate_code_table(&huffsize.0);
                println!("gen:      {:?}", gen);
                // println!("also_gen: {:?}", also_gen);
                //println!("eq:       {}", gen == also_gen)
                println!("---------")
            }
        }
        println!("==========")
    }
}

//Implementation of the huffman decoder as specified by t.81
//Figure C.1
#[allow(dead_code)]
fn generate_size_table(huff_data: &HuffmanTable) -> (Vec<u8>, usize) {
    let mut k = 0;
    let mut i = 1;
    let mut j = 1;
    let mut huffsize = Vec::new();
    let bits = &huff_data.ls;

    while i <= 16 {
        while j <= bits[i - 1] {
            //t81 1=indexes, but rust 0-indexes
            //huffsize[k] = i;
            //since k is strictly monotonically increasing, this is equivalent, and also not an error
            huffsize.push(i as u8);
            k += 1;
            j += 1
        }
        i += 1;
        j = 1;
    }
    //huffsize[k] = 0;
    huffsize.push(0);

    (huffsize, k)
}

//now lets do it with iterators
fn also_generate_size_table(huff_data: &HuffmanTable) -> (Vec<u8>, usize) {
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
fn generate_code_table(huffsize: &[u8]) -> Vec<u8> {
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
fn order_codes(
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

//fdct and idct (hopefully) implemented ~verbatim from A.3.3
//TODO: test these?
#[allow(unused)]
fn fdct(samples: [[u8; 8]; 8]) -> [[f64; 8]; 8] {
    let mut dct = [[0.0; 8]; 8];
    for u in 0u8..7 {
        for v in 0u8..7 {
            let c_u = if u == 0 { FRAC_1_SQRT_2 } else { 1.0 };
            let c_v = if v == 0 { FRAC_1_SQRT_2 } else { 1.0 };
            let mut sum = 0.0;
            for x in 0u8..7 {
                for y in 0u8..7 {
                    let s_yx = f64::from(samples[y as usize][x as usize]);
                    let horiz = ((((2.0 * f64::from(x)) + 1.0) * f64::from(u) * PI) / 16.0).cos();
                    let vert = ((((2.0 * f64::from(y)) + 1.0) * f64::from(v) * PI) / 16.0).cos();
                    sum += s_yx * horiz * vert;
                }
            }
            dct[v as usize][u as usize] = 0.25 * c_u * c_v * sum;
        }
    }

    dct
}

#[allow(unused)]
fn idct(samples: [[f64; 8]; 8]) -> [[u8; 8]; 8] {
    let mut undct = [[0; 8]; 8];
    for x in 0u8..7 {
        for y in 0u8..7 {
            let mut sum = 0.0;
            for u in 0u8..7 {
                for v in 0u8..7 {
                    let c_u = if u == 0 { FRAC_1_SQRT_2 } else { 1.0 };
                    let c_v = if v == 0 { FRAC_1_SQRT_2 } else { 1.0 };
                    let s_vu = samples[v as usize][u as usize];
                    let horiz = ((((2.0 * f64::from(x)) + 1.0) * f64::from(u) * PI) / 16.0).cos();
                    let vert = ((((2.0 * f64::from(y)) + 1.0) * f64::from(v) * PI) / 16.0).cos();
                    sum += c_u * c_v * s_vu * horiz * vert
                }
            }
            //TODO: make lossy rounding explicit? this seems sketchy as-is
            undct[y as usize][x as usize] = (0.25 * sum) as u8;
        }
    }

    undct
}
