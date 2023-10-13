use std::fs::File;
use binrw::BinRead;
use jpeg::JIF;

pub mod jpeg;
mod helpers;

fn main() {
    let filepath = "/tmp/flowers.jpg";
    let mut file = File::open(filepath).unwrap();
    let jif_data = JIF::read_be(&mut file);
    println!("{:#?}", jif_data);
}
