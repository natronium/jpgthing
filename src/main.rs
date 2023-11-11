use binrw::BinRead;
use jpeg::JIF;
use std::{fs::File, io::Write};

pub mod decoder;
mod helpers;
pub mod jpeg;
pub mod jpeg_huff;

fn main() {
    //let filepath = "/tmp/testimages/peppers-50.jpg";
    let filepath = "/home/na/testjpegs/flowers.jpg";
    let mut file = File::open(filepath).unwrap();
    let jif_data = JIF::read_be(&mut file).unwrap();
    //println!("{:#?}", jif_data);

    let frame = jif_data.jpeg_frame;

    let components = decoder::read_frame(&frame).unwrap();

    for component in components {
        let mut f = File::create(format!("/tmp/component_{}.pgm", component.id)).unwrap();
        write!(f, "P5 {} {} 255\n", component.underlying_width, component.underlying_height).unwrap();
        f.write(
            component
                .body
                .iter()
                .map(|n| *n as u8)
                .collect::<Vec<_>>()
                .as_slice(),
        )
        .unwrap();
    }
}
