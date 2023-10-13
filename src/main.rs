use binrw::{binread, parser, BinRead, BinResult, Endian};
use std::fs::File;
use std::io::{Read, Seek, SeekFrom};
use std::iter::from_fn;

fn main() {
    let filepath = "/tmp/flowers.jpg";
    let mut file = File::open(filepath).unwrap();
    let jif_data = JIF::read_be(&mut file);
    println!("{:#?}", jif_data);
}


//stolen and readapted from binrw::helpers::{until_eof, until_eof_with}
//TODO: understand how this works and why i couldn't write a non-clojurey generic version
pub fn until_invalid<Reader, T, Arg, Ret>(
    reader: &mut Reader,
    endian: Endian,
    args: Arg,
) -> BinResult<Ret>
where
    T: for<'a> BinRead<Args<'a> = Arg>,
    Reader: Read + Seek,
    Arg: Clone,
    Ret: FromIterator<T>,
{
    until_invalid_with(T::read_options)(reader, endian, args)
}

pub fn until_invalid_with<Reader, T, Arg, ReadFn, Ret>(
    read: ReadFn,
) -> impl Fn(&mut Reader, Endian, Arg) -> BinResult<Ret>
where
    Reader: Read + Seek,
    Arg: Clone,
    ReadFn: Fn(&mut Reader, Endian, Arg) -> BinResult<T>,
    Ret: FromIterator<T>,
{
    move |reader, endian, args| {
        from_fn(|| match read(reader, endian, args.clone()) {
            ok @ Ok(_) => Some(ok),
            Err(_) => None,
        })
        .fuse()
        .collect()
    }
}

#[binread]
#[derive(Debug)]
#[br(magic = 0xFFD8u16)]
struct JIF {
    jpeg_frame: Frame,
    #[br(magic = 0xFFD9u16, dbg, temp)]
    _eoi: (),
}

#[binread]
#[derive(Debug)]
struct Frame {
    #[br(parse_with = until_invalid)]
    misc_tables: Vec<MiscSegment>,
    header: FrameHeader,
    scan_1: Scan,
    #[br(try)]
    dnl: Option<DNL>,
    #[br(parse_with = until_invalid)]
    remaining_scans: Vec<Scan>,
}

#[binread]
#[derive(Debug)]
enum MiscSegment {
    #[br(magic = 0xFFDBu16)]
    QuantizationTable {
        len: u16,
        #[br(count = len - 2)]
        body: Vec<u8>,
    },
    #[br(magic = 0xFFC4u16)]
    HuffmanTable {
        len: u16,
        #[br(count = len - 2)]
        body: Vec<u8>,
    },
    #[br(magic = 0xFFCCu16)]
    ArithmeticConditioningTable {
        len: u16,
        #[br(count = len - 2)]
        body: Vec<u8>,
    },
    #[br(magic = 0xFFDDu16)]
    RestartIntervalDefinition {
        len: u16,
        #[br(count = len - 2)]
        body: Vec<u8>,
    },
    #[br(magic = 0xFFFEu16)]
    Comment {
        len: u16,
        #[br(count = len - 2)]
        body: Vec<u8>,
    },
    //FFE0-FFEF
    #[br(magic = 0xFFu8)]
    #[br(assert( n >= 0xE0 && n <= 0xEF))]
    AppData {
        n: u8,
        len: u16,
        #[br(count = len - 2)]
        body: Vec<u8>,
    },
}

#[binread]
#[derive(Debug)]
#[br(magic = 0xFFDCu16)]
struct DNL {
    len: u16, //is always = 4. 2 bytes for this field, 2 for the next
    nl: u16,
}

#[binread]
#[derive(Debug)]
#[br(assert(n >= 0xC0 && n <= 0xCF && n != 0xC4 && n != 0xC8 && n!= 0xCC))]
#[br(assert(_marker == 0xFF))]
struct SOFMarker {
    _marker: u8,
    n: u8,
}

#[binread]
#[derive(Debug)]
struct FrameHeader {
    sof_marker: SOFMarker,
    lf: u16,
    p: u8,
    y: u16,
    x: u16,
    nf: u8,
    #[br(count = nf)]
    component_parameters: Vec<FrameComponentParameterSet>,
}

#[binread]
#[derive(Debug)]
struct FrameComponentParameterSet {
    c: u8,
    #[br(temp)]
    _raw_h_v: u8,
    #[br(calc = (_raw_h_v & 0b1111_0000) >> 4)]
    h: u8, //u4
    #[br(calc = _raw_h_v & 0b0000_1111)]
    v: u8, //u4
    tq: u8,
}

#[binread]
#[derive(Debug)]
#[br(assert(_raw_n >= 0xD0 && _raw_n <= 0xD7))]
#[br(assert(_marker == 0xFF))]
struct RSTMarker {
    _marker: u8,
    _raw_n: u8,
    #[br(calc = _raw_n & 0b0000_1111)]
    n: u8,
}

//TODO i'm pretty sure there's a more idiomatic way to do this
//TODO i'm pretty sure there's also a faster way to do this
#[parser(reader)]
fn parse_ecs() -> BinResult<Vec<u8>> {
    let mut ret = Vec::new();
    let mut byte_buf = [0; 1];
    loop {
        reader.read_exact(&mut byte_buf)?;
        let byte = byte_buf[0];
        if byte != 0xff {
            ret.push(byte);
            continue;
        } else if byte == 0xff {
            reader.read_exact(&mut byte_buf)?;
            let next_byte = byte_buf[0];
            if next_byte == 0x00 {
                ret.push(byte);
            } else if next_byte != 0x00 {
                reader.seek(SeekFrom::Current(-2))?;
                break;
            }
        }
    }

    return Ok(ret);
}

#[binread]
#[derive(Debug)]
struct Scan {
    #[br(parse_with = until_invalid)]
    misc_tables: Vec<MiscSegment>,
    header: ScanHeader,
    #[br(parse_with = parse_ecs)]
    raw_ecs_data: Vec<u8>,
}

#[binread]
#[derive(Debug)]
#[br(magic = 0xFFDAu16)]
struct ScanHeader {
    len: u16,
    ns: u8,
    #[br(count = ns)]
    component_parameters: Vec<ScanComponentParameterSet>,
    ss: u8,
    se: u8,
    _raw_ah_al: u8,
    #[br(calc = (_raw_ah_al &0b1111_0000) >> 4)]
    ah: u8, //u4
    #[br(calc = (_raw_ah_al &0b0000_1111) >> 0)]
    al: u8, //u4
}

#[binread]
#[derive(Debug)]
struct ScanComponentParameterSet {
    cs: u8,
    #[br(temp)]
    _raw_td_ta: u8,
    #[br(calc = (_raw_td_ta &0b1111_0000) >> 4)]
    td: u8, //u4
    #[br(calc = (_raw_td_ta &0b0000_1111) >> 0)]
    ta: u8, //u4
}

/*
A Jpeg File is:
SOI, Frame, EOI

A Frame is:
[Tables/Misc], Frame header, Scan_1, [DNL], [Scan_2, ..., Scan_last]

A Frame Header is:
SOF_n, Lf:u16, P:u8, Y:u16, X:u16, Nf:u8, Component-Specification Parameters

Component-Specification Parameters:
C_1:u8, H_1:u4, V_1:u4, Tq_1:u8, ... Nf times

A Scan is:
[Tables/Misc], Scan Header [ECS_0, RST_0, ..., ECS_last-1, RST_last-1], ECS_last

An Etropy-coded segment (ECS) is:
<MCU_1>, <MCU_2>, ..., <MCU_Ri> // when there's a reset interval
<MCU_n>, <MCU_n+1>, ..., <MCU_last> // for ECS_last
*/

/*
JIF:
0xff+ :marker-byte: [Marker contents]
markers:
    Start of Frame:

*/
