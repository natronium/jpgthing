use std::io::Cursor;

use crate::helpers::{parse_ecs, until_invalid};
use binrw::{args, binread, helpers::args_iter, BinRead, Endian};

/*
A (JIF) Jpeg File is:
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

//TODO: 4bit hi and lo types we use them in 5 places already
//TODO: optional 0xff fill bytes which may preceede a marker and shall be discarded
#[binread]
#[derive(Debug)]
#[br(magic = 0xFFD8u16)]
pub struct JIF {
    pub jpeg_frame: Frame,
    #[br(magic = 0xFFD9u16, temp)]
    pub _eoi: (),
}

#[binread]
#[derive(Debug)]
pub struct Frame {
    #[br(parse_with = until_invalid)]
    pub misc_tables: Vec<MiscSegment>,
    pub header: FrameHeader,
    pub scan_1: Scan,
    #[br(try)]
    pub dnl: Option<DNL>,
    #[br(parse_with = until_invalid)]
    pub remaining_scans: Vec<Scan>,
}

#[binread]
#[derive(Debug)]
pub enum MiscSegment {
    #[br(magic = 0xFFDBu16)]
    DefineQuantizationTable {
        len: u16,
        //TODO: is there a nice way to do this with parse_with? some reasonable way to factor this out into a helper?
        #[br(try_map = |raw_bytes: Vec<u8>| {
            let mut limited_cursor = Cursor::new(raw_bytes);
            until_invalid(&mut limited_cursor, Endian::Big, ())
        }, count = len-2)]
        tables: Vec<QuantizationTable>,
    },
    #[br(magic = 0xFFC4u16)]
    DefineHuffmanTable {
        len: u16,
        #[br(try_map = |raw_bytes: Vec<u8>| {
            let mut limited_cursor = Cursor::new(raw_bytes);
            until_invalid(&mut limited_cursor, Endian::Big, ())
        }, count = len - 2)]
        tables: Vec<HuffmanTable>,
    },
    #[br(magic = 0xFFCCu16)]
    DefineArithmeticConditioningTable {
        len: u16,
        #[br(try_map = |raw_bytes: Vec<u8>| {
            let mut limited_cursor = Cursor::new(raw_bytes);
            until_invalid(&mut limited_cursor, Endian::Big, ())
        }, count = len - 2)]
        tables: Vec<ArithmeticConditioningTable>,
    },
    #[br(magic = 0xFFDDu16)]
    RestartIntervalDefinition { len: u16, n: u16 },
    #[br(magic = 0xFFFEu16)]
    Comment {
        len: u16,
        #[br(count = len - 2)]
        body: Vec<u8>,
    },
    //FFE0-FFEF
    #[br(magic = 0xFFu8)]
    #[br(assert( (0xE0..=0xEF).contains(&n)))]
    AppData {
        n: u8,
        len: u16,
        #[br(count = len - 2)]
        body: Vec<u8>,
    },
}

#[binread]
#[derive(Debug)]
pub struct ArithmeticConditioningTable {
    #[br(temp)]
    _raw_tc_tb: u8,
    #[br(calc = (_raw_tc_tb &0b1111_0000) >> 4)]
    pub tc: u8,
    #[br(calc = _raw_tc_tb &0b0000_1111)]
    pub tb: u8,
    pub cs: u8,
}

#[binread]
#[derive(Debug)]
pub struct QuantizationTable {
    #[br(temp)]
    _raw_pq_tq: u8,
    #[br(calc = (_raw_pq_tq &0b1111_0000) >> 4)]
    pub pq: u8,
    #[br(calc = _raw_pq_tq &0b0000_1111)]
    pub tq: u8,
    #[br(args(pq == 1,))]
    pub qs: QuantizationEntries,
}

#[binread]
#[derive(Debug)]
pub struct HuffmanTable {
    #[br(temp)]
    _raw_tc_th: u8,
    #[br(calc = (_raw_tc_th &0b1111_0000) >> 4)]
    pub tc: u8,
    #[br(calc = _raw_tc_th &0b0000_1111)]
    pub th: u8,
    #[br(count = 16)]
    pub ls: Vec<u8>,
    #[br(parse_with = args_iter(ls.iter().map(|&size| -> <Vec<u8> as BinRead>::Args<'_>  {
        args! {count: size.into()}
    })))]
    pub symbol_length_assignment: Vec<Vec<u8>>,
}

#[binread]
#[derive(Debug)]
#[br(import(hi_precision: bool))]
pub enum QuantizationEntries {
    #[br(assert(!hi_precision))]
    Low(#[br(count = 64)] Vec<u8>),
    #[br(assert(hi_precision))]
    Hi(#[br(count = 64)] Vec<u16>),
}

#[binread]
#[derive(Debug)]
#[br(magic = 0xFFDCu16)]
pub struct DNL {
    pub len: u16, //is always = 4. 2 bytes for this field, 2 for the next
    pub nl: u16,
}

#[binread]
#[derive(Debug)]
#[br(assert((0xC0..=0xCF).contains(&n) && n != 0xC4 && n != 0xC8 && n!= 0xCC))]
#[br(assert(_marker == 0xFF))]
pub struct SOFMarker {
    pub _marker: u8,
    pub n: u8,
}

#[binread]
#[derive(Debug)]
pub struct FrameHeader {
    pub sof_marker: SOFMarker,
    pub lf: u16,
    pub p: u8,
    pub y: u16,
    pub x: u16,
    pub nf: u8,
    #[br(count = nf)]
    pub component_parameters: Vec<FrameComponentParameterSet>,
}

#[binread]
#[derive(Debug)]
pub struct FrameComponentParameterSet {
    pub c: u8,
    #[br(temp)]
    _raw_h_v: u8,
    #[br(calc = (_raw_h_v & 0b1111_0000) >> 4)]
    pub h: u8, //u4
    #[br(calc = _raw_h_v & 0b0000_1111)]
    pub v: u8, //u4
    pub tq: u8,
}

#[binread]
#[derive(Debug)]
pub struct Scan {
    #[br(parse_with = until_invalid)]
    pub misc_tables: Vec<MiscSegment>,
    pub header: ScanHeader,
    #[br(parse_with = until_invalid)]
    pub ecs_segments: Vec<(ECSSegment, RSTSegment)>,
    pub ecs_last: ECSSegment,
}

#[binread]
#[derive(Debug)]
pub struct ECSSegment {
    #[br(parse_with = parse_ecs)]
    pub body: Vec<u8>,
}

#[binread]
#[derive(Debug)]
#[br(assert((0xD0..=0xD7).contains(&_raw_n)))]
#[br(assert(_marker == 0xFF))]
pub struct RSTSegment {
    #[br(temp)]
    pub _marker: u8,
    #[br(temp)]
    _raw_n: u8,
    #[br(calc = _raw_n & 0b0000_1111)]
    pub n: u8,
}

#[binread]
#[derive(Debug)]
#[br(magic = 0xFFDAu16)]
pub struct ScanHeader {
    pub len: u16,
    pub ns: u8,
    #[br(count = ns)]
    pub component_parameters: Vec<ScanComponentParameterSet>,
    pub ss: u8,
    pub se: u8,
    _raw_ah_al: u8,
    #[br(calc = (_raw_ah_al &0b1111_0000) >> 4)]
    pub ah: u8, //u4
    #[br(calc = _raw_ah_al &0b0000_1111)]
    pub al: u8, //u4
}

#[binread]
#[derive(Debug)]
pub struct ScanComponentParameterSet {
    pub cs: u8,
    #[br(temp)]
    _raw_td_ta: u8,
    #[br(calc = (_raw_td_ta &0b1111_0000) >> 4)]
    pub td: u8, //u4
    #[br(calc = _raw_td_ta &0b0000_1111)]
    pub ta: u8, //u4
}
