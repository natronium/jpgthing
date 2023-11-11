use std::{
    f64::consts::{FRAC_1_SQRT_2, PI},
    iter::once,
};

use bitvec::prelude::*;

use crate::{
    jpeg::{
        Frame, FrameComponentParameterSet, MiscSegment, QuantizationEntries,
        ScanComponentParameterSet,
    },
    jpeg_huff::{read_ac, read_dc_diff, BitReader, HuffmanInfo},
};

struct HuffmanTables {
    ac: [Option<HuffmanInfo>; 4],
    dc: [Option<HuffmanInfo>; 4],
}

pub struct Component {
    pub width: u16,
    pub height: u16,
    pub sampling_w: u8,
    pub sampling_h: u8,
    pub underlying_width: u16,
    pub underlying_height: u16,
    pub id: u8,         //TODO: ComponentID
    pub body: Vec<u16>, //TODO: ImageSample?
}

impl Component {
    //u8::try_from(component.id).unwrap(), unsigned_samples, frame.header.x, frame.header.y, frame_params.h, frame_params.v
    fn new(
        id: u8,
        data_units: Vec<UnsignedImageSamples>,
        component_params: &Vec<CombinedComponentParameters>,
        frame_width: u16,
        frame_height: u16,
        sampling_h: u8,
        sampling_w: u8,
    ) -> Self {
        //TODO? should i push the calculation out to the caller?
        let max_sampling_w = component_params
            .iter()
            .map(|params| params.horizontal_sampling)
            .max()
            .expect("No component params exist, so can't get max sampling factor");
        let max_sampling_h = component_params
            .iter()
            .map(|params| params.vertical_sampling)
            .max()
            .expect("No component params exist, so can't get max sampling factor");
        let width = (f64::from(frame_width) * (f64::from(sampling_w) / f64::from(max_sampling_w)))
            .ceil() as u16;
        let height = (f64::from(frame_height) * (f64::from(sampling_h) / f64::from(max_sampling_h)))
            .ceil() as u16;

        let mut thingy = ComponentSampleGridThingy::new(width, height, sampling_w, sampling_h);
        for (index, data_unit) in data_units.into_iter().enumerate() {
            thingy.insert_data_unit_at_idx(data_unit, index)
        }

        Component {
            width,
            height,
            sampling_w,
            sampling_h,
            underlying_width: thingy.underlying_width,
            underlying_height: thingy.underlying_height,
            id,
            body: thingy.get_data(),
        }
    }
}

#[allow(unused)]
struct ComponentSampleGridThingy {
    data: Vec<u16>,
    view_width: u16,
    view_height: u16,
    underlying_width: u16,
    underlying_height: u16,
    sampling_w: u8,
    sampling_h: u8,
}
impl ComponentSampleGridThingy {
    fn new(view_width: u16, view_height: u16, sampling_w: u8, sampling_h: u8) -> Self {
        let du_width = 8; // samples/du
        //sampling_w: du/chunk
        let chunk_width = du_width * sampling_w as usize; // samples/du * du/chunk = samples/chunk
        //view_width: samples
        let chunks_per_row = (f64::from(view_width) / f64::from(chunk_width as u32)).ceil() as usize; // samples / samples/chunk = chunks
        let samples_per_row = chunks_per_row * chunk_width; //chunks * samples/chunk = samples

        let du_height = 8;
        let chunk_height = du_height * sampling_h as usize;
        let chunks_per_col = (f64::from(view_height) / f64::from(chunk_height as u32)).ceil() as usize;
        let samples_per_col = chunks_per_col * chunk_height;

        ComponentSampleGridThingy {
            data: vec![0; samples_per_row * samples_per_col],
            view_width,
            view_height,
            underlying_width: samples_per_row
                .try_into()
                .expect("samples_per_row bigger than u16 MAX"),
            underlying_height: samples_per_col
                .try_into()
                .expect("samples_per_row bigger than u16 MAX"),
            sampling_w,
            sampling_h,
        }
    }

    fn idx_for_xy(&self, x: u16, y: u16) -> usize {
        self.underlying_width as usize * y as usize + x as usize
    }

    fn insert_sample_at(&mut self, sample: u16, x: u16, y: u16) {
        let pos = self.idx_for_xy(x, y);
        self.data[pos as usize] = sample;
    }

    fn insert_data_unit_at_idx(&mut self, data_unit: UnsignedImageSamples, index: usize) {
        let dataunits_per_chunk = (self.sampling_w * self.sampling_h) as usize;
        let dataunit_idx_in_chunk = index % dataunits_per_chunk as usize;
        let dataunit_x_in_chunk = dataunit_idx_in_chunk % self.sampling_w as usize;
        let dataunit_y_in_chunk =
            (dataunit_idx_in_chunk - dataunit_x_in_chunk) / self.sampling_w as usize;

        let chunk_idx =
            (f64::from(index as u32) / f64::from(dataunits_per_chunk as u32)).floor() as usize;
        let chunk_width = 8 * self.sampling_w as usize;
        let chunks_per_row =
            (f64::from(self.view_width) / f64::from(chunk_width as u32)).ceil() as usize;
        let chunk_height = 8 * self.sampling_h as usize;

        let chunk_x = chunk_idx % chunks_per_row;
        let chunk_y = (chunk_idx - chunk_x) / chunks_per_row;
        // chunk_idx = chunks_per_row * chunk_y + chunk_x

        let du_x = chunk_x * chunk_width + dataunit_x_in_chunk * 8;
        let du_y = chunk_y * chunk_height + dataunit_y_in_chunk * 8;

        self.insert_data_unit_at_xy(
            data_unit,
            du_x.try_into().unwrap(),
            du_y.try_into().unwrap(),
        );
    }

    fn insert_data_unit_at_xy(
        &mut self,
        UnsignedImageSamples(samples): UnsignedImageSamples,
        x: u16,
        y: u16,
    ) {
        for (du_y, row) in samples.into_iter().enumerate() {
            for (du_x, sample) in row.into_iter().enumerate() {
                self.insert_sample_at(sample, x + du_x as u16, y + du_y as u16);
            }
        }
    }

    fn get_data(self) -> Vec<u16> {
        self.data
    }
}

struct CombinedComponentParameters {
    id: u8,
    //frame
    horizontal_sampling: u8,   //u4
    vertical_sampling: u8,     //u4
    quantization_table_id: u8, //this one's actually a u8 for real (even though spec says only four quantization tables)
    //scan
    dc_huffman_table_id: u8, //u4
    ac_huffman_table_id: u8, //u4
}

impl CombinedComponentParameters {
    fn new(
        frame_components: &FrameComponentParameterSet,
        scan_components: &ScanComponentParameterSet,
    ) -> Self {
        assert!(
            frame_components.c == scan_components.cs,
            "frame and scan component parameters must share id to be combined"
        );
        CombinedComponentParameters {
            id: frame_components.c,
            horizontal_sampling: frame_components.h,
            vertical_sampling: frame_components.v,
            quantization_table_id: frame_components.tq,
            dc_huffman_table_id: scan_components.td,
            ac_huffman_table_id: scan_components.ta,
        }
    }
}

const TABLE_CLASS_DC: u8 = 0;
const TABLE_CLASS_AC: u8 = 1;

pub fn read_frame(frame: &Frame) -> Result<Vec<Component>, &'static str> {
    if frame.remaining_scans.len() != 0 {
        return Err("frame has more than one scan. don't currently know how to handle this");
    }

    let scan = &frame.scan_1;

    let mut misc_segments = frame.misc_tables.clone();
    misc_segments.append(&mut scan.misc_tables.clone());

    let huffman_tables = {
        let mut huffman_tables: HuffmanTables = HuffmanTables {
            ac: Default::default(),
            dc: Default::default(),
        };
        for segment in &misc_segments {
            if let MiscSegment::DefineHuffmanTable { len: _, tables } = segment {
                for table in tables {
                    if table.tc == TABLE_CLASS_DC {
                        huffman_tables.dc[table.th as usize] = Some(HuffmanInfo::new(&table));
                    } else if table.tc == TABLE_CLASS_AC {
                        huffman_tables.ac[table.th as usize] = Some(HuffmanInfo::new(&table));
                    } else {
                        return Err("Huffman table class must be 0 (DC) or 1 (AC)");
                    }
                }
            }
        }
        huffman_tables
    };

    let quantization_tables = {
        let mut quantization_tables: [Option<QuantizationEntries>; 4] = Default::default();
        for segment in &misc_segments {
            if let MiscSegment::DefineQuantizationTable { len: _, tables } = segment {
                for table in tables {
                    quantization_tables[table.tq as usize] = Some(table.qs.clone());
                }
            }
        }
        quantization_tables
    };

    /*
    B.2.2: The value of Nf shall be equal to the number of sets of frame component specification parameters (C_i , H_i , V_i , and Tq_i ) present in the frame header
    B.2.3: Each Cs_j shall match one of the C_i values specified in the frame header, and the ordering in the scan header shall follow the ordering in the frame header
    together means: we don't have any extraneous component params *and* they correspond between frame and scan, so we can zip them (TODO: is this true in practice? in the wild? of other decoders?)
     */
    let component_params: Vec<CombinedComponentParameters> = frame
        .header
        .component_parameters
        .iter()
        .zip(scan.header.component_parameters.iter())
        .map(|(frame_params, scan_params)| {
            CombinedComponentParameters::new(frame_params, scan_params)
        })
        .collect();

    let data_unit_counts: Vec<_> = component_params
        .iter()
        .map(|params| {
            (
                params.id,
                params.horizontal_sampling * params.vertical_sampling,
            )
        })
        .collect();

    let ecs_segments = scan
        .ecs_segments
        .iter()
        .map(|(ecs, _)| ecs)
        .chain(once(&scan.ecs_last));

    let componentwise_dataunits = ecs_segments
        .map(|ecs| {
            let mut raw_component_dataunits: Vec<Vec<RawDataUnit>> =
                (0..frame.header.nf.into()).map(|_| Vec::new()).collect(); //TODO: (rust learning) why can't i just use vec! here? allocation?? (it wants to Clone? but it's fine when i Vec::new myself??)

            let mut pred = vec![0; frame.header.nf as usize];
            let mut ecs_reader = BitReader::new(BitVec::<u8, Msb0>::from_slice(&ecs.body));

            let mut mcu_iter = data_unit_counts.iter().cycle();
            loop {
                let (component_id, du_count) =
                    mcu_iter.next().expect("iterator should be infinite???");

                let params: &CombinedComponentParameters =
                    &component_params[*component_id as usize - 1];
                let dc_table = huffman_tables.dc[params.dc_huffman_table_id as usize]
                    .as_ref()
                    .expect("specified dc huffman table for component does not exist");
                let ac_table = huffman_tables.ac[params.ac_huffman_table_id as usize]
                    .as_ref()
                    .expect("specified ac huffman table for component does not exist");

                for _ in 0..*du_count {
                    let raw_data_unit = read_raw_data_unit(
                        &mut ecs_reader,
                        dc_table,
                        ac_table,
                        pred[*component_id as usize - 1],
                    );
                    pred[*component_id as usize - 1] = raw_data_unit.0[0];
                    raw_component_dataunits[usize::from(*component_id) - 1].push(raw_data_unit);
                }

                if ecs_reader.is_empty() {
                    let quantization_entries = quantization_tables
                        [params.quantization_table_id as usize]
                        .as_ref()
                        .expect("specified quantization table for component does not exist");
                    let componentwise_image_samples: Vec<Vec<UnsignedImageSamples>> =
                        raw_component_dataunits
                            .into_iter()
                            .map(|component_dus| {
                                component_dus
                                    .into_iter()
                                    .map(|du| {
                                        let dequantized =
                                            dequantize_data_unit(du, &quantization_entries);
                                        let unzigzug = unzag_data_unit(dequantized);
                                        let signed_samples = undct_data_unit(unzigzug);
                                        unshift_data_unit(signed_samples, frame.header.p)
                                    })
                                    .collect::<Vec<UnsignedImageSamples>>()
                            })
                            .collect();
                    break componentwise_image_samples;
                }
            }
        })
        //TODO? /!\ assumes that component IDs *always* start at 0, and are contiguous, so correspond to our vector indices.
        //          i don't think this is guaranteed to be true by the spec though!
        .fold(
            (0..frame.header.nf.into())
                .map(|_| Vec::new())
                .collect::<Vec<Vec<UnsignedImageSamples>>>(),
            |mut acc, e| {
                for (id, mut component_dus) in e.into_iter().enumerate() {
                    acc[id].append(&mut component_dus)
                }
                acc
            },
        );

    let components: Vec<Component> = componentwise_dataunits
        .into_iter()
        .enumerate()
        .map(|(id, data_units)| {
            Component::new(
                u8::try_from(id).expect("id should be 4 bits. why does it not fit in a u8???"),
                data_units,
                &component_params,
                frame.header.x,
                frame.header.y,
                component_params[id].horizontal_sampling,
                component_params[id].horizontal_sampling,
            )
        })
        .collect();

    Ok(components)
}

///Quantized, ZigZagged, DCT coefficients, **with DC differences that need to be corrected**
struct RawDataUnit([i16; 64]);

/// read/decompress a data unit (8x8 style only for now) from an entropy-coded segment
/// returns a quantized, zig-zagged array of signed DCT coefficients
fn read_raw_data_unit(
    bitreader: &mut BitReader,
    dc_huffinfo: &HuffmanInfo,
    ac_huffinfo: &HuffmanInfo,
    dc_pred: i16,
) -> RawDataUnit {
    let diff = read_dc_diff(bitreader, dc_huffinfo);
    let dc = dc_pred + diff;
    let ac_coeffs = read_ac(bitreader, ac_huffinfo);

    //TODO: there must be a nicer way to do this? pass the array around to my reader functions?
    let mut coeffs = [0; 64];
    coeffs[0] = dc;
    coeffs[1..].clone_from_slice(&ac_coeffs);

    RawDataUnit(coeffs)
}

struct DequantizedCoeffs([i16; 64]);

/// dequantize a raw data unit
/// returns a **de**quantized, zig-zagged array of signed DCT coefficients
fn dequantize_data_unit(
    RawDataUnit(raw_data_unit): RawDataUnit,
    quantization_table: &QuantizationEntries,
) -> DequantizedCoeffs {
    let qt_vals = match quantization_table {
        QuantizationEntries::Low(vals) => vals.map(|n| n.into()),
        QuantizationEntries::Hi(vals) => vals.clone(),
    };

    //TODO: do we actually need to downconvert from 32 bits after doing our math?
    // worst cases are:
    //   i16::MIN * u16::MAX = -(2^15) * (2^16) = -2_147_483_648 = i32::MIN
    //   i16::MAX * u16::MAX = ((2^15)-1) * (2^16) = 2_147_418_112 < i32::MAX
    // so an i32 will hold the result of our calculation without issue.
    // HOWEVER: i'm not sure if i'm supposed to hang on to the precision here, or
    //  clamp back down into 8/12 bits since that's meant to be the sample precision
    let dequantized_vec: Vec<i16> = raw_data_unit
        .into_iter()
        .zip(qt_vals)
        .map(|(s, q)| -> i16 {
            let s: i32 = s.into();
            let q: i32 = q.into();
            let prod = s * q;
            prod.try_into()
                .unwrap_or(if s < 0 { i16::MIN } else { i16::MAX })
        })
        .collect();

    DequantizedCoeffs(std::array::from_fn(|n| dequantized_vec[n]))
}

//Figure A.6
const ZIGZAG: [[usize; 8]; 8] = [
    [0, 1, 5, 6, 14, 15, 27, 28],
    [2, 4, 7, 13, 16, 26, 29, 42],
    [3, 8, 12, 17, 25, 30, 41, 43],
    [9, 11, 18, 24, 31, 40, 44, 53],
    [10, 19, 23, 32, 39, 45, 52, 54],
    [20, 22, 33, 38, 46, 51, 55, 60],
    [21, 34, 37, 47, 50, 56, 59, 61],
    [35, 36, 48, 49, 57, 58, 62, 63],
];

struct Coeffs2D([[i16; 8]; 8]);

/// de-zig-zag a dequantized data unit into an actual 2D array
/// returns a dequantized **2D** aray of signed DCT coefficients
fn unzag_data_unit(DequantizedCoeffs(dequantized_data_unit): DequantizedCoeffs) -> Coeffs2D {
    let mut unzagged: [[i16; 8]; 8] = [[0; 8]; 8];
    for x in 0..8 {
        for y in 0..8 {
            unzagged[y][x] = dequantized_data_unit[ZIGZAG[y][x]];
        }
    }

    Coeffs2D(unzagged)
}

struct SignedImageSamples([[i16; 8]; 8]);

/// discrete cosine **un**transform our data unit
/// returns a 2D array of signed image samples
fn undct_data_unit(Coeffs2D(data_unit): Coeffs2D) -> SignedImageSamples {
    SignedImageSamples(idct(data_unit.map(|row| row.map(f64::from))))
}

pub struct UnsignedImageSamples(pub [[u16; 8]; 8]);

/// level shift the data unit from signed samples centered on 0 to unsigned samples which *start* at 0
/// returns a 2D array of unsigned image samples
fn unshift_data_unit(
    SignedImageSamples(data_unit): SignedImageSamples,
    sample_precision: u8,
) -> UnsignedImageSamples {
    UnsignedImageSamples(data_unit.map(|row| {
        row.map(|n| {
            let boop = n as i32 + 2i32.pow(sample_precision as u32 - 1);
            u16::try_from(boop).unwrap_or(u16::MAX)
        })
    }))
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

fn idct(samples: [[f64; 8]; 8]) -> [[i16; 8]; 8] {
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
            undct[y as usize][x as usize] = (0.25 * sum) as i16;
        }
    }

    undct
}
