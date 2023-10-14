use std::{
    io::{Read, Seek, SeekFrom},
    iter::from_fn,
};

use binrw::{parser, BinRead, BinResult, Endian};

//stolen and readapted from binrw::helpers::{until_eof, until_eof_with}
//TODO: understand how this works and why i couldn't write a non-clojurey generic version
pub(crate) fn until_invalid<Reader, T, Arg, Ret>(
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

pub(crate) fn until_invalid_with<Reader, T, Arg, ReadFn, Ret>(
    read: ReadFn,
) -> impl Fn(&mut Reader, Endian, Arg) -> BinResult<Ret>
where
    Reader: Read + Seek,
    Arg: Clone,
    ReadFn: Fn(&mut Reader, Endian, Arg) -> BinResult<T>,
    Ret: FromIterator<T>,
{
    move |reader, endian, args| {
        from_fn(|| {
            let initial_pos = reader
                .stream_position()
                .expect("reader should have a stream position");
            match read(reader, endian, args.clone()) {
                ok @ Ok(_) => Some(ok),
                Err(_) => {
                    reader
                        .seek(SeekFrom::Start(initial_pos))
                        .expect("reader should be able to seek back to old position");
                    None
                }
            }
        })
        .fuse()
        .collect()
    }
}

//TODO i'm pretty sure there's a more idiomatic way to do this
//TODO i'm pretty sure there's also a faster way to do this
#[parser(reader)]
pub(crate) fn parse_ecs() -> BinResult<Vec<u8>> {
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
