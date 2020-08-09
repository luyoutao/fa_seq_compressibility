use chrono::Local;
use env_logger::{self, Builder};
use flate2::{read, write, Compression};
use getopts::Options;
use log::{info, LevelFilter};
use std::env;
use std::fs::File;
use std::io::{self, BufRead, BufReader, Read, Write};
use std::path::Path;
use std::process::exit;

static VERSION: &str = "0.1.1";

struct Params {
    infile: String,
    outfile: Option<String>,
    seqlen: u32,
    infh: BufReader<File>,
    outfh: Box<dyn Write>,
    gzout: bool,
}

fn init_logger() {
    Builder::new()
        .format(|buf, record| {
            writeln!(
                buf,
                "[{} {}] {}",
                Local::now().format("%Y-%m-%d %H:%M:%S%.3f %z"),
                record.level(),
                record.args()
            )
        })
        .filter(None, LevelFilter::Info)
        .init();
}

fn usage(arg0: &str, opts: Options) {
    let s = format!(
        "\
Summary:
    Computes GZIP compressibility for genomic regions in every given interval (--seqlen)

Usage:
    {} --inFile hg38.fa --outfile output.bed [--seqlen 50] [--version|-v] [--help|-h]

Output:
    The output has 6 columns:
        1) chromosome name;
        2) start coordinate;
        3) end coordinate;
        4) sequence;
        5) GZIP compressibility ();
        6) genome strand;",
        arg0
    );
    eprintln!("{}", opts.usage(&s));
}

fn proc_args(args: &Vec<String>, mut opts: Options) -> Params {
    opts.optopt("i", "inFile", "", "input file in FASTA format");
    opts.optopt("o", "outFile", "", "output file; if omitted, write to STDOUT; otherwise, if ending with '.gz', will be GZ compressed");
    opts.optopt("l", "seqlen", "", "length (bp) of the intervals");
    opts.optflag("h", "help", "print usage");
    opts.optflag("v", "version", "print version");
    let matches = opts.parse(&args[1..]).unwrap();
    if matches.opt_present("h") {
        usage(&args[0], opts);
        exit(0);
    }
    if matches.opt_present("v") {
        println!("{} v{}", &args[0], &VERSION);
        exit(0);
    }
    let infile = match matches.opt_str("inFile") {
        Some(f) => match Path::new(&f).exists() {
            true => match &*(f
                .split('.')
                .last()
                .expect("Faied to find the file extension!")
                .to_lowercase())
            {
                "fa" => f,
                _ => panic!("{} does not seem to be a FASTA file!", f),
            },
            false => panic!("{} does not exist!", f),
        },
        None => panic!("--inFile is empty!"),
    };
    let infh = File::open(infile.as_str())
        .expect(&format!("Failed to open {} for read!", infile.as_str()));
    let infh = BufReader::new(infh);

    let outfile = matches.opt_str("outFile");

    let outfh = match &outfile {
        Some(s) => Box::new(File::create(s).expect(&format!("Cannot open {} for write!", s)))
            as Box<dyn Write>,
        None => Box::new(io::stdout()) as Box<dyn Write>,
    };

    let gzout = match &outfile {
        Some(s) => match &*(s
            .split('.')
            .last()
            .expect("Unknown file extension!")
            .to_lowercase())
        {
            "gz" => true,
            _ => false,
        },
        None => false,
    };

    let seqlen = match matches.opt_str("seqlen") {
        Some(s) => match s.parse::<u32>() {
            Ok(l) => l,
            Err(e) => panic!("--seqlen cannot be parsed: {}!", e),
        },
        _ => panic!("--seqlen incorrect!"),
    };

    let params = Params {
        infile: infile,
        outfile: outfile,
        seqlen: seqlen,
        infh: infh,
        outfh: outfh,
        gzout: gzout,
    };
    return params;
}

fn main() {
    init_logger();
    let args: Vec<String> = env::args().collect();
    let params = proc_args(&args, Options::new());
    proc_args(&args, Options::new());
    info!(
        "{{ infile = {}, outfile = {}, seqlen = {}, gzout = {}, VERSION = {} }}",
        &params.infile,
        match &params.outfile {
            Some(s) => s,
            None => "",
        },
        &params.seqlen,
        match &params.gzout {
            true => "true",
            false => "false",
        },
        VERSION
    );
    let infh = params.infh;
    let mut outfh = params.outfh;
    let seqlen = params.seqlen;
    let gzout = params.gzout;

    let mut chr: String = String::new();
    let mut seq = String::new();
    let mut i: u32 = 0;
    let mut _buf = vec![0u8; seqlen as usize];
    info!("Start processing FASTA...");
    for line in infh.lines() {
        let l = line.expect("Cannot read the current line!");
        match &l[..1] {
            ">" => {
                if !seq.is_empty() {
                    let s = &seq[..];
                    let mut gz = read::GzEncoder::new(s.as_bytes(), Compression::best());
                    let gzlen = match gz.read(&mut _buf) {
                        Ok(x) => x,
                        Err(e) => panic!(format!("Compression failed! {}", e)),
                    };
                    let r = seq.as_bytes().len() as f32 / (gzlen - 10) as f32;
                    if !gzout {
                        outfh
                            .write(
                                format!(
                                    "{}\t{}\t{}\t{}\t{}\t{}\n",
                                    chr,
                                    i * seqlen,
                                    (i + 1) * seqlen,
                                    s,
                                    r,
                                    "+"
                                )
                                .as_bytes(),
                            )
                            .unwrap();
                    } else {
                        let mut e = write::GzEncoder::new(Vec::new(), Compression::default());
                        e.write_all(
                            format!(
                                "{}\t{}\t{}\t{}\t{}\t{}\n",
                                chr,
                                i * seqlen,
                                (i + 1) * seqlen,
                                s,
                                r,
                                "+"
                            )
                            .as_bytes(),
                        )
                        .unwrap();
                        outfh.write(&e.finish().unwrap()).unwrap();
                    }
                    seq.clear();
                }
                chr = l[1..].to_string();
                i = 0;
            }
            _ => {
                seq.push_str(&l.to_uppercase());
                if seq.len() >= seqlen as usize {
                    let s = &seq[..seqlen as usize];
                    let mut gz = read::GzEncoder::new(s.as_bytes(), Compression::best());
                    let gzlen = match gz.read(&mut _buf) {
                        Ok(x) => x,
                        Err(e) => panic!(format!("Compression failed! {}", e)),
                    };
                    // after compression GACTTGCAGTGGGGGGA becomes
                    //          [1F,8B,08,00,00,00,00,00,02,FF,73,77,74,0E,09,71,77,76,74,0F,71,07,03,47,00]
                    //           -----------header------------                                              -----footer(CRC32)-----
                    // so we need to subtract 10-byte header
                    let r = seq.as_bytes().len() as f32 / (gzlen - 10) as f32;
                    if !gzout {
                        outfh
                            .write(
                                format!(
                                    "{}\t{}\t{}\t{}\t{}\t{}\n",
                                    chr,
                                    i * seqlen,
                                    (i + 1) * seqlen,
                                    s,
                                    r,
                                    "+"
                                )
                                .as_bytes(),
                            )
                            .unwrap(); // &_buf[0..gzlen] \n{:02X?}
                    } else {
                        let mut e = write::GzEncoder::new(Vec::new(), Compression::default());
                        e.write_all(
                            format!(
                                "{}\t{}\t{}\t{}\t{}\t{}\n",
                                chr,
                                i * seqlen,
                                (i + 1) * seqlen,
                                s,
                                r,
                                "+"
                            )
                            .as_bytes(),
                        )
                        .unwrap();
                        outfh.write(&e.finish().unwrap()).unwrap();
                    }
                    seq = seq.chars().skip(seqlen as usize).collect();
                    i += 1;
                }
            }
        };
    }
    info!("All done!");
}
