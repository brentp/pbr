#[global_allocator]
static GLOBAL: mimalloc::MiMalloc = mimalloc::MiMalloc;

mod cached_faidx;
mod processor;
mod lua_filter;

use anyhow::Result;
use clap::Parser;
use processor::BasicProcessor;
use mlua::prelude::*;
use mlua::Function;
use perbase_lib::{
    par_granges,
    position::pileup_position::PileupPosition,
};
use std::path::PathBuf;

fn register_pile(lua: &Lua) -> mlua::Result<()> {
    lua.register_userdata_type::<PileupPosition>(|reg| {
        reg.add_field_method_get("depth", |_, this| Ok(this.depth));
        reg.add_field_method_get("a", |_, this| Ok(this.a));
        reg.add_field_method_get("c", |_, this| Ok(this.c));
        reg.add_field_method_get("g", |_, this| Ok(this.g));
        reg.add_field_method_get("t", |_, this| Ok(this.t));
        reg.add_field_method_get("n", |_, this| Ok(this.n));
        reg.add_field_method_get("fail", |_, this| Ok(this.fail));
        reg.add_field_method_get("ins", |_, this| Ok(this.ins));
        reg.add_field_method_get("del", |_, this| Ok(this.del));
        reg.add_field_method_get("ref_skip", |_, this| Ok(this.ref_skip));
        reg.add_field_method_get("pos", |_, this| Ok(this.pos));
    })
}

#[derive(Parser, Default, Debug)]
#[command(author, version, about, long_about = None)]
struct Args {
    #[arg(help = "Path to the bamfile")]
    bam_path: PathBuf,
    #[clap(short = 'e', long, default_value = "return true", help = "Lua expression to evaluate")]
    expression: String,
    #[clap(short, long, default_value = "2", help = "Number of threads to use")]
    threads: usize,
    #[clap(
        short,
        long,
        default_value_t = 100000,
        help = "maximum depth in the pileup"
    )]
    max_depth: u32,
    #[clap(short, long, help = "optional path to the BED of include regions")]
    bedfile: Option<PathBuf>,
    #[clap(short, long, help = "optional path to the reference fasta file")]
    fasta: Option<PathBuf>,
    #[clap(short = 'x', long, help = "optional path to BED of exclude regions")]
    exclude: Option<PathBuf>,

    #[clap(
        long,
        help = "adjust depth to not double count overlapping mates",
        long_help = "note that for now this is much slower than the default"
    )]
    mate_fix: bool,

    #[clap(short, long, help = "optional expression required for the pileup")]
    pile_expression: Option<String>,

    #[clap(
        short = 'k',
        long,
        default_value_t = 0,
        help = "number of flanking bases to fetch on each side of the reference base"
    )]
    flanking: usize,
}

fn main() -> Result<()> {
    let opts = Args::parse();

    let basic_processor = BasicProcessor {
        bamfile: PathBuf::from(&opts.bam_path),
        expression: String::from("") + opts.expression.as_str(),
        max_depth: opts.max_depth,
        exclude_regions: opts.exclude,
        mate_fix: opts.mate_fix,
        fasta_path: opts.fasta.clone(),
        flanking: opts.flanking,
    };

    let par_granges_runner = par_granges::ParGranges::new(
        opts.bam_path,      // pass in bam
        opts.fasta,         // optional ref fasta
        opts.bedfile,       // bedfile to narrow regions
        None,               // optional bcf/vcf file to specify positions of interest
        true,               // Merge any overlapping regions in the BED file
        Some(opts.threads), // optional allowed number of threads, defaults to max
        None,               // optional chunksize modification
        None,               // optional modifier on the size of the channel for sending Positions
        basic_processor,
    );

    let pile_lua = Lua::new();
    register_pile(&pile_lua)?;

    let pile_expression: Option<Function> = if let Some(expression) = opts.pile_expression {
        Some(pile_lua.load(expression.as_str()).into_function()?)
    } else {
        None
    };

    // Run the processor
    let receiver = par_granges_runner.process()?;
    println!("#chrom\tpos0\tref_base\tdepth\ta\tc\tg\tt\tn");
    // Pull the in-order results from the receiver channel
    receiver
        .into_iter()
        .filter(|(p, _)| p.depth > 0)
        // filter on the pile expression
        .filter(|(p, _)| {
            if let Some(pile_expression) = &pile_expression {
                let r = pile_lua.scope(|scope| {
                    let globals = pile_lua.globals();
                    let ud = scope.create_any_userdata_ref(p)?;
                    globals.set("pile", ud).expect("error setting pile");

                    pile_expression.call::<bool>(())
                });
                match r {
                    Ok(r) => r,
                    Err(e) => {
                        eprintln!("Error evaluating expression: {}", e);
                        std::process::exit(1);
                    }
                }
            } else {
                true
            }
        })
        .for_each(|(p, ref_seq)| {
            println!(
                "{chrom}\t{pos}\t{ref_base}\t{depth}\t{a}\t{c}\t{g}\t{t}\t{n}",
                chrom = p.ref_seq,
                pos = p.pos,
                depth = p.depth,
                ref_base = ref_seq.map(|seq| {
                    seq.chars().nth(opts.flanking).unwrap_or('.').to_string()
                }).unwrap_or_else(|| p.ref_base.unwrap_or('.').to_string()),
                a = p.a,
                c = p.c,
                g = p.g,
                t = p.t,
                n = p.n
            );
        });

    Ok(())
}

#[cfg(test)]
mod tests {
    use super::*;
    use mlua::Lua;
    use rust_htslib::bam;
    use rust_htslib::bam::pileup::Pileup;
    use rust_htslib::bam::record::Record;
    use rust_htslib::bam::{header::HeaderRecord, Header, HeaderView, IndexedReader, Read};
    use tempfile::NamedTempFile;
    use perbase_lib::read_filter::ReadFilter;

    #[test]
    fn test_read_bq() -> Result<()> {
        // Create a header with chr1
        let mut header = Header::new();
        let mut sq = HeaderRecord::new(b"SQ");
        sq.push_tag(b"SN", "chr1");
        sq.push_tag(b"LN", &1000000u32);
        header.push_record(&sq);
        let header_view = HeaderView::from_header(&header);

        // Create a test BAM record using SAM format
        let record = Record::from_sam(
            &header_view,
            b"test_read\t0\tchr1\t100\t30\t4M\t*\t0\t0\tACGT\t&&&&\tRG:Z:test", // Use standard ASCII for qual
        )
        .expect("Failed to create record from SAM");

        // Write record to a temporary BAM file
        let tmp = NamedTempFile::new()?;
        let path = tmp.path();
        {
            let mut writer = bam::Writer::from_path(path, &header, bam::Format::Bam)?;
            writer.write(&record)?;
        }
        bam::index::build(path, None, bam::index::Type::Bai, 1)?;

        // Create pileup
        let mut reader = IndexedReader::from_path(path)?;
        reader.fetch(("chr1", 100, 101))?; // Fetch the position of the record
        let mut pileups = reader.pileup();
        let pileup: Pileup = pileups.next().unwrap().expect("Failed to get pileup");

        // Get the first alignment from the pileup
        let alignment = pileup
            .alignments()
            .next()
            .expect("No alignment found in pileup");

        let lua = Lua::new();
        let rf = lua_filter::LuaReadFilter::new(
            "return read.bq > 0 and read.distance_from_5prime == 0 and read.distance_from_3prime > 0",
            &lua,
        )?; // Example expression

        // Test the bq functionality using the alignment from the pileup
        let result = rf.filter_read(&alignment.record(), Some(&alignment));
        assert!(result);

        Ok(())
    }

    #[test]
    fn test_pileup_position() -> mlua::Result<()> {
        let pileup_position = PileupPosition {
            depth: 10,
            a: 1,
            c: 2,
            g: 3,
            t: 4,
            n: 5,
            fail: 6,
            ins: 7,
            del: 8,
            ref_skip: 9,
            pos: 10,
            ..Default::default()
        };

        let lua = Lua::new();
        register_pile(&lua)?;
        let globals = lua.globals();
        for (expected, expression) in [
            (true, "pile.g > 3"),
            (true, "pile.a > 0"),
            (false, "pile.a > 10"),
            (false, "pile.ref_skip == 100"),
            (true, "pile.ref_skip == 9"),
        ] {
            eprintln!("Testing expression: {}", expression);
            lua.scope(|scope| {
                let p = scope
                    .create_any_userdata_ref(&pileup_position)
                    .expect("error creating user data");
                globals.set("pile", p)?;
                let f = lua
                    .load(&(String::from("return ") + expression))
                    .into_function()?;
                let result: bool = f.call(())?;
                Ok(result == expected)
            })
            .expect("error evaluating expression");
        }
        Ok(())
    }
}
