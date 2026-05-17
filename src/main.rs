use clap::Parser;
use csv::WriterBuilder;
use flate2::read::GzDecoder;
use half::f16;
use indexmap::IndexMap;
use serde::{Deserialize, Serialize};
use std::{
    collections::{HashMap, HashSet},
    error::Error,
    fs::File,
    io::BufReader,
    path::{Path, PathBuf},
};

#[derive(Parser, Debug)]
#[command(author, version, about, long_about = None)]
struct Args {
    /// Path to the variants file
    #[arg(long)]
    variants: PathBuf,

    /// Path to the tissue samples map file
    #[arg(long)]
    tissue_samples_map: PathBuf,

    /// Path to the transcript TPMs file
    #[arg(long)]
    transcript_tpms: PathBuf,

    /// Path to the output file (will be created if it does not exist)
    #[arg(long)]
    output: PathBuf,
}

#[derive(Deserialize, Serialize)]
struct RawVariant {
    #[serde(rename = "CHROM")]
    chr: String,

    #[serde(rename = "POS")]
    pos: usize,

    #[serde(rename = "REF")]
    reference: String,

    #[serde(rename = "ALT")]
    alternative: String,

    #[serde(rename = "vepFeature")]
    transcript: String,

    #[serde(rename = "vepSYMBOL")]
    gene: String,

    #[serde(rename = "vepConsequence")]
    consequence: String,

    #[serde(rename = "vepLoF")]
    loftee: String,
}

struct Transcript {
    id: String,
    gene: String,
    consequence: String,
    loftee: String,
}

struct Variant {
    chr: String,
    pos: usize,
    reference: String,
    alternative: String,
    transcripts: Vec<Transcript>,
}

struct AnnotatedVariant {
    chr: String,
    pos: usize,
    reference: String,
    alternative: String,
    transcripts: Vec<Transcript>,
    pext: Vec<f16>,
}

#[derive(Deserialize, Serialize)]
struct SampleTissue {
    #[serde(rename = "SAMPLE")]
    sample: String,

    #[serde(rename = "TISSUE")]
    tissue: String,
}

pub struct Tissue {
    pub tissue: String,
    pub samples: Vec<f16>,
}

pub type TsvMap = HashMap<String, Vec<Tissue>>;

type VariantKey = (String, usize, String, String);

struct SampleTPM {
    gene: String,
    consequence: String,
    tissue: String,
    loftee: String,
    tpm: Vec<f32>,
}

struct PextVariance {
    gene: String,
    consequence: String,
    loftee: String,
    variance: f32,
}

fn variant_key(row: &RawVariant) -> VariantKey {
    (
        row.chr.clone(),
        row.pos,
        row.reference.clone(),
        row.alternative.clone(),
    )
}

fn read_variants(path: PathBuf) -> Result<Vec<Variant>, Box<dyn Error>> {
    let mut reader = csv::ReaderBuilder::new()
        .delimiter(b'\t')
        .has_headers(true)
        .from_reader(BufReader::new(File::open(path)?));

    // IndexMap preserves the first-seen order of each variant.
    let mut map: IndexMap<VariantKey, Variant> = IndexMap::new();

    for result in reader.deserialize::<RawVariant>() {
        let row = result?;
        let key = variant_key(&row);

        let transcript = Transcript {
            id: row.transcript.clone(),
            gene: row.gene.clone(),
            consequence: row.consequence.clone(),
            loftee: row.loftee.clone(),
        };

        map.entry(key)
            .or_insert_with(|| Variant {
                chr: row.chr,
                pos: row.pos,
                reference: row.reference,
                alternative: row.alternative,
                transcripts: Vec::new(),
            })
            .transcripts
            .push(transcript);
    }

    Ok(map.into_values().collect())
}

fn read_sample_tissue_mapping(path: PathBuf) -> Result<Vec<SampleTissue>, Box<dyn Error>> {
    let mut reader = csv::ReaderBuilder::new()
        .delimiter(b'\t')
        .has_headers(true)
        .from_reader(BufReader::new(File::open(path)?));

    let rows: Vec<SampleTissue> = reader.deserialize().collect::<Result<Vec<_>, _>>()?;

    Ok(rows)
}

fn read_tsv(
    path: PathBuf,
    sample_tissues: &[SampleTissue],
) -> Result<TsvMap, Box<dyn std::error::Error>> {
    let file = File::open(path)?;
    let gz = GzDecoder::new(BufReader::new(file));

    let mut rdr = csv::ReaderBuilder::new()
        .delimiter(b'\t')
        .has_headers(true)
        .from_reader(gz);

    // Pre-compute unique tissue names in column order so we can build
    // the Vec<Tissue> skeleton once per row rather than re-hashing every cell
    let tissue_names: Vec<&str> = sample_tissues.iter().map(|st| st.tissue.as_str()).collect();

    // Stable, deduplicated list of tissue names preserving first-seen order
    let mut seen = std::collections::HashSet::new();
    let unique_tissues: Vec<&str> = tissue_names
        .iter()
        .copied()
        .filter(|t| seen.insert(*t))
        .collect();

    // Map tissue name → its index in the per-row Vec<Tissue>
    let tissue_idx: HashMap<&str, usize> = unique_tissues
        .iter()
        .enumerate()
        .map(|(i, &t)| (t, i))
        .collect();

    let mut map: TsvMap = HashMap::new();
    let mut record = csv::ByteRecord::new();

    while rdr.read_byte_record(&mut record)? {
        let key = std::str::from_utf8(&record[0])?.to_owned();

        // Build a fresh Vec<Tissue> skeleton for this row
        let mut tissues: Vec<Tissue> = unique_tissues
            .iter()
            .map(|&t| Tissue {
                tissue: t.to_owned(),
                samples: Vec::new(),
            })
            .collect();

        // Distribute each float into the correct Tissue bucket
        for (col_idx, st) in sample_tissues.iter().enumerate() {
            let raw = &record[col_idx + 1];
            let val = std::str::from_utf8(raw)?
                .parse::<f32>()
                .map(f16::from_f32)
                .unwrap_or(f16::NAN);

            let idx = tissue_idx[st.tissue.as_str()];
            tissues[idx].samples.push(val);
        }

        map.insert(key, tissues);
    }

    Ok(map)
}

fn calculate_variance(v: Vec<f32>) -> f32 {
    let len: f32 = v.len() as f32;
    let mean = v.iter().sum::<f32>() / len;
    return v.iter().map(|x| (x - mean).powi(2)).sum::<f32>() / len;
}

fn calculate_pext(variant: &Variant, matrix: &TsvMap) -> Result<Vec<PextVariance>, Box<dyn Error>> {
    let unique_genes: Vec<String> = variant
        .transcripts
        .iter()
        .map(|t| t.gene.to_string())
        .collect::<HashSet<_>>()
        .into_iter()
        .collect();

    let mut tpms: Vec<SampleTPM> = Vec::new();
    for gene in &unique_genes {
        for transcript in &variant.transcripts {
            if &transcript.gene != gene {
                continue;
            }

            if let Some(tissues) = matrix.get(&transcript.id) {
                for tissue in tissues {
                    tpms.push(SampleTPM {
                        gene: gene.clone(),
                        consequence: transcript.consequence.clone(),
                        tissue: tissue.tissue.clone(),
                        loftee: transcript.loftee.clone(),
                        tpm: tissue.samples.iter().map(|x| x.to_f32()).collect(),
                    });
                }
            }
        }
    }

    let unique_consequences: Vec<String> = tpms
        .iter()
        .map(|tpm| tpm.consequence.to_string())
        .collect::<HashSet<_>>()
        .into_iter()
        .collect();

    let unique_loftee: Vec<String> = tpms
        .iter()
        .map(|tpm| tpm.loftee.to_string())
        .collect::<HashSet<_>>()
        .into_iter()
        .collect();

    let unique_tissues: Vec<String> = tpms
        .iter()
        .map(|tpm| tpm.tissue.to_string())
        .collect::<HashSet<_>>()
        .into_iter()
        .collect();

    let mut variances: Vec<PextVariance> = Vec::new();
    for gene in &unique_genes {
        for tissue in &unique_tissues {
            let tissue_tpms: Vec<&Vec<f32>> = tpms
                .iter()
                .filter(|tpm| &tpm.gene == gene && &tpm.tissue == tissue)
                .map(|tpm| &tpm.tpm)
                .collect();

            let n_samples = tissue_tpms.get(0).ok_or("Found no TPMs.").iter().len();

            let total_tissue_tpms: Vec<f32> =
                tissue_tpms
                    .iter()
                    .fold(vec![0.0; n_samples], |mut acc, row| {
                        for (a, &x) in acc.iter_mut().zip(row.iter()) {
                            *a += x;
                        }
                        acc
                    });

            for consequence in &unique_consequences {
                for loftee in &unique_loftee {
                    let annotation_tpms: Vec<&Vec<f32>> = tpms
                        .iter()
                        .filter(|tpm| {
                            &tpm.gene == gene
                                && &tpm.consequence == consequence
                                && &tpm.loftee == loftee
                                && &tpm.tissue == tissue
                        })
                        .map(|tpm| &tpm.tpm)
                        .collect();

                    let n_samples = annotation_tpms.get(0).ok_or("Found no TPMs.").iter().len();

                    let total_annotation_tpms: Vec<f32> =
                        annotation_tpms
                            .iter()
                            .fold(vec![0.0; n_samples], |mut acc, row| {
                                for (a, &x) in acc.iter_mut().zip(row.iter()) {
                                    *a += x;
                                }
                                acc
                            });

                    let pext_scores: Vec<f32> = total_annotation_tpms
                        .iter()
                        .zip(&total_tissue_tpms)
                        .map(|(x, d)| x / d)
                        .collect();

                    let variance = calculate_variance(pext_scores);
                    variances.push(PextVariance {
                        gene: gene.clone(),
                        consequence: consequence.clone(),
                        loftee: loftee.clone(),
                        variance: variance,
                    });
                }
            }
        }
    }

    // let prop_tissue_tpms: Vec<Vec<f32>> = tissue_tpms
    // .iter()
    // .map(|row| {
    //     row.iter()
    //         .zip(total_tissue_tpms.iter())
    //         .map(|(&x, &d)| x / d)
    //         .collect()
    // })
    // .collect();

    Ok(variances)
}

fn main() -> Result<(), Box<dyn std::error::Error>> {
    let args = Args::parse();

    let variants: Vec<Variant> = read_variants(args.variants)?;
    let samples: Vec<SampleTissue> = read_sample_tissue_mapping(args.tissue_samples_map)?;
    let matrix = read_tsv(args.transcript_tpms, &samples)?;

    let mut writer = WriterBuilder::new()
        .delimiter(b'\t')
        .from_path(args.output)?;

    for variant in variants {
        let variances = calculate_pext(&variant, &matrix)?;

        for transcript in &variant.transcripts {
            let tissue_variances: Vec<_> = variances
                .iter()
                .filter(|v| {
                    v.gene == transcript.gene
                        && v.consequence == transcript.consequence
                        && v.loftee == transcript.loftee
                })
                .map(|v| v.variance.to_string())
                .collect();

            let joined_tissue_variances: String = tissue_variances.join("\t");

            writer.write_record(&[
                &variant.chr,
                &variant.pos.to_string(),
                &variant.reference,
                &variant.alternative,
                &transcript.gene,
                &transcript.loftee,
                &transcript.consequence,
                &joined_tissue_variances,
            ])?;
        }
    }

    Ok(())
}
