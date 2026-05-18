use clap::{Parser, ValueEnum};
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
    path::PathBuf,
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

    #[arg(long)]
    calculation: CalculationType,

    /// Path to the output file (will be created if it does not exist)
    #[arg(long)]
    output: PathBuf,
}

#[derive(Debug, Clone, PartialEq, ValueEnum)]
enum CalculationType {
    Simple,
    Annotated,
    Variance,
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

    #[serde(rename = "vepGene")]
    gene: String,

    #[serde(rename = "vepConsequence")]
    consequence: String,

    #[serde(rename = "vepLoF")]
    loftee: String,
}

struct TranscriptAnnotation {
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
    transcripts: Vec<TranscriptAnnotation>,
}

#[derive(Deserialize, Serialize)]
struct SampleTissue {
    #[serde(rename = "SAMPLE")]
    sample: String,

    #[serde(rename = "TISSUE")]
    tissue: String,
}

struct Tissue {
    tissue: String,
    samples: Vec<f16>,
}

struct Transcript {
    transcript_id: String,
    tissues: Vec<Tissue>,
}

type TsvMap = HashMap<String, Vec<Transcript>>;

type VariantKey = (String, usize, String, String);

struct AnnotatedTranscript {
    gene: String,
    consequence: String,
    tissue: String,
    loftee: String,
    tpm: Vec<f32>,
}

struct SimpleTranscript {
    gene: String,
    tissue: String,
    present: bool,
    tpm: Vec<f32>,
}

struct AnnotatedPextScore {
    gene: String,
    consequence: String,
    loftee: String,
    tissue: String,
    score: f16,
}

struct SimplePextScore {
    gene: String,
    tissue: String,
    score: f16,
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

        let transcript = TranscriptAnnotation {
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
        let gid = std::str::from_utf8(&record[1])?
            .split(".")
            .next()
            .ok_or("Could not split gene id")?
            .to_owned();

        let tid = std::str::from_utf8(&record[0])?
            .split(".")
            .next()
            .ok_or("Could not split transcript id")?
            .to_owned();

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
            let raw = &record[col_idx + 2];

            let val = std::str::from_utf8(raw)?
                .parse::<f32>()
                .map(f16::from_f32)?;

            let idx = tissue_idx[st.tissue.as_str()];
            tissues[idx].samples.push(val);
        }

        let t = Transcript {
            transcript_id: tid,
            tissues: tissues,
        };

        if let Some(transcripts) = map.get_mut(&gid) {
            transcripts.push(t);
            continue;
        }

        let mut transcripts: Vec<Transcript> = Vec::new();
        transcripts.push(t);

        map.insert(gid, transcripts);
    }

    Ok(map)
}

fn calculate_variance(v: Vec<f32>) -> f32 {
    let len: f32 = v.len() as f32;
    let mean = v.iter().sum::<f32>() / len;
    return v.iter().map(|x| (x - mean).powi(2)).sum::<f32>() / len;
}

fn get_unique_genes(variant: &Variant) -> Vec<String> {
    return variant
        .transcripts
        .iter()
        .map(|t| t.gene.to_string())
        .collect::<HashSet<_>>()
        .into_iter()
        .collect();
}

fn get_unique_consequences(tpms: &Vec<AnnotatedTranscript>) -> Vec<String> {
    return tpms
        .iter()
        .map(|tpm| tpm.consequence.to_string())
        .filter(|csq| csq != "")
        .collect::<HashSet<_>>()
        .into_iter()
        .collect();
}

fn get_unique_loftee(tpms: &Vec<AnnotatedTranscript>) -> Vec<String> {
    return tpms
        .iter()
        .map(|tpm| tpm.loftee.to_string())
        .filter(|loftee| loftee != "")
        .collect::<HashSet<_>>()
        .into_iter()
        .collect();
}

fn get_unique_tissues(tpms: &Vec<AnnotatedTranscript>) -> Vec<String> {
    return tpms
        .iter()
        .map(|tpm| tpm.tissue.to_string())
        .collect::<HashSet<_>>()
        .into_iter()
        .collect();
}

fn get_annotated_transcripts_for_variant(
    variant: &Variant,
    matrix: &TsvMap,
) -> Vec<AnnotatedTranscript> {
    let unique_genes = get_unique_genes(variant);

    let mut tpms: Vec<AnnotatedTranscript> = Vec::new();
    for gene in &unique_genes {
        if let Some(transcripts) = matrix.get(gene) {
            for transcript in transcripts {
                if let Some(variant_transcript) = variant
                    .transcripts
                    .iter()
                    .find(|t| transcript.transcript_id == t.id)
                {
                    for tissue in &transcript.tissues {
                        tpms.push(AnnotatedTranscript {
                            gene: gene.clone(),
                            consequence: variant_transcript.consequence.clone(),
                            tissue: tissue.tissue.clone(),
                            loftee: variant_transcript.loftee.clone(),
                            tpm: tissue.samples.iter().map(|x| x.to_f32()).collect(),
                        });
                    }

                    continue;
                }

                for tissue in &transcript.tissues {
                    tpms.push(AnnotatedTranscript {
                        gene: gene.clone(),
                        consequence: "".to_string(),
                        tissue: tissue.tissue.clone(),
                        loftee: "".to_string(),
                        tpm: tissue.samples.iter().map(|x| x.to_f32()).collect(),
                    });
                }
            }
        }
    }

    return tpms;
}

fn get_simple_transcripts_for_variant(variant: &Variant, matrix: &TsvMap) -> Vec<SimpleTranscript> {
    let unique_genes = get_unique_genes(variant);

    let mut tpms: Vec<SimpleTranscript> = Vec::new();
    for gene in &unique_genes {
        if let Some(transcripts) = matrix.get(gene) {
            for transcript in transcripts {
                if let Some(_) = variant
                    .transcripts
                    .iter()
                    .find(|t| transcript.transcript_id == t.id)
                {
                    for tissue in &transcript.tissues {
                        tpms.push(SimpleTranscript {
                            gene: gene.clone(),
                            tissue: tissue.tissue.clone(),
                            present: true,
                            tpm: tissue.samples.iter().map(|x| x.to_f32()).collect(),
                        });
                    }

                    continue;
                }

                for tissue in &transcript.tissues {
                    tpms.push(SimpleTranscript {
                        gene: gene.clone(),
                        tissue: tissue.tissue.clone(),
                        present: false,
                        tpm: tissue.samples.iter().map(|x| x.to_f32()).collect(),
                    });
                }
            }
        }
    }

    return tpms;
}

// TODO: Use a math library for vector math instead of all of this
fn matrix_vertical_sum(matrix: &Vec<&Vec<f32>>, n_columns: usize) -> Vec<f32> {
    return matrix.iter().fold(vec![0.0; n_columns], |mut acc, row| {
        for (a, &x) in acc.iter_mut().zip(row.iter()) {
            *a += x;
        }
        acc
    });
}

fn calculate_pext_variance(
    variant: &Variant,
    matrix: &TsvMap,
) -> Result<Vec<AnnotatedPextScore>, Box<dyn Error>> {
    let tpms: Vec<AnnotatedTranscript> = get_annotated_transcripts_for_variant(variant, matrix);

    let unique_genes = get_unique_genes(variant);
    let unique_consequences = get_unique_consequences(&tpms);
    let unique_loftee = get_unique_loftee(&tpms);
    let unique_tissues = get_unique_tissues(&tpms);

    let mut variances: Vec<AnnotatedPextScore> = Vec::new();
    for gene in &unique_genes {
        for tissue in &unique_tissues {
            let tissue_tpms: Vec<&Vec<f32>> = tpms
                .iter()
                .filter(|tpm| &tpm.gene == gene && &tpm.tissue == tissue)
                .map(|tpm| &tpm.tpm)
                .collect();

            let n_samples = tissue_tpms.get(0).ok_or("Found no TPMs.")?.len();
            let total_tissue_tpms: Vec<f32> = matrix_vertical_sum(&tissue_tpms, n_samples);

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

                    // If this combination of annotations doesn't exist
                    if annotation_tpms.len() == 0 {
                        continue;
                    }

                    let total_annotation_tpms: Vec<f32> =
                        matrix_vertical_sum(&annotation_tpms, n_samples);

                    // TODO: Again, math library!
                    let pext_scores: Vec<f32> = total_annotation_tpms
                        .iter()
                        .zip(&total_tissue_tpms)
                        .map(|(x, d)| x / d)
                        .collect();

                    let variance = f16::from_f32(calculate_variance(pext_scores));

                    if variance.is_nan() {
                        variances.push(AnnotatedPextScore {
                            gene: gene.clone(),
                            consequence: consequence.clone(),
                            loftee: loftee.clone(),
                            tissue: tissue.clone(),
                            score: f16::from_f32(0.0),
                        });
                        continue;
                    }

                    variances.push(AnnotatedPextScore {
                        gene: gene.clone(),
                        consequence: consequence.clone(),
                        loftee: loftee.clone(),
                        tissue: tissue.clone(),
                        score: variance,
                    });
                }
            }
        }
    }

    Ok(variances)
}

fn calculate_pext(
    variant: &Variant,
    matrix: &TsvMap,
) -> Result<Vec<SimplePextScore>, Box<dyn Error>> {
    let tpms: Vec<SimpleTranscript> = get_simple_transcripts_for_variant(variant, matrix);

    let unique_genes = get_unique_genes(variant);
    let unique_tissues: Vec<String> = tpms
        .iter()
        .map(|tpm| tpm.tissue.to_string())
        .collect::<HashSet<_>>()
        .into_iter()
        .collect();

    let mut variances: Vec<SimplePextScore> = Vec::new();
    for gene in &unique_genes {
        for tissue in &unique_tissues {
            let total_transcripts: Vec<&Vec<f32>> = tpms
                .iter()
                .filter(|tpm| &tpm.gene == gene && &tpm.tissue == tissue)
                .map(|tpm| &tpm.tpm)
                .collect();

            let n_samples = total_transcripts.get(0).ok_or("Found no TPMs.")?.len();
            let total_tissue_tpms: Vec<f32> = matrix_vertical_sum(&total_transcripts, n_samples);

            let present_transcripts: Vec<&Vec<f32>> = tpms
                .iter()
                .filter(|tpm| &tpm.gene == gene && &tpm.tissue == tissue && tpm.present == true)
                .map(|tpm| &tpm.tpm)
                .collect();

            let present_tissue_tpms: Vec<f32> =
                matrix_vertical_sum(&present_transcripts, n_samples);

            // TODO: Again, math library!
            let pext_scores: Vec<f32> = present_tissue_tpms
                .iter()
                .zip(&total_tissue_tpms)
                .map(|(x, d)| x / d)
                .collect();

            let pext: f16 = f16::from_f32(pext_scores.iter().sum());

            if pext.is_nan() {
                variances.push(SimplePextScore {
                    gene: gene.clone(),
                    tissue: tissue.clone(),
                    score: f16::from_f32(0.0),
                });
                continue;
            }

            variances.push(SimplePextScore {
                gene: gene.clone(),
                tissue: tissue.clone(),
                score: pext,
            });
        }
    }

    Ok(variances)
}

fn calculate_annotated_pext(
    variant: &Variant,
    matrix: &TsvMap,
) -> Result<Vec<AnnotatedPextScore>, Box<dyn Error>> {
    let tpms: Vec<AnnotatedTranscript> = get_annotated_transcripts_for_variant(variant, matrix);

    let unique_genes = get_unique_genes(variant);
    let unique_consequences = get_unique_consequences(&tpms);
    let unique_loftee = get_unique_loftee(&tpms);
    let unique_tissues = get_unique_tissues(&tpms);

    let mut variances: Vec<AnnotatedPextScore> = Vec::new();
    for gene in &unique_genes {
        for tissue in &unique_tissues {
            let tissue_tpms: Vec<&Vec<f32>> = tpms
                .iter()
                .filter(|tpm| &tpm.gene == gene && &tpm.tissue == tissue)
                .map(|tpm| &tpm.tpm)
                .collect();

            let n_samples = tissue_tpms.get(0).ok_or("Found no TPMs.")?.len();
            let total_tissue_tpms: Vec<f32> = matrix_vertical_sum(&tissue_tpms, n_samples);

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

                    // If this combination of annotations doesn't exist
                    if annotation_tpms.len() == 0 {
                        continue;
                    }

                    let total_annotation_tpms: Vec<f32> =
                        matrix_vertical_sum(&annotation_tpms, n_samples);

                    // TODO: Again, math library!
                    let pext_scores: Vec<f32> = total_annotation_tpms
                        .iter()
                        .zip(&total_tissue_tpms)
                        .map(|(x, d)| x / d)
                        .collect();

                    let pext: f16 = f16::from_f32(pext_scores.iter().sum());

                    if pext.is_nan() {
                        variances.push(AnnotatedPextScore {
                            gene: gene.clone(),
                            consequence: consequence.clone(),
                            loftee: loftee.clone(),
                            tissue: tissue.clone(),
                            score: f16::from_f32(0.0),
                        });
                        continue;
                    }

                    variances.push(AnnotatedPextScore {
                        gene: gene.clone(),
                        consequence: consequence.clone(),
                        loftee: loftee.clone(),
                        tissue: tissue.clone(),
                        score: pext,
                    });
                }
            }
        }
    }

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
        if args.calculation == CalculationType::Variance {
            let score = calculate_pext_variance(&variant, &matrix)?;

            for variance in score {
                writer.write_record(&[
                    &variant.chr,
                    &variant.pos.to_string(),
                    &variant.reference,
                    &variant.alternative,
                    &variance.gene,
                    &variance.loftee,
                    &variance.consequence,
                    &variance.tissue,
                    &variance.score.to_string(),
                ])?;
            }
        } else if args.calculation == CalculationType::Simple {
            let score = calculate_pext(&variant, &matrix)?;

            for variance in score {
                writer.write_record(&[
                    &variant.chr,
                    &variant.pos.to_string(),
                    &variant.reference,
                    &variant.alternative,
                    &variance.gene,
                    &variance.tissue,
                    &variance.score.to_string(),
                ])?;
            }
        } else if args.calculation == CalculationType::Annotated {
            let score = calculate_annotated_pext(&variant, &matrix)?;

            for variance in score {
                writer.write_record(&[
                    &variant.chr,
                    &variant.pos.to_string(),
                    &variant.reference,
                    &variant.alternative,
                    &variance.gene,
                    &variance.loftee,
                    &variance.consequence,
                    &variance.tissue,
                    &variance.score.to_string(),
                ])?;
            }
        }
    }

    Ok(())
}
