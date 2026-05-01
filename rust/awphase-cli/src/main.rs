use std::fs;
use std::path::{Path, PathBuf};
use std::time::Instant;

use anyhow::{Context, Result};
use awphase_core::bridge::block_summary::{
    boundary_scores_from_blocks, BlockSummary, BoundaryScoreInput,
};
use awphase_core::bridge::blockify::make_blocks_from_local_calls;
use awphase_core::bridge::boundary::score_boundary;
use awphase_core::bridge::debug::BoundaryScoreDebug;
use awphase_core::bridge::project::project_stitched_phase_to_calls;
use awphase_core::bridge::stitch::stitch_blocks;
use awphase_core::bridge::summarize::summarize_blocks_from_local_calls;
use awphase_core::candidate::pbwt_v1::PanelDonorScoringEngine;
use awphase_core::confidence::interface::{BoundaryDecision, Decision};
use awphase_core::config::BenchmarkConfig;
use awphase_core::decoder::local::ReadDonorLocalDecoder;
use awphase_core::decoder::two_pass::run_two_pass_local_decode;
use awphase_core::io::block_summary::load_block_summaries_json;
use awphase_core::io::panel::load_panel_haplotypes;
use awphase_core::io::reads::load_reads_json;
use awphase_core::io::site_priors::load_site_priors_json;
use awphase_core::io::vcf::load_variants_json;
use awphase_core::report::RunReport;
use awphase_core::summary::RunSummary;
use awphase_core::types::{AncestryPrior, SiteDonorPriors, WindowId};
use clap::{Parser, Subcommand};

#[derive(Parser)]
#[command(name = "awphase")]
#[command(version = "0.1.0")]
#[command(about = "awphase starter CLI")]
struct Cli {
    #[arg(long, default_value = "configs/benchmark_chr20.yaml")]
    config: String,

    #[command(subcommand)]
    command: Commands,
}

#[derive(Subcommand)]
enum Commands {
    Candidates {
        #[arg(long)]
        region: Option<String>,
        #[arg(long)]
        ancestry_off: bool,
    },
    Local {
        #[arg(long)]
        region: Option<String>,
        #[arg(long)]
        ancestry_off: bool,
    },
    Bridge {
        #[arg(long)]
        same_score: Option<f32>,
        #[arg(long)]
        flip_score: Option<f32>,
        #[arg(long)]
        informative_sites: Option<usize>,
        #[arg(long)]
        ancestry_off: bool,
    },
    Benchmark {
        #[arg(long)]
        region: Option<String>,
        #[arg(long)]
        ancestry_off: bool,
    },
}

fn parse_region(region: &str) -> Result<WindowId> {
    let mut parts = region.split(':');
    let chrom = parts
        .next()
        .context("region missing chromosome, expected format chr:start-end")?
        .to_string();
    let range = parts
        .next()
        .context("region missing coordinates, expected format chr:start-end")?;
    let mut bounds = range.split('-');
    let start = bounds
        .next()
        .context("region missing start coordinate")?
        .replace(',', "")
        .parse::<u64>()?;
    let end = bounds
        .next()
        .context("region missing end coordinate")?
        .replace(',', "")
        .parse::<u64>()?;

    Ok(WindowId { chrom, start, end })
}

fn load_config(path: &str) -> Result<BenchmarkConfig> {
    let text =
        fs::read_to_string(path).with_context(|| format!("failed to read config file: {path}"))?;
    let cfg: BenchmarkConfig = serde_yaml::from_str(&text)
        .with_context(|| format!("failed to parse YAML config: {path}"))?;
    Ok(cfg)
}

fn ancestry_from_config(cfg: &BenchmarkConfig) -> AncestryPrior {
    AncestryPrior {
        groups: cfg.ancestry.groups.clone(),
        weights: cfg.ancestry.weights.clone(),
    }
}

fn ancestry_disabled(cfg: &BenchmarkConfig) -> AncestryPrior {
    let n = cfg.ancestry.groups.len().max(1);
    let uniform = 1.0_f32 / (n as f32);

    AncestryPrior {
        groups: cfg.ancestry.groups.clone(),
        weights: vec![uniform; cfg.ancestry.groups.len()],
    }
}

fn site_priors_path(cfg: &BenchmarkConfig, ancestry_off: bool) -> String {
    if ancestry_off {
        if cfg.sample_id.contains("messy") {
            "data/scenarios/messy_middle_block/hg002_chr20_site_priors_ancestry_off.json"
                .to_string()
        } else if cfg.sample_id.contains("noisy") {
            "data/scenarios/noisy_but_recoverable/hg002_chr20_site_priors_ancestry_off.json"
                .to_string()
        } else if cfg.sample_id.contains("clean") {
            "data/scenarios/clean_blocks/hg002_chr20_site_priors_ancestry_off.json".to_string()
        } else {
            "data/derived/hg002_chr20_site_priors_ancestry_off.json".to_string()
        }
    } else {
        cfg.inputs.target_site_priors.clone()
    }
}

fn block_summary_path(cfg: &BenchmarkConfig, ancestry_off: bool) -> Option<String> {
    if ancestry_off {
        Some("data/derived/hg002_chr20_blocks_ancestry_off.json".to_string())
    } else {
        cfg.inputs.target_block_summaries.clone()
    }
}

fn panel_path_for_config(cfg: &BenchmarkConfig) -> String {
    let json_path = if !cfg.inputs.panel_path.is_empty() {
        cfg.inputs.panel_path.clone()
    } else if cfg.sample_id.contains("messy") {
        "data/scenarios/messy_middle_block/hg002_chr20_panel_candidates.json".to_string()
    } else if cfg.sample_id.contains("noisy") {
        "data/scenarios/noisy_but_recoverable/hg002_chr20_panel_candidates.json".to_string()
    } else if cfg.sample_id.contains("clean") {
        "data/scenarios/clean_blocks/hg002_chr20_panel_candidates.json".to_string()
    } else {
        "data/derived/panel_chr20_candidates.json".to_string()
    };

    let bin_path = if let Some(stripped) = json_path.strip_suffix(".json") {
        format!("{stripped}.bin")
    } else {
        format!("{json_path}.bin")
    };

    if Path::new(&bin_path).exists() {
        bin_path
    } else {
        json_path
    }
}

fn ensure_results_dir(sample_id: &str, ancestry_off: bool) -> Result<PathBuf> {
    let suffix = if ancestry_off {
        "ancestry_off"
    } else {
        "ancestry_on"
    };
    let dir = PathBuf::from(format!("results/runs/{sample_id}/{suffix}"));
    fs::create_dir_all(&dir)
        .with_context(|| format!("failed to create results directory: {}", dir.display()))?;
    Ok(dir)
}

fn ensure_bridge_dir(sample_id: &str) -> Result<PathBuf> {
    let dir = PathBuf::from(format!("results/runs/{sample_id}/bridge"));
    fs::create_dir_all(&dir)
        .with_context(|| format!("failed to create bridge directory: {}", dir.display()))?;
    Ok(dir)
}

fn derive_boundary_scores(blocks: &[BlockSummary]) -> Vec<BoundaryScoreInput> {
    blocks
        .windows(2)
        .map(|pair| boundary_scores_from_blocks(&pair[0], &pair[1]))
        .collect()
}

fn boundaries_from_scores(
    score_inputs: &[BoundaryScoreInput],
    cfg: &BenchmarkConfig,
) -> Vec<BoundaryDecision> {
    score_inputs
        .iter()
        .map(|score_in| {
            score_boundary(
                &score_in.left_block,
                &score_in.right_block,
                score_in.same_score,
                score_in.flip_score,
                cfg.bridge.min_confidence,
                cfg.bridge.min_informative_sites,
                score_in.informative_sites,
            )
        })
        .collect()
}

fn combine_boundary_debug(
    score_inputs: &[BoundaryScoreInput],
    decisions: &[BoundaryDecision],
) -> Vec<BoundaryScoreDebug> {
    score_inputs
        .iter()
        .zip(decisions.iter())
        .map(|(s, d)| BoundaryScoreDebug {
            left_block: s.left_block.clone(),
            right_block: s.right_block.clone(),
            same_score: s.same_score,
            flip_score: s.flip_score,
            informative_sites: s.informative_sites,
            decision: format!("{:?}", d.decision),
            confidence: d.confidence.value,
            reason: d.confidence.reason.clone(),
        })
        .collect()
}

fn maybe_debug_block_summaries(
    cfg: &BenchmarkConfig,
    ancestry_off: bool,
) -> Result<Option<Vec<BlockSummary>>> {
    match block_summary_path(cfg, ancestry_off) {
        Some(path) => Ok(Some(load_block_summaries_json(&path)?)),
        None => Ok(None),
    }
}

struct LoadedInputs {
    variants: Vec<awphase_core::types::VariantSite>,
    reads: Vec<awphase_core::types::ReadObservation>,
    site_priors: SiteDonorPriors,
    panel_haplotypes: Vec<awphase_core::candidate::interface::PanelHaplotype>,
}

fn load_inputs(
    cfg: &BenchmarkConfig,
    window: &WindowId,
    ancestry_off: bool,
) -> Result<LoadedInputs> {
    let variants = load_variants_json(&cfg.inputs.target_variants, window)?;
    let reads = load_reads_json(&cfg.inputs.target_reads, window)?;
    let site_priors = load_site_priors_json(&site_priors_path(cfg, ancestry_off), window)?;
    let panel_haplotypes = load_panel_haplotypes(&panel_path_for_config(cfg))?;

    Ok(LoadedInputs {
        variants,
        reads,
        site_priors,
        panel_haplotypes,
    })
}

fn main() -> Result<()> {
    let cli = Cli::parse();
    let cfg = load_config(&cli.config)?;

    match cli.command {
        Commands::Candidates {
            region,
            ancestry_off,
        } => {
            let ancestry = if ancestry_off {
                ancestry_disabled(&cfg)
            } else {
                ancestry_from_config(&cfg)
            };

            let region_str = region.unwrap_or_else(|| cfg.region.clone());
            let window = parse_region(&region_str)?;
            let inputs = load_inputs(&cfg, &window, ancestry_off)?;
            let engine = PanelDonorScoringEngine {
                panel_path: panel_path_for_config(&cfg),
                top_k: cfg.candidate_engine.top_k_donors,
                max_composite_states: cfg.candidate_engine.max_composite_states,
            };
            let decoder = ReadDonorLocalDecoder;

            let two_pass = run_two_pass_local_decode(
                &engine,
                &decoder,
                &window,
                &inputs.variants,
                &inputs.reads,
                &ancestry,
                &inputs.site_priors,
                &inputs.panel_haplotypes,
            )?;

            println!("{}", serde_json::to_string_pretty(&two_pass.candidates)?);
        }
        Commands::Local {
            region,
            ancestry_off,
        } => {
            let ancestry = if ancestry_off {
                ancestry_disabled(&cfg)
            } else {
                ancestry_from_config(&cfg)
            };

            let region_str = region.unwrap_or_else(|| cfg.region.clone());
            let window = parse_region(&region_str)?;
            let inputs = load_inputs(&cfg, &window, ancestry_off)?;
            let engine = PanelDonorScoringEngine {
                panel_path: panel_path_for_config(&cfg),
                top_k: cfg.candidate_engine.top_k_donors,
                max_composite_states: cfg.candidate_engine.max_composite_states,
            };
            let decoder = ReadDonorLocalDecoder;

            let two_pass = run_two_pass_local_decode(
                &engine,
                &decoder,
                &window,
                &inputs.variants,
                &inputs.reads,
                &ancestry,
                &inputs.site_priors,
                &inputs.panel_haplotypes,
            )?;

            println!("{}", serde_json::to_string_pretty(&two_pass.final_result)?);
        }
        Commands::Bridge {
            same_score,
            flip_score,
            informative_sites,
            ancestry_off,
        } => {
            let decision = if let (Some(same), Some(flip), Some(info)) =
                (same_score, flip_score, informative_sites)
            {
                score_boundary(
                    "block_1",
                    "block_2",
                    same,
                    flip,
                    cfg.bridge.min_confidence,
                    cfg.bridge.min_informative_sites,
                    info,
                )
            } else {
                let blocks = maybe_debug_block_summaries(&cfg, ancestry_off)?
                    .context("no debug block summaries available for bridge command")?;
                let score_inputs = derive_boundary_scores(&blocks);
                let derived = boundaries_from_scores(&score_inputs, &cfg);
                derived
                    .into_iter()
                    .next()
                    .context("no block boundaries available")?
            };

            let out_dir = ensure_bridge_dir(&cfg.sample_id)?;
            let out_json = out_dir.join("boundary_decision.json");
            fs::write(&out_json, serde_json::to_string_pretty(&decision)?)
                .with_context(|| format!("failed to write {}", out_json.display()))?;

            println!("{}", serde_json::to_string_pretty(&decision)?);
        }
        Commands::Benchmark {
            region,
            ancestry_off,
        } => {
            let ancestry = if ancestry_off {
                ancestry_disabled(&cfg)
            } else {
                ancestry_from_config(&cfg)
            };

            let timer = Instant::now();
            let region_str = region.unwrap_or_else(|| cfg.region.clone());
            let window = parse_region(&region_str)?;
            let inputs = load_inputs(&cfg, &window, ancestry_off)?;
            let engine = PanelDonorScoringEngine {
                panel_path: panel_path_for_config(&cfg),
                top_k: cfg.candidate_engine.top_k_donors,
                max_composite_states: cfg.candidate_engine.max_composite_states,
            };
            let decoder = ReadDonorLocalDecoder;

            let two_pass = run_two_pass_local_decode(
                &engine,
                &decoder,
                &window,
                &inputs.variants,
                &inputs.reads,
                &ancestry,
                &inputs.site_priors,
                &inputs.panel_haplotypes,
            )?;

            let result = two_pass.final_result;
            let candidates = two_pass.candidates;
            let local_debug = two_pass.final_debug;

            let block_ranges = make_blocks_from_local_calls(&result.calls);
            let derived_blocks = summarize_blocks_from_local_calls(&result.calls, &block_ranges);
            let boundary_scores = derive_boundary_scores(&derived_blocks);
            let boundaries = boundaries_from_scores(&boundary_scores, &cfg);
            let boundary_debug = combine_boundary_debug(&boundary_scores, &boundaries);
            let stitched = stitch_blocks(
                &derived_blocks,
                &boundaries
                    .iter()
                    .map(|b| b.decision.clone())
                    .collect::<Vec<_>>(),
            );
            let final_stitched_calls =
                project_stitched_phase_to_calls(&result.calls, &block_ranges, &stitched.blocks);

            let number_of_sites = final_stitched_calls.calls.len();
            let number_phased = final_stitched_calls
                .calls
                .iter()
                .filter(|c| c.stitched_phase_state != 0)
                .count();
            let number_abstained = number_of_sites.saturating_sub(number_phased);

            let informative_boundaries = boundaries
                .iter()
                .filter(|b| b.confidence.informative_sites >= cfg.bridge.min_informative_sites)
                .count();
            let joined_boundaries = boundaries
                .iter()
                .filter(|b| matches!(b.decision, Decision::Join))
                .count();
            let flipped_boundaries = boundaries
                .iter()
                .filter(|b| matches!(b.decision, Decision::Flip))
                .count();
            let abstained_boundaries = boundaries
                .iter()
                .filter(|b| matches!(b.decision, Decision::Abstain))
                .count();

            let summary = RunSummary {
                sample_id: cfg.sample_id.clone(),
                region: region_str.clone(),
                decoder_mode: if ancestry_off {
                    "read_plus_donor_ancestry_off".to_string()
                } else {
                    "read_plus_donor_ancestry_on".to_string()
                },
                number_of_sites,
                number_phased,
                number_abstained,
                informative_boundaries,
                joined_boundaries,
                abstained_boundaries,
                flipped_boundaries,
                runtime_seconds: timer.elapsed().as_secs_f32(),
            };

            let report = RunReport {
                summary: summary.clone(),
                derived_block_summaries: derived_blocks.clone(),
                derived_boundary_scores: boundary_debug.clone(),
                stitched_phase: stitched.clone(),
                final_stitched_calls: final_stitched_calls.clone(),
            };

            let out_dir = ensure_results_dir(&cfg.sample_id, ancestry_off)?;
            let local_json = out_dir.join("local_decode_result.json");
            let summary_json = out_dir.join("benchmark_summary.json");
            let candidates_json = out_dir.join("candidate_set.json");
            let calls_tsv = out_dir.join("local_calls.tsv");
            let boundaries_json = out_dir.join("boundary_decisions.json");
            let block_summary_json = out_dir.join("derived_block_summaries.json");
            let boundary_debug_json = out_dir.join("derived_boundary_scores.json");
            let run_report_json = out_dir.join("run_report.json");
            let final_calls_json = out_dir.join("final_stitched_calls.json");
            let final_calls_tsv = out_dir.join("final_stitched_calls.tsv");
            let local_debug_json = out_dir.join("local_window_states.json");

            fs::write(&candidates_json, serde_json::to_string_pretty(&candidates)?)
                .with_context(|| format!("failed to write {}", candidates_json.display()))?;
            fs::write(&local_json, serde_json::to_string_pretty(&result)?)
                .with_context(|| format!("failed to write {}", local_json.display()))?;
            fs::write(&summary_json, serde_json::to_string_pretty(&summary)?)
                .with_context(|| format!("failed to write {}", summary_json.display()))?;
            fs::write(&boundaries_json, serde_json::to_string_pretty(&boundaries)?)
                .with_context(|| format!("failed to write {}", boundaries_json.display()))?;
            fs::write(
                &block_summary_json,
                serde_json::to_string_pretty(&derived_blocks)?,
            )
            .with_context(|| format!("failed to write {}", block_summary_json.display()))?;
            fs::write(
                &boundary_debug_json,
                serde_json::to_string_pretty(&boundary_debug)?,
            )
            .with_context(|| format!("failed to write {}", boundary_debug_json.display()))?;
            fs::write(&run_report_json, serde_json::to_string_pretty(&report)?)
                .with_context(|| format!("failed to write {}", run_report_json.display()))?;
            fs::write(
                &final_calls_json,
                serde_json::to_string_pretty(&final_stitched_calls)?,
            )
            .with_context(|| format!("failed to write {}", final_calls_json.display()))?;
            fs::write(
                &local_debug_json,
                serde_json::to_string_pretty(&local_debug)?,
            )
            .with_context(|| format!("failed to write {}", local_debug_json.display()))?;

            let mut tsv =
                String::from("pos\tblock_id\tphase_state\tlocal_phase_state\tconfidence\n");
            for call in &result.calls {
                let block_id = block_ranges
                    .iter()
                    .find(|r| call.pos >= r.start && call.pos <= r.end)
                    .map(|r| r.block_id.as_str())
                    .unwrap_or("unassigned");

                tsv.push_str(&format!(
                    "{}\t{}\t{}\t{}\t{:.6}\n",
                    call.pos, block_id, call.phase_state, call.phase_state, call.confidence
                ));
            }
            fs::write(&calls_tsv, tsv)
                .with_context(|| format!("failed to write {}", calls_tsv.display()))?;

            let mut final_tsv = String::from(
                "pos\tblock_id\tlocal_phase_state\tstitched_phase_state\tconfidence\n",
            );
            for call in &final_stitched_calls.calls {
                final_tsv.push_str(&format!(
                    "{}\t{}\t{}\t{}\t{:.6}\n",
                    call.pos,
                    call.block_id,
                    call.local_phase_state,
                    call.stitched_phase_state,
                    call.confidence
                ));
            }
            fs::write(&final_calls_tsv, final_tsv)
                .with_context(|| format!("failed to write {}", final_calls_tsv.display()))?;

            println!("{}", serde_json::to_string_pretty(&summary)?);
        }
    }

    Ok(())
}
