use std::cmp::Ordering;
use std::collections::{BTreeMap, BTreeSet, HashMap, HashSet};
use std::fs;
use std::path::Path;
use std::process::Command;

use anyhow::{bail, Context, Result};
use clap::{Parser, ValueEnum};
use csv::{ReaderBuilder, Terminator, WriterBuilder};
use serde::Serialize;
use serde_json::Value;

#[derive(Copy, Clone, Debug, Eq, PartialEq, ValueEnum)]
enum Mode {
    Pbwt,
    PbwtHmm,
    PbwtV2,
    PbwtHmmV2,
}

#[derive(Copy, Clone, Debug, Eq, PartialEq, ValueEnum)]
enum PbwtRetrieval {
    ExactPrefix,
    NeighborDebug,
}

#[derive(Parser, Debug)]
#[command(about = "PBWT-backed panel projection from an AWPhase read backbone")]
struct Args {
    #[arg(long, value_enum, default_value_t = Mode::Pbwt)]
    mode: Mode,

    #[arg(long)]
    panel_bcf: String,
    #[arg(long)]
    backbone_local_calls_tsv: String,
    #[arg(long)]
    variant_json: String,
    #[arg(long)]
    chrom: String,
    #[arg(long)]
    start: u64,
    #[arg(long)]
    end: u64,

    #[arg(long)]
    genetic_map: Option<String>,

    #[arg(long)]
    out_tsv: String,
    #[arg(long)]
    out_candidates_tsv: String,
    #[arg(long)]
    out_summary_json: String,

    #[arg(long, default_value_t = 20_000)]
    max_anchor_dist_bp: u64,
    #[arg(long, default_value_t = 12)]
    max_anchors_each_side: usize,
    #[arg(long, default_value_t = 64)]
    top_k: usize,
    #[arg(long, default_value_t = 24)]
    pbwt_seed_count: usize,
    #[arg(long, default_value_t = 24)]
    pbwt_neighbor_radius: usize,
    #[arg(long, value_enum, default_value_t = PbwtRetrieval::ExactPrefix)]
    pbwt_retrieval: PbwtRetrieval,
    #[arg(long, default_value_t = 1)]
    min_pbwt_match_anchors: usize,
    #[arg(long, default_value_t = 2)]
    min_anchors: usize,
    #[arg(long, default_value_t = 12)]
    min_donors: usize,

    #[arg(long, default_value_t = 0.82)]
    min_confidence: f64,
    #[arg(long, default_value_t = 1.0)]
    min_margin: f64,
    #[arg(long, default_value_t = 1.0)]
    hmm_min_margin: f64,
    #[arg(long, default_value_t = 0.25)]
    min_best_vs_second_margin: f64,

    #[arg(long, default_value_t = 4.0)]
    hmm_mismatch_penalty: f64,
    #[arg(long, default_value_t = 1.0)]
    hmm_recomb_scale: f64,
}

#[derive(Debug, Clone)]
struct Variant {
    pos: u64,
    ref_allele: String,
    alt_allele: String,
}

#[derive(Debug, Clone)]
struct CallRow {
    fields: BTreeMap<String, String>,
}

#[derive(Debug, Clone)]
struct Anchor {
    pos: u64,
    state: i8,
    block_id: String,
}

#[derive(Debug, Clone)]
struct PanelSite {
    alleles: Vec<i8>,
}

#[derive(Debug)]
struct PanelData {
    n_haps: usize,
    sites: HashMap<u64, PanelSite>,
    positions_found: usize,
    allele_mismatch: usize,
    nonbiallelic: usize,
}

#[derive(Debug, Clone)]
struct PbwtSiteOrder {
    order: Vec<usize>,
    rank: Vec<usize>,
    divergence: Vec<usize>,
}

#[derive(Debug)]
struct PbwtIndex {
    by_pos: HashMap<u64, PbwtSiteOrder>,
    site_index_by_pos: HashMap<u64, usize>,
}

#[derive(Debug, Clone)]
struct Prediction {
    pos: u64,
    block_id: String,
    pred_state: i8,
    accepted: bool,
    support_plus: f64,
    support_minus: f64,
    margin: f64,
    confidence: f64,
    donors: usize,
    anchors: usize,
    best_vs_second_margin: f64,
    n_block_candidates: usize,
    method: String,
    donor_entropy: f64,
    effective_donors: f64,
    pbwt_left_match_anchors: usize,
    pbwt_right_match_anchors: usize,
    retrieval: String,
}

#[derive(Debug, Clone)]
struct BlockCandidate {
    pos: u64,
    block_id: String,
    anchors: Vec<Anchor>,
}

#[derive(Debug, Clone)]
struct ScoredDonor {
    hap: usize,
    score: f64,
    matched_anchors: usize,
    left_match_anchors: usize,
    right_match_anchors: usize,
}

#[derive(Debug)]
struct GeneticMap {
    bp: Vec<u64>,
    cm: Vec<f64>,
}

#[derive(Serialize)]
struct Summary<'a> {
    mode: &'a str,
    panel_bcf: &'a str,
    backbone_local_calls_tsv: &'a str,
    variant_json: &'a str,
    chrom: &'a str,
    start: u64,
    end: u64,
    backbone_anchors: usize,
    candidate_unphased_positions: usize,
    candidate_block_pairs: usize,
    panel_positions_needed: usize,
    panel_positions_found: usize,
    panel_positions_missing: usize,
    panel_allele_mismatch_records: usize,
    panel_nonbiallelic_records: usize,
    pbwt_index_positions: usize,
    pbwt_retrieval: &'a str,
    hmm_decoder: &'a str,
    accepted_sites: usize,
    out_tsv: &'a str,
}

fn as_array_root(value: &Value) -> Result<&Vec<Value>> {
    if let Some(arr) = value.as_array() {
        return Ok(arr);
    }
    if let Some(arr) = value.get("variants").and_then(Value::as_array) {
        return Ok(arr);
    }
    bail!("variant JSON must be an array or an object with a variants array")
}

fn load_variants(path: &str, chrom: &str, start: u64, end: u64) -> Result<Vec<Variant>> {
    let text = fs::read_to_string(path)
        .with_context(|| format!("failed to read variant JSON: {path}"))?;
    let value: Value = serde_json::from_str(&text)
        .with_context(|| format!("failed to parse variant JSON: {path}"))?;

    let mut out = Vec::new();
    for row in as_array_root(&value)? {
        let Some(pos) = row
            .get("pos")
            .and_then(Value::as_u64)
            .or_else(|| row.get("pos").and_then(Value::as_str)?.parse().ok())
        else {
            continue;
        };

        if pos < start || pos > end {
            continue;
        }

        let row_chrom = row
            .get("chrom")
            .or_else(|| row.get("contig"))
            .and_then(Value::as_str)
            .unwrap_or(chrom);
        if !row_chrom.is_empty() && row_chrom != chrom {
            continue;
        }

        let ref_allele = row
            .get("ref_allele")
            .or_else(|| row.get("ref"))
            .and_then(Value::as_str)
            .unwrap_or("")
            .to_ascii_uppercase();
        let alt_allele = row
            .get("alt_allele")
            .or_else(|| row.get("alt"))
            .and_then(Value::as_str)
            .unwrap_or("")
            .to_ascii_uppercase();

        if ref_allele.is_empty() || alt_allele.is_empty() {
            continue;
        }

        out.push(Variant {
            pos,
            ref_allele,
            alt_allele,
        });
    }
    out.sort_by_key(|v| v.pos);
    out.dedup_by_key(|v| v.pos);
    Ok(out)
}

fn read_calls(path: &str) -> Result<(Vec<String>, Vec<CallRow>)> {
    let mut reader = ReaderBuilder::new()
        .delimiter(b'\t')
        .from_path(path)
        .with_context(|| format!("failed to read calls TSV: {path}"))?;
    let headers: Vec<String> = reader.headers()?.iter().map(str::to_string).collect();
    let mut rows = Vec::new();

    for record in reader.records() {
        let record = record?;
        let mut fields = BTreeMap::new();
        for (key, value) in headers.iter().zip(record.iter()) {
            fields.insert(key.clone(), value.to_string());
        }
        rows.push(CallRow { fields });
    }

    Ok((headers, rows))
}

fn parse_i8(value: Option<&String>) -> i8 {
    value
        .and_then(|v| v.parse::<f64>().ok())
        .map(|v| v as i8)
        .unwrap_or(0)
}

fn parse_u64(value: Option<&String>) -> Option<u64> {
    value
        .and_then(|v| v.parse::<f64>().ok())
        .map(|v| v as u64)
}

fn call_state(row: &CallRow) -> i8 {
    if row.fields.contains_key("local_phase_state") {
        parse_i8(row.fields.get("local_phase_state"))
    } else {
        parse_i8(row.fields.get("phase_state"))
    }
}

fn set_call_state(row: &mut CallRow, state: i8) {
    if row.fields.contains_key("local_phase_state") {
        row.fields
            .insert("local_phase_state".to_string(), state.to_string());
    }
    if row.fields.contains_key("phase_state") {
        row.fields
            .insert("phase_state".to_string(), state.to_string());
    }
}

fn clean_block_id(value: Option<&String>) -> Option<String> {
    let b = value?.trim();
    if b.is_empty() || matches!(b, "0" | "." | "NA" | "na" | "none" | "unassigned") {
        None
    } else {
        Some(b.to_string())
    }
}

fn anchor_expected_allele(state: i8, plus_path: bool) -> i8 {
    if plus_path {
        if state == 1 {
            1
        } else {
            0
        }
    } else if state == 1 {
        0
    } else {
        1
    }
}

fn candidate_positions(rows: &[CallRow], variant_pos: &HashSet<u64>) -> Vec<u64> {
    rows.iter()
        .filter_map(|row| {
            let pos = parse_u64(row.fields.get("pos"))?;
            if variant_pos.contains(&pos) && call_state(row) == 0 {
                Some(pos)
            } else {
                None
            }
        })
        .collect()
}

fn load_anchors(rows: &[CallRow], variant_pos: &HashSet<u64>) -> Vec<Anchor> {
    let mut anchors = Vec::new();
    for row in rows {
        let Some(pos) = parse_u64(row.fields.get("pos")) else {
            continue;
        };
        if !variant_pos.contains(&pos) {
            continue;
        }
        let state = call_state(row);
        if state == 0 {
            continue;
        }
        let Some(block_id) = clean_block_id(row.fields.get("block_id")) else {
            continue;
        };
        anchors.push(Anchor {
            pos,
            state,
            block_id,
        });
    }
    anchors.sort_by_key(|a| a.pos);
    anchors
}

fn nearby_anchors(
    pos: u64,
    anchors: &[Anchor],
    max_dist: u64,
    max_each_side: usize,
) -> Vec<Anchor> {
    let idx = anchors.partition_point(|a| a.pos < pos);
    let mut out = Vec::new();

    let mut count = 0;
    let mut j = idx;
    while j > 0 && count < max_each_side {
        j -= 1;
        let dist = pos.saturating_sub(anchors[j].pos);
        if dist > max_dist {
            break;
        }
        out.push(anchors[j].clone());
        count += 1;
    }

    out.reverse();

    let mut count = 0;
    let mut j = idx;
    while j < anchors.len() && count < max_each_side {
        let dist = anchors[j].pos.saturating_sub(pos);
        if dist > max_dist {
            break;
        }
        out.push(anchors[j].clone());
        count += 1;
        j += 1;
    }

    out
}

fn block_candidates_for_position(
    pos: u64,
    anchors: &[Anchor],
    max_dist: u64,
    max_each_side: usize,
    min_anchors: usize,
) -> Vec<BlockCandidate> {
    let mut by_block: HashMap<String, Vec<Anchor>> = HashMap::new();
    for anchor in nearby_anchors(pos, anchors, max_dist, max_each_side) {
        by_block
            .entry(anchor.block_id.clone())
            .or_default()
            .push(anchor);
    }

    let mut out: Vec<_> = by_block
        .into_iter()
        .filter_map(|(block_id, mut block_anchors)| {
            block_anchors.sort_by_key(|a| a.pos);
            if block_anchors.len() >= min_anchors {
                Some(BlockCandidate {
                    pos,
                    block_id,
                    anchors: block_anchors,
                })
            } else {
                None
            }
        })
        .collect();
    out.sort_by(|a, b| a.block_id.cmp(&b.block_id));
    out
}

fn parse_gt(gt: &str) -> (i8, i8) {
    if !gt.contains('|') {
        return (-1, -1);
    }
    let mut parts = gt.split('|');
    let a = parse_panel_allele(parts.next().unwrap_or(""));
    let b = parse_panel_allele(parts.next().unwrap_or(""));
    (a, b)
}

fn parse_panel_allele(value: &str) -> i8 {
    match value {
        "0" => 0,
        "1" => 1,
        _ => -1,
    }
}

fn load_panel_bcf(
    panel_bcf: &str,
    chrom: &str,
    start: u64,
    end: u64,
    needed_variants: &HashMap<u64, Variant>,
) -> Result<PanelData> {
    let region = format!("{chrom}:{start}-{end}");
    let output = Command::new("bcftools")
        .args([
            "query",
            "-r",
            &region,
            "-f",
            "%POS\t%REF\t%ALT[\t%GT]\n",
            panel_bcf,
        ])
        .output()
        .with_context(|| "failed to execute bcftools query")?;

    if !output.status.success() {
        bail!(
            "bcftools query failed: {}",
            String::from_utf8_lossy(&output.stderr)
        );
    }

    let stdout = String::from_utf8(output.stdout).context("bcftools output was not UTF-8")?;
    let mut sites = HashMap::new();
    let mut n_haps = 0;
    let mut positions_found = 0;
    let mut allele_mismatch = 0;
    let mut nonbiallelic = 0;

    for line in stdout.lines() {
        if line.trim().is_empty() {
            continue;
        }
        let parts: Vec<&str> = line.split('\t').collect();
        if parts.len() < 4 {
            continue;
        }
        let Ok(pos) = parts[0].parse::<u64>() else {
            continue;
        };
        let Some(variant) = needed_variants.get(&pos) else {
            continue;
        };
        if parts[2].contains(',') {
            nonbiallelic += 1;
            continue;
        }
        if parts[1].to_ascii_uppercase() != variant.ref_allele
            || parts[2].to_ascii_uppercase() != variant.alt_allele
        {
            allele_mismatch += 1;
            continue;
        }

        let sample_count = parts.len() - 3;
        if n_haps == 0 {
            n_haps = sample_count * 2;
        }
        if sample_count * 2 != n_haps {
            continue;
        }

        let mut alleles = Vec::with_capacity(n_haps);
        for gt in &parts[3..] {
            let (a, b) = parse_gt(gt);
            alleles.push(a);
            alleles.push(b);
        }

        sites.insert(pos, PanelSite { alleles });
        positions_found += 1;
    }

    Ok(PanelData {
        n_haps,
        sites,
        positions_found,
        allele_mismatch,
        nonbiallelic,
    })
}

fn build_pbwt_index(panel: &PanelData, positions: &[u64]) -> PbwtIndex {
    let mut order: Vec<usize> = (0..panel.n_haps).collect();
    let mut divergence = vec![0_usize; panel.n_haps];
    let mut by_pos = HashMap::new();
    let mut site_index_by_pos = HashMap::new();

    for (site_index, pos) in positions.iter().enumerate() {
        let Some(site) = panel.sites.get(pos) else {
            continue;
        };

        let mut zeros = Vec::new();
        let mut ones = Vec::new();
        let mut missing = Vec::new();
        let mut zero_divergence = Vec::new();
        let mut one_divergence = Vec::new();
        let mut missing_divergence = Vec::new();
        let mut zero_start = site_index + 1;
        let mut one_start = site_index + 1;
        let mut missing_start = site_index + 1;

        for (idx, hap) in order.iter().enumerate() {
            let div = divergence[idx];
            match site.alleles[*hap] {
                0 => {
                    zeros.push(*hap);
                    zero_divergence.push(zero_start);
                    zero_start = div;
                }
                1 => {
                    ones.push(*hap);
                    one_divergence.push(one_start);
                    one_start = div;
                }
                _ => {
                    missing.push(*hap);
                    missing_divergence.push(missing_start);
                    missing_start = div;
                }
            }
        }

        order.clear();
        order.extend(zeros);
        order.extend(ones);
        order.extend(missing);
        divergence.clear();
        divergence.extend(zero_divergence);
        divergence.extend(one_divergence);
        divergence.extend(missing_divergence);

        let mut rank = vec![0_usize; panel.n_haps];
        for (idx, hap) in order.iter().enumerate() {
            rank[*hap] = idx;
        }

        by_pos.insert(
            *pos,
            PbwtSiteOrder {
                order: order.clone(),
                rank,
                divergence: divergence.clone(),
            },
        );
        site_index_by_pos.insert(*pos, site_index);
    }

    PbwtIndex {
        by_pos,
        site_index_by_pos,
    }
}

fn anchor_weight(candidate_pos: u64, anchor_pos: u64) -> f64 {
    let dist = candidate_pos.abs_diff(anchor_pos) as f64;
    1.0 / (1.0 + dist / 5_000.0)
}

fn score_haplotypes(
    candidate_pos: u64,
    anchors: &[Anchor],
    plus_path: bool,
    panel: &PanelData,
) -> Vec<f64> {
    let mut scores = vec![0.0_f64; panel.n_haps];
    for anchor in anchors {
        let Some(site) = panel.sites.get(&anchor.pos) else {
            continue;
        };
        let expected = anchor_expected_allele(anchor.state, plus_path);
        let weight = anchor_weight(candidate_pos, anchor.pos);
        for (hap, allele) in site.alleles.iter().enumerate() {
            if *allele == expected {
                scores[hap] += weight;
            } else if *allele == 0 || *allele == 1 {
                scores[hap] -= 0.25 * weight;
            }
        }
    }
    scores
}

fn sort_haps_by_score(scores: &[f64]) -> Vec<ScoredDonor> {
    let mut donors: Vec<_> = scores
        .iter()
        .enumerate()
        .map(|(hap, score)| ScoredDonor {
            hap,
            score: *score,
            matched_anchors: 0,
            left_match_anchors: 0,
            right_match_anchors: 0,
        })
        .collect();
    donors.sort_by(|a, b| {
        b.score
            .partial_cmp(&a.score)
            .unwrap_or(Ordering::Equal)
            .then_with(|| a.hap.cmp(&b.hap))
    });
    donors
}

fn pbwt_retrieve(
    scores: &[f64],
    anchors: &[Anchor],
    pbwt: &PbwtIndex,
    top_k: usize,
    seed_count: usize,
    radius: usize,
) -> Vec<ScoredDonor> {
    let ranked = sort_haps_by_score(scores);
    let mut candidate_haps = HashSet::new();

    for seed in ranked
        .iter()
        .filter(|d| d.score > 0.0)
        .take(seed_count.max(top_k))
    {
        candidate_haps.insert(seed.hap);
        for anchor in anchors {
            let Some(order) = pbwt.by_pos.get(&anchor.pos) else {
                continue;
            };
            let rank = order.rank[seed.hap];
            let lo = rank.saturating_sub(radius);
            let hi = (rank + radius + 1).min(order.order.len());
            for hap in &order.order[lo..hi] {
                candidate_haps.insert(*hap);
            }
        }
    }

    if candidate_haps.is_empty() {
        return ranked.into_iter().take(top_k).collect();
    }

    let mut donors: Vec<_> = candidate_haps
        .into_iter()
        .map(|hap| ScoredDonor {
            hap,
            score: scores[hap],
            matched_anchors: 0,
            left_match_anchors: 0,
            right_match_anchors: 0,
        })
        .collect();
    donors.sort_by(|a, b| {
        b.score
            .partial_cmp(&a.score)
            .unwrap_or(Ordering::Equal)
            .then_with(|| a.hap.cmp(&b.hap))
    });
    donors.truncate(top_k);
    donors
}

fn path_match_stats(
    candidate_pos: u64,
    anchors: &[Anchor],
    plus_path: bool,
    panel: &PanelData,
    hap: usize,
) -> ScoredDonor {
    let mut score = 0.0;
    let mut matched_anchors = 0;

    for anchor in anchors {
        let Some(site) = panel.sites.get(&anchor.pos) else {
            continue;
        };
        let expected = anchor_expected_allele(anchor.state, plus_path);
        let weight = anchor_weight(candidate_pos, anchor.pos);
        match site.alleles[hap] {
            allele if allele == expected => {
                score += weight;
                matched_anchors += 1;
            }
            0 | 1 => {
                score -= 0.25 * weight;
            }
            _ => {}
        }
    }

    let mut left_match_anchors = 0;
    for anchor in anchors.iter().rev().filter(|a| a.pos <= candidate_pos) {
        let Some(site) = panel.sites.get(&anchor.pos) else {
            continue;
        };
        let expected = anchor_expected_allele(anchor.state, plus_path);
        match site.alleles[hap] {
            allele if allele == expected => left_match_anchors += 1,
            0 | 1 => break,
            _ => continue,
        }
    }

    let mut right_match_anchors = 0;
    for anchor in anchors.iter().filter(|a| a.pos > candidate_pos) {
        let Some(site) = panel.sites.get(&anchor.pos) else {
            continue;
        };
        let expected = anchor_expected_allele(anchor.state, plus_path);
        match site.alleles[hap] {
            allele if allele == expected => right_match_anchors += 1,
            0 | 1 => break,
            _ => continue,
        }
    }

    let exact_bonus = 0.75 * (left_match_anchors + right_match_anchors) as f64
        + 0.10 * matched_anchors as f64;

    ScoredDonor {
        hap,
        score: score + exact_bonus,
        matched_anchors,
        left_match_anchors,
        right_match_anchors,
    }
}

fn pbwt_long_match_neighbors(
    seed: &ScoredDonor,
    candidate_pos: u64,
    anchors: &[Anchor],
    pbwt: &PbwtIndex,
    min_match_anchors: usize,
) -> Vec<usize> {
    let Some(anchor) = anchors
        .iter()
        .rev()
        .find(|a| a.pos <= candidate_pos && pbwt.by_pos.contains_key(&a.pos))
    else {
        return vec![seed.hap];
    };
    let Some(order) = pbwt.by_pos.get(&anchor.pos) else {
        return vec![seed.hap];
    };
    let Some(anchor_site_index) = pbwt.site_index_by_pos.get(&anchor.pos).copied() else {
        return vec![seed.hap];
    };

    let shared_anchor_count = seed.left_match_anchors.max(min_match_anchors);
    let min_site_index = anchor_site_index.saturating_sub(shared_anchor_count.saturating_sub(1));
    let rank = order.rank[seed.hap];
    let mut out = vec![seed.hap];

    let mut j = rank;
    while j > 0 {
        if order.divergence[j] > min_site_index {
            break;
        }
        j -= 1;
        out.push(order.order[j]);
    }

    let mut j = rank + 1;
    while j < order.order.len() {
        if order.divergence[j] > min_site_index {
            break;
        }
        out.push(order.order[j]);
        j += 1;
    }

    out
}

fn exact_prefix_retrieve(
    candidate_pos: u64,
    anchors: &[Anchor],
    plus_path: bool,
    panel: &PanelData,
    pbwt: &PbwtIndex,
    top_k: usize,
    seed_count: usize,
    min_match_anchors: usize,
) -> Vec<ScoredDonor> {
    let mut ranked: Vec<_> = (0..panel.n_haps)
        .map(|hap| path_match_stats(candidate_pos, anchors, plus_path, panel, hap))
        .collect();
    ranked.sort_by(|a, b| {
        (b.left_match_anchors + b.right_match_anchors)
            .cmp(&(a.left_match_anchors + a.right_match_anchors))
            .then_with(|| b.matched_anchors.cmp(&a.matched_anchors))
            .then_with(|| b.score.partial_cmp(&a.score).unwrap_or(Ordering::Equal))
            .then_with(|| a.hap.cmp(&b.hap))
    });

    let mut candidate_haps = HashSet::new();
    for seed in ranked
        .iter()
        .filter(|d| d.score > 0.0)
        .filter(|d| d.left_match_anchors + d.right_match_anchors >= min_match_anchors)
        .take(seed_count.max(top_k))
    {
        for hap in pbwt_long_match_neighbors(seed, candidate_pos, anchors, pbwt, min_match_anchors)
        {
            candidate_haps.insert(hap);
        }
    }

    if candidate_haps.is_empty() {
        return ranked.into_iter().filter(|d| d.score > 0.0).take(top_k).collect();
    }

    let mut donors: Vec<_> = candidate_haps
        .into_iter()
        .map(|hap| path_match_stats(candidate_pos, anchors, plus_path, panel, hap))
        .filter(|d| d.score > 0.0)
        .collect();
    donors.sort_by(|a, b| {
        (b.left_match_anchors + b.right_match_anchors)
            .cmp(&(a.left_match_anchors + a.right_match_anchors))
            .then_with(|| b.matched_anchors.cmp(&a.matched_anchors))
            .then_with(|| b.score.partial_cmp(&a.score).unwrap_or(Ordering::Equal))
            .then_with(|| a.hap.cmp(&b.hap))
    });
    donors.truncate(top_k);
    donors
}

fn retrieve_path_donors(
    candidate_pos: u64,
    anchors: &[Anchor],
    plus_path: bool,
    panel: &PanelData,
    pbwt: &PbwtIndex,
    args: &Args,
    retrieval: PbwtRetrieval,
) -> Vec<ScoredDonor> {
    match retrieval {
        PbwtRetrieval::NeighborDebug => {
            let scores = score_haplotypes(candidate_pos, anchors, plus_path, panel);
            pbwt_retrieve(
                &scores,
                anchors,
                pbwt,
                args.top_k,
                args.pbwt_seed_count,
                args.pbwt_neighbor_radius,
            )
        }
        PbwtRetrieval::ExactPrefix => exact_prefix_retrieve(
            candidate_pos,
            anchors,
            plus_path,
            panel,
            pbwt,
            args.top_k,
            args.pbwt_seed_count,
            args.min_pbwt_match_anchors,
        ),
    }
}

fn donor_entropy(donors: &[ScoredDonor]) -> f64 {
    let weights: Vec<f64> = donors.iter().map(|d| d.score.max(0.0)).collect();
    entropy_from_weights(&weights)
}

fn effective_donor_count(donors: &[ScoredDonor]) -> f64 {
    donor_entropy(donors).exp()
}

fn max_left_match(donors: &[ScoredDonor]) -> usize {
    donors.iter().map(|d| d.left_match_anchors).max().unwrap_or(0)
}

fn max_right_match(donors: &[ScoredDonor]) -> usize {
    donors.iter().map(|d| d.right_match_anchors).max().unwrap_or(0)
}

fn pbwt_prediction_for_block(
    bc: &BlockCandidate,
    panel: &PanelData,
    pbwt: &PbwtIndex,
    args: &Args,
    method: &str,
    retrieval: PbwtRetrieval,
) -> Option<Prediction> {
    let candidate_site = panel.sites.get(&bc.pos)?;
    let plus_donors = retrieve_path_donors(
        bc.pos,
        &bc.anchors,
        true,
        panel,
        pbwt,
        args,
        retrieval,
    );
    let minus_donors = retrieve_path_donors(
        bc.pos,
        &bc.anchors,
        false,
        panel,
        pbwt,
        args,
        retrieval,
    );

    let mut support_plus = 0.0;
    let mut support_minus = 0.0;
    let mut donors_used = 0;

    for donor in plus_donors.iter().filter(|d| d.score > 0.0) {
        let allele = candidate_site.alleles[donor.hap];
        if allele == 1 {
            support_plus += donor.score;
            donors_used += 1;
        } else if allele == 0 {
            support_minus += donor.score;
            donors_used += 1;
        }
    }

    for donor in minus_donors.iter().filter(|d| d.score > 0.0) {
        let allele = candidate_site.alleles[donor.hap];
        if allele == 0 {
            support_plus += donor.score;
            donors_used += 1;
        } else if allele == 1 {
            support_minus += donor.score;
            donors_used += 1;
        }
    }

    let total = support_plus + support_minus;
    if total <= 0.0 {
        return None;
    }

    let confidence = support_plus.max(support_minus) / total;
    let margin = (support_plus - support_minus).abs();
    let pred_state = if support_plus >= support_minus { 1 } else { -1 };

    Some(Prediction {
        pos: bc.pos,
        block_id: bc.block_id.clone(),
        pred_state,
        accepted: false,
        support_plus,
        support_minus,
        margin,
        confidence,
        donors: donors_used,
        anchors: bc.anchors.len(),
        best_vs_second_margin: 0.0,
        n_block_candidates: 0,
        method: method.to_string(),
        donor_entropy: (donor_entropy(&plus_donors) + donor_entropy(&minus_donors)) / 2.0,
        effective_donors: effective_donor_count(&plus_donors) + effective_donor_count(&minus_donors),
        pbwt_left_match_anchors: max_left_match(&plus_donors).max(max_left_match(&minus_donors)),
        pbwt_right_match_anchors: max_right_match(&plus_donors).max(max_right_match(&minus_donors)),
        retrieval: format!("{:?}", retrieval),
    })
}

fn choose_predictions(
    candidates_by_pos: &HashMap<u64, Vec<Prediction>>,
    args: &Args,
) -> HashMap<u64, Prediction> {
    let mut chosen = HashMap::new();

    for (pos, preds) in candidates_by_pos {
        if preds.is_empty() {
            continue;
        }
        let mut ranked = preds.clone();
        ranked.sort_by(|a, b| {
            b.margin
                .partial_cmp(&a.margin)
                .unwrap_or(Ordering::Equal)
                .then_with(|| {
                    b.confidence
                        .partial_cmp(&a.confidence)
                        .unwrap_or(Ordering::Equal)
                })
                .then_with(|| b.donors.cmp(&a.donors))
                .then_with(|| b.anchors.cmp(&a.anchors))
                .then_with(|| a.block_id.cmp(&b.block_id))
        });

        let second_margin = ranked.get(1).map(|p| p.margin).unwrap_or(0.0);
        let best_vs_second = ranked[0].margin - second_margin;
        let mut best = ranked[0].clone();
        best.best_vs_second_margin = best_vs_second;
        best.n_block_candidates = ranked.len();
        best.accepted = best.donors >= args.min_donors
            && best.anchors >= args.min_anchors
            && best.confidence >= args.min_confidence
            && best.margin
                >= if best.method.contains("hmm") {
                    args.hmm_min_margin
                } else {
                    args.min_margin
                }
            && best_vs_second >= args.min_best_vs_second_margin;
        if best.accepted {
            chosen.insert(*pos, best);
        }
    }

    chosen
}

fn pbwt_predictions(
    block_candidates: &[BlockCandidate],
    panel: &PanelData,
    pbwt: &PbwtIndex,
    args: &Args,
    method: &str,
    retrieval: PbwtRetrieval,
) -> HashMap<u64, Vec<Prediction>> {
    let mut by_pos: HashMap<u64, Vec<Prediction>> = HashMap::new();
    for bc in block_candidates {
        if let Some(pred) = pbwt_prediction_for_block(bc, panel, pbwt, args, method, retrieval) {
            by_pos.entry(bc.pos).or_default().push(pred);
        }
    }
    by_pos
}

fn load_genetic_map(path: Option<&str>) -> Result<Option<GeneticMap>> {
    let Some(path) = path else {
        return Ok(None);
    };
    let text = fs::read_to_string(path)
        .with_context(|| format!("failed to read genetic map: {path}"))?;
    let mut pairs = Vec::new();

    for line in text.lines() {
        if line.trim().is_empty() || line.starts_with('#') {
            continue;
        }
        let parts: Vec<&str> = line.split_whitespace().collect();
        if parts.len() < 4 {
            continue;
        }
        let Ok(cm) = parts[parts.len() - 2].parse::<f64>() else {
            continue;
        };
        let Ok(bp) = parts[parts.len() - 1].parse::<f64>() else {
            continue;
        };
        pairs.push((bp as u64, cm));
    }

    if pairs.is_empty() {
        bail!("no genetic map rows parsed from {path}");
    }

    pairs.sort_by_key(|x| x.0);
    Ok(Some(GeneticMap {
        bp: pairs.iter().map(|x| x.0).collect(),
        cm: pairs.iter().map(|x| x.1).collect(),
    }))
}

fn interp_cm(pos: u64, genetic_map: Option<&GeneticMap>) -> f64 {
    let Some(map) = genetic_map else {
        return pos as f64 / 1_000_000.0;
    };
    if pos <= map.bp[0] {
        return map.cm[0];
    }
    if pos >= *map.bp.last().unwrap() {
        return *map.cm.last().unwrap();
    }
    let hi = map.bp.partition_point(|bp| *bp < pos);
    let lo = hi - 1;
    let bp0 = map.bp[lo] as f64;
    let bp1 = map.bp[hi] as f64;
    let cm0 = map.cm[lo];
    let cm1 = map.cm[hi];
    if bp1 <= bp0 {
        cm0
    } else {
        cm0 + (cm1 - cm0) * ((pos as f64 - bp0) / (bp1 - bp0))
    }
}

fn logsumexp(values: &[f64]) -> f64 {
    let max_value = values.iter().copied().fold(f64::NEG_INFINITY, f64::max);
    if !max_value.is_finite() {
        return max_value;
    }
    let sum: f64 = values.iter().map(|v| (*v - max_value).exp()).sum();
    max_value + sum.ln()
}

fn normalize_log_weights(logw: &[f64]) -> Vec<f64> {
    let z = logsumexp(logw);
    if !z.is_finite() {
        return vec![0.0; logw.len()];
    }
    logw.iter().map(|v| (*v - z).exp()).collect()
}

fn entropy_from_weights(weights: &[f64]) -> f64 {
    let total: f64 = weights.iter().filter(|v| **v > 0.0).sum();
    if total <= 0.0 {
        return 0.0;
    }
    weights
        .iter()
        .filter(|v| **v > 0.0)
        .map(|v| {
            let p = *v / total;
            -p * p.ln()
        })
        .sum()
}

fn effective_count_from_weights(weights: &[f64]) -> f64 {
    entropy_from_weights(weights).exp()
}

fn hmm_transition_log(
    same_state: bool,
    prev_pos: u64,
    curr_pos: u64,
    genetic_map: Option<&GeneticMap>,
    recomb_scale: f64,
    k: usize,
) -> f64 {
    if k <= 1 {
        return 0.0;
    }
    let d_cm = (interp_cm(curr_pos, genetic_map) - interp_cm(prev_pos, genetic_map)).abs();
    let r = (0.01 * d_cm * recomb_scale).clamp(1e-8, 0.20);
    if same_state {
        (1.0 - r).ln()
    } else {
        (r / ((k - 1) as f64)).ln()
    }
}

fn donor_emission(
    hap: usize,
    pos: u64,
    anchor_state: Option<i8>,
    plus_path: bool,
    panel: &PanelData,
    mismatch_penalty: f64,
) -> f64 {
    let Some(state) = anchor_state else {
        return 0.0;
    };
    let Some(site) = panel.sites.get(&pos) else {
        return -0.5 * mismatch_penalty;
    };
    let expected = anchor_expected_allele(state, plus_path);
    match site.alleles[hap] {
        allele if allele == expected => 0.0,
        0 | 1 => -mismatch_penalty,
        _ => -0.5 * mismatch_penalty,
    }
}

fn donor_emission_log_prob(
    hap: usize,
    pos: u64,
    anchor_state: Option<i8>,
    plus_path: bool,
    panel: &PanelData,
    mismatch_penalty: f64,
) -> f64 {
    let Some(state) = anchor_state else {
        return 0.0;
    };
    let Some(site) = panel.sites.get(&pos) else {
        return (0.5_f64).ln();
    };
    let expected = anchor_expected_allele(state, plus_path);
    let mismatch_prob = (1.0 / (1.0 + mismatch_penalty.exp())).clamp(1e-5, 0.45);
    let match_prob = 1.0 - mismatch_prob;
    match site.alleles[hap] {
        allele if allele == expected => match_prob.ln(),
        0 | 1 => mismatch_prob.ln(),
        _ => (0.5_f64).ln(),
    }
}

struct HmmPosterior {
    by_pos: HashMap<u64, Vec<f64>>,
}

fn run_forward_backward_hmm(
    donors: &[ScoredDonor],
    positions: &[u64],
    anchor_state_by_pos: &HashMap<u64, i8>,
    plus_path: bool,
    panel: &PanelData,
    genetic_map: Option<&GeneticMap>,
    args: &Args,
) -> HmmPosterior {
    if donors.is_empty() || positions.is_empty() {
        return HmmPosterior {
            by_pos: HashMap::new(),
        };
    }

    let k = donors.len();
    let n = positions.len();
    let mut emissions = vec![vec![0.0_f64; k]; n];
    for (t, pos) in positions.iter().enumerate() {
        let anchor_state = anchor_state_by_pos.get(pos).copied();
        for (j, donor) in donors.iter().enumerate() {
            emissions[t][j] = donor_emission_log_prob(
                donor.hap,
                *pos,
                anchor_state,
                plus_path,
                panel,
                args.hmm_mismatch_penalty,
            );
        }
    }

    let mut alpha = vec![vec![f64::NEG_INFINITY; k]; n];
    let init_logw: Vec<f64> = donors
        .iter()
        .map(|d| d.score.max(0.0) + 1.0)
        .map(f64::ln)
        .collect();
    let init_z = logsumexp(&init_logw);
    for j in 0..k {
        alpha[0][j] = init_logw[j] - init_z + emissions[0][j];
    }
    let z0 = logsumexp(&alpha[0]);
    for v in &mut alpha[0] {
        *v -= z0;
    }

    for t in 1..n {
        for j in 0..k {
            let mut vals = Vec::with_capacity(k);
            for i in 0..k {
                vals.push(
                    alpha[t - 1][i]
                        + hmm_transition_log(
                            i == j,
                            positions[t - 1],
                            positions[t],
                            genetic_map,
                            args.hmm_recomb_scale,
                            k,
                        ),
                );
            }
            alpha[t][j] = emissions[t][j] + logsumexp(&vals);
        }
        let z = logsumexp(&alpha[t]);
        for v in &mut alpha[t] {
            *v -= z;
        }
    }

    let mut beta = vec![vec![0.0_f64; k]; n];
    for t in (0..n.saturating_sub(1)).rev() {
        for i in 0..k {
            let mut vals = Vec::with_capacity(k);
            for j in 0..k {
                vals.push(
                    hmm_transition_log(
                        i == j,
                        positions[t],
                        positions[t + 1],
                        genetic_map,
                        args.hmm_recomb_scale,
                        k,
                    ) + emissions[t + 1][j]
                        + beta[t + 1][j],
                );
            }
            beta[t][i] = logsumexp(&vals);
        }
        let z = logsumexp(&beta[t]);
        for v in &mut beta[t] {
            *v -= z;
        }
    }

    let mut by_pos = HashMap::new();
    for (t, pos) in positions.iter().enumerate() {
        let logw: Vec<f64> = (0..k).map(|j| alpha[t][j] + beta[t][j]).collect();
        by_pos.insert(*pos, normalize_log_weights(&logw));
    }

    HmmPosterior { by_pos }
}

fn run_forward_hmm(
    donors: &[ScoredDonor],
    positions: &[u64],
    anchor_state_by_pos: &HashMap<u64, i8>,
    plus_path: bool,
    panel: &PanelData,
    genetic_map: Option<&GeneticMap>,
    args: &Args,
) -> HashMap<u64, Vec<f64>> {
    if donors.is_empty() || positions.is_empty() {
        return HashMap::new();
    }

    let k = donors.len();
    let mut logw: Vec<f64> = donors
        .iter()
        .map(|d| d.score.max(0.0) + 1.0)
        .map(f64::ln)
        .collect();
    let z = logsumexp(&logw);
    for v in &mut logw {
        *v -= z;
    }

    let mut posteriors = HashMap::new();
    for (pi, pos) in positions.iter().enumerate() {
        if pi > 0 {
            let prev_pos = positions[pi - 1];
            let prev = logw.clone();
            let mut next = vec![f64::NEG_INFINITY; k];

            for j in 0..k {
                let mut vals = Vec::with_capacity(k);
                for (i, prev_score) in prev.iter().enumerate() {
                    vals.push(
                        *prev_score
                            + hmm_transition_log(
                                i == j,
                                prev_pos,
                                *pos,
                                genetic_map,
                                args.hmm_recomb_scale,
                                k,
                            ),
                    );
                }
                next[j] = logsumexp(&vals);
            }
            logw = next;
        }

        let anchor_state = anchor_state_by_pos.get(pos).copied();
        for (j, donor) in donors.iter().enumerate() {
            logw[j] += donor_emission(
                donor.hap,
                *pos,
                anchor_state,
                plus_path,
                panel,
                args.hmm_mismatch_penalty,
            );
        }

        let z = logsumexp(&logw);
        for v in &mut logw {
            *v -= z;
        }
        posteriors.insert(*pos, normalize_log_weights(&logw));
    }

    posteriors
}

fn weighted_donor_entropy(weights: &[f64]) -> f64 {
    entropy_from_weights(weights)
}

fn weighted_effective_donors(weights: &[f64]) -> f64 {
    effective_count_from_weights(weights)
}

fn weighted_left_match(donors: &[ScoredDonor], weights: &[f64]) -> usize {
    donors
        .iter()
        .zip(weights.iter())
        .max_by(|(_, wa), (_, wb)| wa.partial_cmp(wb).unwrap_or(Ordering::Equal))
        .map(|(d, _)| d.left_match_anchors)
        .unwrap_or(0)
}

fn weighted_right_match(donors: &[ScoredDonor], weights: &[f64]) -> usize {
    donors
        .iter()
        .zip(weights.iter())
        .max_by(|(_, wa), (_, wb)| wa.partial_cmp(wb).unwrap_or(Ordering::Equal))
        .map(|(d, _)| d.right_match_anchors)
        .unwrap_or(0)
}

fn hmm_predictions(
    block_candidates: &[BlockCandidate],
    panel: &PanelData,
    pbwt: &PbwtIndex,
    genetic_map: Option<&GeneticMap>,
    args: &Args,
) -> HashMap<u64, Vec<Prediction>> {
    let mut by_block: HashMap<String, Vec<BlockCandidate>> = HashMap::new();
    for bc in block_candidates {
        by_block
            .entry(bc.block_id.clone())
            .or_default()
            .push(bc.clone());
    }

    let mut by_pos: HashMap<u64, Vec<Prediction>> = HashMap::new();
    for (block_id, bcs) in by_block {
        let mut anchor_map: BTreeMap<u64, Anchor> = BTreeMap::new();
        let mut candidate_positions = BTreeSet::new();
        for bc in &bcs {
            candidate_positions.insert(bc.pos);
            for anchor in &bc.anchors {
                anchor_map.insert(anchor.pos, anchor.clone());
            }
        }
        let anchors: Vec<Anchor> = anchor_map.values().cloned().collect();
        if anchors.len() < args.min_anchors {
            continue;
        }

        let center_pos = *candidate_positions
            .iter()
            .nth(candidate_positions.len() / 2)
            .unwrap_or(&anchors[anchors.len() / 2].pos);

        let plus_scores = score_haplotypes(center_pos, &anchors, true, panel);
        let minus_scores = score_haplotypes(center_pos, &anchors, false, panel);
        let plus_donors = pbwt_retrieve(
            &plus_scores,
            &anchors,
            pbwt,
            args.top_k,
            args.pbwt_seed_count,
            args.pbwt_neighbor_radius,
        );
        let minus_donors = pbwt_retrieve(
            &minus_scores,
            &anchors,
            pbwt,
            args.top_k,
            args.pbwt_seed_count,
            args.pbwt_neighbor_radius,
        );

        if plus_donors.len() < args.min_donors || minus_donors.len() < args.min_donors {
            continue;
        }

        let anchor_state_by_pos: HashMap<u64, i8> =
            anchors.iter().map(|a| (a.pos, a.state)).collect();
        let mut sequence_positions: Vec<u64> = anchor_state_by_pos.keys().copied().collect();
        sequence_positions.extend(candidate_positions.iter().copied());
        sequence_positions.sort_unstable();
        sequence_positions.dedup();

        let plus_post = run_forward_hmm(
            &plus_donors,
            &sequence_positions,
            &anchor_state_by_pos,
            true,
            panel,
            genetic_map,
            args,
        );
        let minus_post = run_forward_hmm(
            &minus_donors,
            &sequence_positions,
            &anchor_state_by_pos,
            false,
            panel,
            genetic_map,
            args,
        );

        for pos in candidate_positions {
            let Some(site) = panel.sites.get(&pos) else {
                continue;
            };
            let Some(plus_weights) = plus_post.get(&pos) else {
                continue;
            };
            let Some(minus_weights) = minus_post.get(&pos) else {
                continue;
            };

            let mut plus_alt = 0.0;
            let mut plus_valid = 0.0;
            for (donor, w) in plus_donors.iter().zip(plus_weights) {
                match site.alleles[donor.hap] {
                    1 => {
                        plus_alt += *w;
                        plus_valid += *w;
                    }
                    0 => {
                        plus_valid += *w;
                    }
                    _ => {}
                }
            }

            let mut minus_ref = 0.0;
            let mut minus_valid = 0.0;
            for (donor, w) in minus_donors.iter().zip(minus_weights) {
                match site.alleles[donor.hap] {
                    0 => {
                        minus_ref += *w;
                        minus_valid += *w;
                    }
                    1 => {
                        minus_valid += *w;
                    }
                    _ => {}
                }
            }

            let mut terms = Vec::new();
            if plus_valid > 0.0 {
                terms.push(plus_alt / plus_valid);
            }
            if minus_valid > 0.0 {
                terms.push(minus_ref / minus_valid);
            }
            if terms.is_empty() {
                continue;
            }

            let prob_plus = terms.iter().sum::<f64>() / (terms.len() as f64);
            let support_plus = prob_plus;
            let support_minus = 1.0 - prob_plus;
            let confidence = support_plus.max(support_minus);
            let margin = (support_plus - support_minus).abs();
            let pred_state = if support_plus >= support_minus { 1 } else { -1 };

            by_pos.entry(pos).or_default().push(Prediction {
                pos,
                block_id: block_id.clone(),
                pred_state,
                accepted: false,
                support_plus,
                support_minus,
                margin,
                confidence,
                donors: plus_donors.len() + minus_donors.len(),
                anchors: anchors.len(),
                best_vs_second_margin: 0.0,
                n_block_candidates: 0,
                method: "pbwt_hmm".to_string(),
                donor_entropy: (weighted_donor_entropy(plus_weights)
                    + weighted_donor_entropy(minus_weights))
                    / 2.0,
                effective_donors: weighted_effective_donors(plus_weights)
                    + weighted_effective_donors(minus_weights),
                pbwt_left_match_anchors: max_left_match(&plus_donors).max(max_left_match(&minus_donors)),
                pbwt_right_match_anchors: max_right_match(&plus_donors).max(max_right_match(&minus_donors)),
                retrieval: format!("{:?}", PbwtRetrieval::NeighborDebug),
            });
        }
    }

    by_pos
}

fn hmm_predictions_v2(
    block_candidates: &[BlockCandidate],
    panel: &PanelData,
    pbwt: &PbwtIndex,
    genetic_map: Option<&GeneticMap>,
    args: &Args,
) -> HashMap<u64, Vec<Prediction>> {
    let mut by_block: HashMap<String, Vec<BlockCandidate>> = HashMap::new();
    for bc in block_candidates {
        by_block
            .entry(bc.block_id.clone())
            .or_default()
            .push(bc.clone());
    }

    let mut by_pos: HashMap<u64, Vec<Prediction>> = HashMap::new();
    for (block_id, bcs) in by_block {
        let mut anchor_map: BTreeMap<u64, Anchor> = BTreeMap::new();
        let mut candidate_positions = BTreeSet::new();
        for bc in &bcs {
            candidate_positions.insert(bc.pos);
            for anchor in &bc.anchors {
                anchor_map.insert(anchor.pos, anchor.clone());
            }
        }
        let anchors: Vec<Anchor> = anchor_map.values().cloned().collect();
        if anchors.len() < args.min_anchors {
            continue;
        }

        let center_pos = *candidate_positions
            .iter()
            .nth(candidate_positions.len() / 2)
            .unwrap_or(&anchors[anchors.len() / 2].pos);

        let plus_donors = retrieve_path_donors(
            center_pos,
            &anchors,
            true,
            panel,
            pbwt,
            args,
            args.pbwt_retrieval,
        );
        let minus_donors = retrieve_path_donors(
            center_pos,
            &anchors,
            false,
            panel,
            pbwt,
            args,
            args.pbwt_retrieval,
        );

        if plus_donors.len() < args.min_donors || minus_donors.len() < args.min_donors {
            continue;
        }

        let anchor_state_by_pos: HashMap<u64, i8> =
            anchors.iter().map(|a| (a.pos, a.state)).collect();
        let mut sequence_positions: Vec<u64> = anchor_state_by_pos.keys().copied().collect();
        sequence_positions.extend(candidate_positions.iter().copied());
        sequence_positions.sort_unstable();
        sequence_positions.dedup();

        let plus_post = run_forward_backward_hmm(
            &plus_donors,
            &sequence_positions,
            &anchor_state_by_pos,
            true,
            panel,
            genetic_map,
            args,
        );
        let minus_post = run_forward_backward_hmm(
            &minus_donors,
            &sequence_positions,
            &anchor_state_by_pos,
            false,
            panel,
            genetic_map,
            args,
        );

        for pos in candidate_positions {
            let Some(site) = panel.sites.get(&pos) else {
                continue;
            };
            let Some(plus_weights) = plus_post.by_pos.get(&pos) else {
                continue;
            };
            let Some(minus_weights) = minus_post.by_pos.get(&pos) else {
                continue;
            };

            let mut plus_alt = 0.0;
            let mut plus_valid = 0.0;
            for (donor, w) in plus_donors.iter().zip(plus_weights) {
                match site.alleles[donor.hap] {
                    1 => {
                        plus_alt += *w;
                        plus_valid += *w;
                    }
                    0 => {
                        plus_valid += *w;
                    }
                    _ => {}
                }
            }

            let mut minus_ref = 0.0;
            let mut minus_valid = 0.0;
            for (donor, w) in minus_donors.iter().zip(minus_weights) {
                match site.alleles[donor.hap] {
                    0 => {
                        minus_ref += *w;
                        minus_valid += *w;
                    }
                    1 => {
                        minus_valid += *w;
                    }
                    _ => {}
                }
            }

            let mut terms = Vec::new();
            if plus_valid > 0.0 {
                terms.push(plus_alt / plus_valid);
            }
            if minus_valid > 0.0 {
                terms.push(minus_ref / minus_valid);
            }
            if terms.is_empty() {
                continue;
            }

            let prob_plus = terms.iter().sum::<f64>() / (terms.len() as f64);
            let support_plus = prob_plus;
            let support_minus = 1.0 - prob_plus;
            let confidence = support_plus.max(support_minus);
            let margin = (support_plus - support_minus).abs();
            let pred_state = if support_plus >= support_minus { 1 } else { -1 };

            by_pos.entry(pos).or_default().push(Prediction {
                pos,
                block_id: block_id.clone(),
                pred_state,
                accepted: false,
                support_plus,
                support_minus,
                margin,
                confidence,
                donors: plus_donors.len() + minus_donors.len(),
                anchors: anchors.len(),
                best_vs_second_margin: 0.0,
                n_block_candidates: 0,
                method: "pbwt_hmm_v2".to_string(),
                donor_entropy: (weighted_donor_entropy(plus_weights)
                    + weighted_donor_entropy(minus_weights))
                    / 2.0,
                effective_donors: weighted_effective_donors(plus_weights)
                    + weighted_effective_donors(minus_weights),
                pbwt_left_match_anchors: weighted_left_match(&plus_donors, plus_weights)
                    .max(weighted_left_match(&minus_donors, minus_weights)),
                pbwt_right_match_anchors: weighted_right_match(&plus_donors, plus_weights)
                    .max(weighted_right_match(&minus_donors, minus_weights)),
                retrieval: format!("{:?}", args.pbwt_retrieval),
            });
        }
    }

    by_pos
}

fn write_calls(
    path: &str,
    original_headers: &[String],
    rows: &[CallRow],
    chosen: &HashMap<u64, Prediction>,
    mode: Mode,
) -> Result<()> {
    if let Some(parent) = Path::new(path).parent() {
        fs::create_dir_all(parent)
            .with_context(|| format!("failed to create {}", parent.display()))?;
    }

    let prefix = match mode {
        Mode::Pbwt => "phase8_pbwt",
        Mode::PbwtHmm => "phase8_pbwt_hmm",
        Mode::PbwtV2 => "phase8_pbwt_v2",
        Mode::PbwtHmmV2 => "phase8_pbwt_hmm_v2",
    };
    let extra = [
        format!("{prefix}_filled"),
        format!("{prefix}_confidence"),
        format!("{prefix}_margin"),
        format!("{prefix}_support_plus"),
        format!("{prefix}_support_minus"),
        format!("{prefix}_donors"),
        format!("{prefix}_anchors"),
        format!("{prefix}_method"),
        format!("{prefix}_donor_entropy"),
        format!("{prefix}_effective_donors"),
        format!("{prefix}_pbwt_left_match_anchors"),
        format!("{prefix}_pbwt_right_match_anchors"),
        format!("{prefix}_retrieval"),
    ];

    let mut headers = original_headers.to_vec();
    for field in &extra {
        if !headers.contains(field) {
            headers.push(field.clone());
        }
    }

    let mut writer = WriterBuilder::new()
        .delimiter(b'\t')
        .terminator(Terminator::CRLF)
        .from_path(path)
        .with_context(|| format!("failed to write {path}"))?;
    writer.write_record(&headers)?;

    for row in rows {
        let mut out = row.clone();
        let pos = parse_u64(out.fields.get("pos")).unwrap_or(0);
        if let Some(pred) = chosen.get(&pos) {
            set_call_state(&mut out, pred.pred_state);
            out.fields
                .insert("block_id".to_string(), pred.block_id.clone());
            out.fields
                .insert(extra[0].clone(), "1".to_string());
            out.fields
                .insert(extra[1].clone(), format!("{:.6}", pred.confidence));
            out.fields
                .insert(extra[2].clone(), format!("{:.6}", pred.margin));
            out.fields
                .insert(extra[3].clone(), format!("{:.6}", pred.support_plus));
            out.fields
                .insert(extra[4].clone(), format!("{:.6}", pred.support_minus));
            out.fields
                .insert(extra[5].clone(), pred.donors.to_string());
            out.fields
                .insert(extra[6].clone(), pred.anchors.to_string());
            out.fields
                .insert(extra[7].clone(), pred.method.clone());
            out.fields
                .insert(extra[8].clone(), format!("{:.6}", pred.donor_entropy));
            out.fields
                .insert(extra[9].clone(), format!("{:.6}", pred.effective_donors));
            out.fields
                .insert(extra[10].clone(), pred.pbwt_left_match_anchors.to_string());
            out.fields
                .insert(extra[11].clone(), pred.pbwt_right_match_anchors.to_string());
            out.fields
                .insert(extra[12].clone(), pred.retrieval.clone());
        } else {
            out.fields.entry(extra[0].clone()).or_insert("0".to_string());
            for field in extra.iter().skip(1) {
                out.fields.entry(field.clone()).or_default();
            }
        }

        let record: Vec<String> = headers
            .iter()
            .map(|h| out.fields.get(h).cloned().unwrap_or_default())
            .collect();
        writer.write_record(record)?;
    }
    writer.flush()?;
    Ok(())
}

fn write_candidates(path: &str, preds: &[Prediction]) -> Result<()> {
    if let Some(parent) = Path::new(path).parent() {
        fs::create_dir_all(parent)
            .with_context(|| format!("failed to create {}", parent.display()))?;
    }
    let mut writer = WriterBuilder::new()
        .delimiter(b'\t')
        .terminator(Terminator::CRLF)
        .from_path(path)
        .with_context(|| format!("failed to write {path}"))?;
    writer.write_record([
        "pos",
        "block_id",
        "pred_state",
        "accepted",
        "support_plus",
        "support_minus",
        "margin",
        "confidence",
        "donors",
        "anchors",
        "best_vs_second_margin",
        "n_block_candidates",
        "method",
        "donor_entropy",
        "effective_donors",
        "pbwt_left_match_anchors",
        "pbwt_right_match_anchors",
        "retrieval",
    ])?;
    for p in preds {
        writer.write_record([
            p.pos.to_string(),
            p.block_id.clone(),
            p.pred_state.to_string(),
            if p.accepted { "1" } else { "0" }.to_string(),
            format!("{:.6}", p.support_plus),
            format!("{:.6}", p.support_minus),
            format!("{:.6}", p.margin),
            format!("{:.6}", p.confidence),
            p.donors.to_string(),
            p.anchors.to_string(),
            format!("{:.6}", p.best_vs_second_margin),
            p.n_block_candidates.to_string(),
            p.method.clone(),
            format!("{:.6}", p.donor_entropy),
            format!("{:.6}", p.effective_donors),
            p.pbwt_left_match_anchors.to_string(),
            p.pbwt_right_match_anchors.to_string(),
            p.retrieval.clone(),
        ])?;
    }
    writer.flush()?;
    Ok(())
}

fn main() -> Result<()> {
    let args = Args::parse();
    let variants = load_variants(&args.variant_json, &args.chrom, args.start, args.end)?;
    let variant_by_pos: HashMap<u64, Variant> =
        variants.iter().map(|v| (v.pos, v.clone())).collect();
    let variant_pos: HashSet<u64> = variants.iter().map(|v| v.pos).collect();

    let (headers, rows) = read_calls(&args.backbone_local_calls_tsv)?;
    let candidates = candidate_positions(&rows, &variant_pos);
    let anchors = load_anchors(&rows, &variant_pos);

    let mut block_candidates = Vec::new();
    for pos in &candidates {
        block_candidates.extend(block_candidates_for_position(
            *pos,
            &anchors,
            args.max_anchor_dist_bp,
            args.max_anchors_each_side,
            args.min_anchors,
        ));
    }

    let mut needed_positions: BTreeSet<u64> = BTreeSet::new();
    for bc in &block_candidates {
        needed_positions.insert(bc.pos);
        for anchor in &bc.anchors {
            needed_positions.insert(anchor.pos);
        }
    }

    let needed_variants: HashMap<u64, Variant> = needed_positions
        .iter()
        .filter_map(|pos| variant_by_pos.get(pos).map(|v| (*pos, v.clone())))
        .collect();

    let panel = load_panel_bcf(
        &args.panel_bcf,
        &args.chrom,
        args.start,
        args.end,
        &needed_variants,
    )?;
    let panel_positions: Vec<u64> = needed_positions
        .iter()
        .copied()
        .filter(|pos| panel.sites.contains_key(pos))
        .collect();
    let pbwt = build_pbwt_index(&panel, &panel_positions);

    let retrieval = match args.mode {
        Mode::Pbwt | Mode::PbwtHmm => PbwtRetrieval::NeighborDebug,
        Mode::PbwtV2 | Mode::PbwtHmmV2 => args.pbwt_retrieval,
    };
    let pbwt_method = match args.mode {
        Mode::Pbwt | Mode::PbwtHmm => "pbwt",
        Mode::PbwtV2 | Mode::PbwtHmmV2 => "pbwt_v2",
    };
    let pbwt_candidates_by_pos =
        pbwt_predictions(&block_candidates, &panel, &pbwt, &args, pbwt_method, retrieval);
    let (candidates_by_pos, chosen) = match args.mode {
        Mode::Pbwt | Mode::PbwtV2 => {
            let chosen = choose_predictions(&pbwt_candidates_by_pos, &args);
            (pbwt_candidates_by_pos, chosen)
        }
        Mode::PbwtHmm => {
            let genetic_map = load_genetic_map(args.genetic_map.as_deref())?;
            let hmm_candidates_by_pos = hmm_predictions(
                &block_candidates,
                &panel,
                &pbwt,
                genetic_map.as_ref(),
                &args,
            );
            let mut combined_candidates = pbwt_candidates_by_pos.clone();
            for (pos, preds) in hmm_candidates_by_pos {
                combined_candidates.entry(pos).or_default().extend(preds);
            }

            let mut chosen = choose_predictions(&pbwt_candidates_by_pos, &args);
            let hmm_chosen = choose_predictions(
                &combined_candidates
                    .iter()
                    .map(|(pos, preds)| {
                        (
                            *pos,
                            preds
                                .iter()
                                .filter(|pred| pred.method == "pbwt_hmm")
                                .cloned()
                                .collect::<Vec<_>>(),
                        )
                    })
                    .filter(|(_, preds)| !preds.is_empty())
                    .collect::<HashMap<_, _>>(),
                &args,
            );
            for (pos, pred) in hmm_chosen {
                chosen.insert(pos, pred);
            }

            (combined_candidates, chosen)
        }
        Mode::PbwtHmmV2 => {
            let genetic_map = load_genetic_map(args.genetic_map.as_deref())?;
            let hmm_candidates_by_pos =
                hmm_predictions_v2(&block_candidates, &panel, &pbwt, genetic_map.as_ref(), &args);
            let mut combined_candidates = pbwt_candidates_by_pos.clone();
            for (pos, preds) in hmm_candidates_by_pos {
                combined_candidates.entry(pos).or_default().extend(preds);
            }

            let mut chosen = choose_predictions(&pbwt_candidates_by_pos, &args);
            let hmm_chosen = choose_predictions(
                &combined_candidates
                    .iter()
                    .map(|(pos, preds)| {
                        (
                            *pos,
                            preds
                                .iter()
                                .filter(|pred| pred.method == "pbwt_hmm_v2")
                                .cloned()
                                .collect::<Vec<_>>(),
                        )
                    })
                    .filter(|(_, preds)| !preds.is_empty())
                    .collect::<HashMap<_, _>>(),
                &args,
            );
            for (pos, pred) in hmm_chosen {
                chosen.insert(pos, pred);
            }

            (combined_candidates, chosen)
        }
    };

    let mut candidate_rows = Vec::new();
    for (pos, preds) in &candidates_by_pos {
        let chosen_pos = chosen.get(pos);
        for p in preds {
            let mut pp = p.clone();
            if let Some(best) = chosen_pos {
                if best.block_id == pp.block_id
                    && best.pred_state == pp.pred_state
                    && (best.margin - pp.margin).abs() < 1e-9
                {
                    pp.accepted = true;
                    pp.best_vs_second_margin = best.best_vs_second_margin;
                    pp.n_block_candidates = best.n_block_candidates;
                }
            }
            candidate_rows.push(pp);
        }
    }
    candidate_rows.sort_by(|a, b| {
        a.pos
            .cmp(&b.pos)
            .then_with(|| a.block_id.cmp(&b.block_id))
            .then_with(|| b.margin.partial_cmp(&a.margin).unwrap_or(Ordering::Equal))
    });

    write_calls(&args.out_tsv, &headers, &rows, &chosen, args.mode)?;
    write_candidates(&args.out_candidates_tsv, &candidate_rows)?;

    let mode_name = match args.mode {
        Mode::Pbwt => "pbwt",
        Mode::PbwtHmm => "pbwt_hmm",
        Mode::PbwtV2 => "pbwt_v2",
        Mode::PbwtHmmV2 => "pbwt_hmm_v2",
    };
    let summary = Summary {
        mode: mode_name,
        panel_bcf: &args.panel_bcf,
        backbone_local_calls_tsv: &args.backbone_local_calls_tsv,
        variant_json: &args.variant_json,
        chrom: &args.chrom,
        start: args.start,
        end: args.end,
        backbone_anchors: anchors.len(),
        candidate_unphased_positions: candidates.len(),
        candidate_block_pairs: block_candidates.len(),
        panel_positions_needed: needed_positions.len(),
        panel_positions_found: panel.positions_found,
        panel_positions_missing: needed_positions
            .len()
            .saturating_sub(panel.positions_found),
        panel_allele_mismatch_records: panel.allele_mismatch,
        panel_nonbiallelic_records: panel.nonbiallelic,
        pbwt_index_positions: pbwt.by_pos.len(),
        pbwt_retrieval: match retrieval {
            PbwtRetrieval::ExactPrefix => "exact_prefix",
            PbwtRetrieval::NeighborDebug => "neighbor_debug",
        },
        hmm_decoder: match args.mode {
            Mode::PbwtHmmV2 => "forward_backward",
            Mode::PbwtHmm => "forward_only",
            Mode::Pbwt | Mode::PbwtV2 => "off",
        },
        accepted_sites: chosen.len(),
        out_tsv: &args.out_tsv,
    };
    let text = serde_json::to_string_pretty(&summary)?;
    fs::write(&args.out_summary_json, format!("{text}\n"))
        .with_context(|| format!("failed to write {}", args.out_summary_json))?;
    println!("{text}");
    Ok(())
}
