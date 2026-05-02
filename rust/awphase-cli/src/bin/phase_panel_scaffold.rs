use std::cmp::Ordering;
use std::collections::{BTreeMap, HashMap, HashSet};
use std::fs;
use std::path::Path;
use std::process::Command;

use anyhow::{bail, Context, Result};
use clap::Parser;
use csv::{ReaderBuilder, Terminator, WriterBuilder};
use serde::Serialize;
use serde_json::Value;

#[derive(Parser, Debug)]
#[command(about = "Panel-assisted block scaffolding for AWPhase short-read phase blocks")]
struct Args {
    #[arg(long)]
    panel_bcf: String,
    #[arg(long)]
    input_tsv: String,
    #[arg(long)]
    variant_json: String,
    #[arg(long)]
    chrom: String,
    #[arg(long)]
    start: u64,
    #[arg(long)]
    end: u64,
    #[arg(long)]
    genetic_map: String,
    #[arg(long)]
    out_tsv: String,
    #[arg(long)]
    out_joins_tsv: String,
    #[arg(long)]
    out_summary_json: String,

    #[arg(long, default_value_t = 2)]
    min_anchors_per_block: usize,
    #[arg(long, default_value_t = 48)]
    max_anchors_per_block: usize,
    #[arg(long, default_value_t = 64)]
    top_k: usize,
    #[arg(long, default_value_t = 32)]
    min_donors: usize,
    #[arg(long, default_value_t = 0.55)]
    min_confidence: f64,
    #[arg(long, default_value_t = 8.0)]
    min_margin: f64,
    #[arg(long, default_value_t = 0.50)]
    max_gap_cm: f64,
    #[arg(long, default_value_t = 2_000_000)]
    max_gap_bp: u64,
    #[arg(long, default_value_t = 0.05)]
    anchor_decay_cm: f64,
    #[arg(long, default_value_t = 0.25)]
    mismatch_penalty: f64,
    #[arg(long, default_value_t = 0.20)]
    match_bonus: f64,
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
}

#[derive(Debug, Clone)]
struct BlockInfo {
    block_id: String,
    start: u64,
    end: u64,
    anchors: Vec<Anchor>,
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

#[derive(Debug)]
struct GeneticMap {
    bp: Vec<u64>,
    cm: Vec<f64>,
}

#[derive(Debug, Clone)]
struct PathSupport {
    support: f64,
    donors: usize,
}

#[derive(Debug, Clone)]
struct BoundaryDecision {
    left_block: String,
    right_block: String,
    left_start: u64,
    left_end: u64,
    right_start: u64,
    right_end: u64,
    gap_bp: u64,
    gap_cm: f64,
    left_anchors: usize,
    right_anchors: usize,
    same_support: f64,
    flip_support: f64,
    same_donors: usize,
    flip_donors: usize,
    decision: String,
    accepted: bool,
    confidence: f64,
    margin: f64,
    reason: String,
}

#[derive(Debug, Clone)]
struct ScaffoldAssignment {
    scaffold_id: String,
    orientation: i8,
    joined: bool,
}

#[derive(Serialize)]
struct Summary<'a> {
    panel_bcf: &'a str,
    input_tsv: &'a str,
    variant_json: &'a str,
    chrom: &'a str,
    start: u64,
    end: u64,
    input_blocks: usize,
    scaffold_blocks: usize,
    joins_considered: usize,
    joins_accepted: usize,
    joins_same: usize,
    joins_flip: usize,
    joins_abstained: usize,
    panel_positions_needed: usize,
    panel_positions_found: usize,
    panel_positions_missing: usize,
    panel_allele_mismatch_records: usize,
    panel_nonbiallelic_records: usize,
    min_anchors_per_block: usize,
    max_anchors_per_block: usize,
    top_k: usize,
    min_donors: usize,
    min_confidence: f64,
    min_margin: f64,
    max_gap_cm: f64,
    max_gap_bp: u64,
    out_tsv: &'a str,
    out_joins_tsv: &'a str,
}

fn parse_i8(value: Option<&String>) -> i8 {
    value
        .and_then(|v| v.parse::<f64>().ok())
        .map(|v| v as i8)
        .unwrap_or(0)
}

fn parse_u64(value: Option<&String>) -> Option<u64> {
    value.and_then(|v| v.parse::<f64>().ok()).map(|v| v as u64)
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
    let text =
        fs::read_to_string(path).with_context(|| format!("failed to read variant JSON: {path}"))?;
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

fn parse_panel_allele(value: &str) -> i8 {
    match value {
        "0" => 0,
        "1" => 1,
        _ => -1,
    }
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

fn load_genetic_map(path: &str) -> Result<GeneticMap> {
    let text =
        fs::read_to_string(path).with_context(|| format!("failed to read genetic map: {path}"))?;
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
    Ok(GeneticMap {
        bp: pairs.iter().map(|x| x.0).collect(),
        cm: pairs.iter().map(|x| x.1).collect(),
    })
}

fn interp_cm(pos: u64, map: &GeneticMap) -> f64 {
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

fn load_blocks(rows: &[CallRow], variant_pos: &HashSet<u64>) -> Vec<BlockInfo> {
    let mut by_block: HashMap<String, Vec<Anchor>> = HashMap::new();
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
        by_block
            .entry(block_id)
            .or_default()
            .push(Anchor { pos, state });
    }

    let mut blocks = Vec::new();
    for (block_id, mut anchors) in by_block {
        anchors.sort_by_key(|a| a.pos);
        let Some(first) = anchors.first() else {
            continue;
        };
        let Some(last) = anchors.last() else {
            continue;
        };
        blocks.push(BlockInfo {
            block_id,
            start: first.pos,
            end: last.pos,
            anchors,
        });
    }
    blocks.sort_by(|a, b| {
        a.start
            .cmp(&b.start)
            .then_with(|| a.end.cmp(&b.end))
            .then_with(|| a.block_id.cmp(&b.block_id))
    });
    blocks
}

fn boundary_pos(left: &BlockInfo, right: &BlockInfo) -> u64 {
    if left.end <= right.start {
        left.end + (right.start - left.end) / 2
    } else {
        right.start + (left.end - right.start) / 2
    }
}

fn gap_bp(left: &BlockInfo, right: &BlockInfo) -> u64 {
    right.start.saturating_sub(left.end)
}

fn gap_cm(left: &BlockInfo, right: &BlockInfo, map: &GeneticMap) -> f64 {
    if right.start <= left.end {
        0.0
    } else {
        (interp_cm(right.start, map) - interp_cm(left.end, map)).abs()
    }
}

fn select_anchors(
    block: &BlockInfo,
    boundary: u64,
    panel: &PanelData,
    map: &GeneticMap,
    max_anchors: usize,
) -> Vec<Anchor> {
    let boundary_cm = interp_cm(boundary, map);
    let mut anchors: Vec<_> = block
        .anchors
        .iter()
        .filter(|a| panel.sites.contains_key(&a.pos))
        .cloned()
        .collect();
    anchors.sort_by(|a, b| {
        let da = (interp_cm(a.pos, map) - boundary_cm).abs();
        let db = (interp_cm(b.pos, map) - boundary_cm).abs();
        da.partial_cmp(&db)
            .unwrap_or(Ordering::Equal)
            .then_with(|| a.pos.cmp(&b.pos))
    });
    anchors.truncate(max_anchors);
    anchors.sort_by_key(|a| a.pos);
    anchors
}

fn anchor_weight(anchor_pos: u64, boundary: u64, map: &GeneticMap, args: &Args) -> f64 {
    let d_cm = (interp_cm(anchor_pos, map) - interp_cm(boundary, map)).abs();
    1.0 / (1.0 + d_cm / args.anchor_decay_cm.max(1e-6))
}

fn score_anchor_set(
    hap: usize,
    anchors: &[Anchor],
    plus_path: bool,
    boundary: u64,
    panel: &PanelData,
    map: &GeneticMap,
    args: &Args,
) -> (f64, usize) {
    let mut score = 0.0;
    let mut matches = 0;
    for anchor in anchors {
        let Some(site) = panel.sites.get(&anchor.pos) else {
            continue;
        };
        let expected = anchor_expected_allele(anchor.state, plus_path);
        let weight = anchor_weight(anchor.pos, boundary, map, args);
        match site.alleles[hap] {
            allele if allele == expected => {
                score += weight;
                matches += 1;
            }
            0 | 1 => {
                score -= args.mismatch_penalty * weight;
            }
            _ => {}
        }
    }
    (score, matches)
}

fn support_for_path(
    left: &[Anchor],
    right: &[Anchor],
    left_plus_path: bool,
    same_orientation: bool,
    boundary: u64,
    panel: &PanelData,
    map: &GeneticMap,
    args: &Args,
) -> PathSupport {
    let right_plus_path = if same_orientation {
        left_plus_path
    } else {
        !left_plus_path
    };
    let mut donor_scores = Vec::new();
    for hap in 0..panel.n_haps {
        let (left_score, left_matches) =
            score_anchor_set(hap, left, left_plus_path, boundary, panel, map, args);
        let (right_score, right_matches) =
            score_anchor_set(hap, right, right_plus_path, boundary, panel, map, args);
        if left_score <= 0.0 || right_score <= 0.0 {
            continue;
        }
        let shared_bonus = args.match_bonus * left_matches.min(right_matches) as f64;
        let score = left_score + right_score + shared_bonus;
        if score > 0.0 {
            donor_scores.push(score);
        }
    }
    donor_scores.sort_by(|a, b| b.partial_cmp(a).unwrap_or(Ordering::Equal));
    donor_scores.truncate(args.top_k);
    PathSupport {
        support: donor_scores.iter().sum(),
        donors: donor_scores.len(),
    }
}

fn relation_support(
    left: &[Anchor],
    right: &[Anchor],
    same_orientation: bool,
    boundary: u64,
    panel: &PanelData,
    map: &GeneticMap,
    args: &Args,
) -> PathSupport {
    let plus = support_for_path(
        left,
        right,
        true,
        same_orientation,
        boundary,
        panel,
        map,
        args,
    );
    let minus = support_for_path(
        left,
        right,
        false,
        same_orientation,
        boundary,
        panel,
        map,
        args,
    );
    PathSupport {
        support: plus.support + minus.support,
        donors: plus.donors + minus.donors,
    }
}

fn score_boundary_pair(
    left: &BlockInfo,
    right: &BlockInfo,
    panel: &PanelData,
    map: &GeneticMap,
    args: &Args,
) -> BoundaryDecision {
    let boundary = boundary_pos(left, right);
    let left_anchors = select_anchors(left, boundary, panel, map, args.max_anchors_per_block);
    let right_anchors = select_anchors(right, boundary, panel, map, args.max_anchors_per_block);
    let gap_bp = gap_bp(left, right);
    let gap_cm = gap_cm(left, right, map);

    let reason;
    let mut decision = "abstain".to_string();
    let mut accepted = false;
    let mut confidence = 0.0;
    let mut margin = 0.0;
    let mut same = PathSupport {
        support: 0.0,
        donors: 0,
    };
    let mut flip = same.clone();

    if left_anchors.len() < args.min_anchors_per_block
        || right_anchors.len() < args.min_anchors_per_block
    {
        reason = "insufficient_panel_anchors".to_string();
    } else if gap_bp > args.max_gap_bp || gap_cm > args.max_gap_cm {
        reason = "gap_exceeds_limit".to_string();
    } else {
        same = relation_support(
            &left_anchors,
            &right_anchors,
            true,
            boundary,
            panel,
            map,
            args,
        );
        flip = relation_support(
            &left_anchors,
            &right_anchors,
            false,
            boundary,
            panel,
            map,
            args,
        );
        let total = same.support + flip.support;
        if total <= 0.0 {
            reason = "no_panel_support".to_string();
        } else {
            let best_donors;
            if same.support >= flip.support {
                decision = "same".to_string();
                confidence = same.support / total;
                best_donors = same.donors;
            } else {
                decision = "flip".to_string();
                confidence = flip.support / total;
                best_donors = flip.donors;
            }
            margin = (same.support - flip.support).abs();
            accepted = best_donors >= args.min_donors
                && confidence >= args.min_confidence
                && margin >= args.min_margin;
            reason = if accepted {
                "accepted".to_string()
            } else if best_donors < args.min_donors {
                "insufficient_donors".to_string()
            } else if confidence < args.min_confidence {
                "low_confidence".to_string()
            } else {
                "low_margin".to_string()
            };
        }
    }

    BoundaryDecision {
        left_block: left.block_id.clone(),
        right_block: right.block_id.clone(),
        left_start: left.start,
        left_end: left.end,
        right_start: right.start,
        right_end: right.end,
        gap_bp,
        gap_cm,
        left_anchors: left_anchors.len(),
        right_anchors: right_anchors.len(),
        same_support: same.support,
        flip_support: flip.support,
        same_donors: same.donors,
        flip_donors: flip.donors,
        decision,
        accepted,
        confidence,
        margin,
        reason,
    }
}

fn build_assignments(
    blocks: &[BlockInfo],
    decisions: &[BoundaryDecision],
) -> HashMap<String, ScaffoldAssignment> {
    let mut out = HashMap::new();
    if blocks.is_empty() {
        return out;
    }

    let mut scaffold_idx = 1usize;
    let mut current_scaffold = format!("phase8_scaffold_{scaffold_idx}");
    let mut prev_orientation = 1i8;
    out.insert(
        blocks[0].block_id.clone(),
        ScaffoldAssignment {
            scaffold_id: current_scaffold.clone(),
            orientation: prev_orientation,
            joined: false,
        },
    );

    for (idx, block) in blocks.iter().enumerate().skip(1) {
        let decision = decisions.get(idx - 1);
        let joined = decision.map(|d| d.accepted).unwrap_or(false);
        let orientation = if joined {
            match decision.map(|d| d.decision.as_str()) {
                Some("same") => prev_orientation,
                Some("flip") => -prev_orientation,
                _ => 1,
            }
        } else {
            scaffold_idx += 1;
            current_scaffold = format!("phase8_scaffold_{scaffold_idx}");
            1
        };
        out.insert(
            block.block_id.clone(),
            ScaffoldAssignment {
                scaffold_id: current_scaffold.clone(),
                orientation,
                joined,
            },
        );
        prev_orientation = orientation;
    }

    out
}

fn write_scaffolded_calls(
    path: &str,
    original_headers: &[String],
    rows: &[CallRow],
    assignments: &HashMap<String, ScaffoldAssignment>,
) -> Result<()> {
    if let Some(parent) = Path::new(path).parent() {
        fs::create_dir_all(parent)
            .with_context(|| format!("failed to create {}", parent.display()))?;
    }
    let extra = [
        "phase8_block_scaffold_source_block".to_string(),
        "phase8_block_scaffold_id".to_string(),
        "phase8_block_scaffold_orientation".to_string(),
        "phase8_block_scaffold_joined".to_string(),
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
        let state = call_state(row);
        let source_block = clean_block_id(row.fields.get("block_id"));
        if state != 0 {
            if let Some(block_id) = source_block {
                if let Some(assign) = assignments.get(&block_id) {
                    set_call_state(&mut out, state * assign.orientation);
                    out.fields
                        .insert("block_id".to_string(), assign.scaffold_id.clone());
                    out.fields.insert(extra[0].clone(), block_id.clone());
                    out.fields
                        .insert(extra[1].clone(), assign.scaffold_id.clone());
                    out.fields
                        .insert(extra[2].clone(), assign.orientation.to_string());
                    out.fields.insert(
                        extra[3].clone(),
                        if assign.joined { "1" } else { "0" }.to_string(),
                    );
                }
            }
        }
        for field in &extra {
            out.fields.entry(field.clone()).or_default();
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

fn write_joins(path: &str, decisions: &[BoundaryDecision]) -> Result<()> {
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
        "left_block",
        "right_block",
        "left_start",
        "left_end",
        "right_start",
        "right_end",
        "gap_bp",
        "gap_cm",
        "left_anchors",
        "right_anchors",
        "same_support",
        "flip_support",
        "same_donors",
        "flip_donors",
        "decision",
        "accepted",
        "confidence",
        "margin",
        "reason",
    ])?;
    for d in decisions {
        writer.write_record([
            d.left_block.clone(),
            d.right_block.clone(),
            d.left_start.to_string(),
            d.left_end.to_string(),
            d.right_start.to_string(),
            d.right_end.to_string(),
            d.gap_bp.to_string(),
            format!("{:.8}", d.gap_cm),
            d.left_anchors.to_string(),
            d.right_anchors.to_string(),
            format!("{:.6}", d.same_support),
            format!("{:.6}", d.flip_support),
            d.same_donors.to_string(),
            d.flip_donors.to_string(),
            d.decision.clone(),
            if d.accepted { "1" } else { "0" }.to_string(),
            format!("{:.8}", d.confidence),
            format!("{:.6}", d.margin),
            d.reason.clone(),
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
    let (headers, rows) = read_calls(&args.input_tsv)?;
    let blocks = load_blocks(&rows, &variant_pos);
    let genetic_map = load_genetic_map(&args.genetic_map)?;

    let mut needed_positions = HashSet::new();
    for block in &blocks {
        for anchor in &block.anchors {
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

    let mut decisions = Vec::new();
    for pair in blocks.windows(2) {
        decisions.push(score_boundary_pair(
            &pair[0],
            &pair[1],
            &panel,
            &genetic_map,
            &args,
        ));
    }
    let assignments = build_assignments(&blocks, &decisions);
    write_scaffolded_calls(&args.out_tsv, &headers, &rows, &assignments)?;
    write_joins(&args.out_joins_tsv, &decisions)?;

    let joins_accepted = decisions.iter().filter(|d| d.accepted).count();
    let joins_same = decisions
        .iter()
        .filter(|d| d.accepted && d.decision == "same")
        .count();
    let joins_flip = decisions
        .iter()
        .filter(|d| d.accepted && d.decision == "flip")
        .count();
    let scaffold_blocks: HashSet<_> = assignments
        .values()
        .map(|a| a.scaffold_id.clone())
        .collect();
    let summary = Summary {
        panel_bcf: &args.panel_bcf,
        input_tsv: &args.input_tsv,
        variant_json: &args.variant_json,
        chrom: &args.chrom,
        start: args.start,
        end: args.end,
        input_blocks: blocks.len(),
        scaffold_blocks: scaffold_blocks.len(),
        joins_considered: decisions.len(),
        joins_accepted,
        joins_same,
        joins_flip,
        joins_abstained: decisions.len().saturating_sub(joins_accepted),
        panel_positions_needed: needed_positions.len(),
        panel_positions_found: panel.positions_found,
        panel_positions_missing: needed_positions.len().saturating_sub(panel.positions_found),
        panel_allele_mismatch_records: panel.allele_mismatch,
        panel_nonbiallelic_records: panel.nonbiallelic,
        min_anchors_per_block: args.min_anchors_per_block,
        max_anchors_per_block: args.max_anchors_per_block,
        top_k: args.top_k,
        min_donors: args.min_donors,
        min_confidence: args.min_confidence,
        min_margin: args.min_margin,
        max_gap_cm: args.max_gap_cm,
        max_gap_bp: args.max_gap_bp,
        out_tsv: &args.out_tsv,
        out_joins_tsv: &args.out_joins_tsv,
    };
    let text = serde_json::to_string_pretty(&summary)?;
    fs::write(&args.out_summary_json, format!("{text}\n"))
        .with_context(|| format!("failed to write {}", args.out_summary_json))?;
    println!("{text}");
    Ok(())
}
