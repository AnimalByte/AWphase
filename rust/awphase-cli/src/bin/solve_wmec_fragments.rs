use std::cmp::Ordering;
use std::collections::{BinaryHeap, HashMap, HashSet};
use std::fs;
use std::path::Path;

use anyhow::{bail, Context, Result};
use clap::Parser;
use csv::{ReaderBuilder, Terminator, WriterBuilder};
use serde::Serialize;
use serde_json::Value;

#[derive(Parser, Debug)]
#[command(about = "Solve AWPhase Phase6C WMEC fragments")]
struct Args {
    #[arg(long)]
    fragments_tsv: String,
    #[arg(long)]
    variant_json: String,
    #[arg(long)]
    out_local_calls_tsv: String,
    #[arg(long)]
    out_components_tsv: String,
    #[arg(long)]
    out_summary_json: String,

    #[arg(long, default_value_t = 18)]
    max_exact_sites: usize,
    #[arg(long, default_value_t = 256)]
    max_component_sites: usize,
    #[arg(long, default_value_t = 2)]
    min_component_sites: usize,
    #[arg(long, default_value_t = 1)]
    min_fragments_per_component: usize,
    #[arg(long, default_value_t = 5)]
    local_refine_iters: usize,
}

#[derive(Debug, Clone)]
struct VariantRow {
    pos: u64,
}

#[derive(Debug, Clone)]
struct RawFragment {
    fragment_id: String,
    items: Vec<(u64, i8, f64)>,
}

#[derive(Debug, Clone)]
struct ProjectedFragment {
    sites: Vec<usize>,
    alleles: Vec<i8>,
    weights: Vec<f64>,
}

#[derive(Debug, Clone)]
struct SolvedSite {
    pos: u64,
    block_id: String,
    phase_state: i8,
    local_phase_state: i8,
    confidence: String,
    phase6_solver: String,
    phase6_component_id: String,
    phase6_component_sites: usize,
    phase6_component_fragments: usize,
    phase6_objective: String,
    phase6_margin: String,
    phase6_site_support: String,
    phase6_site_conflict: String,
}

#[derive(Debug, Clone)]
struct ComponentRow {
    component_id: String,
    start_pos: u64,
    end_pos: u64,
    span_bp: u64,
    n_sites: usize,
    n_fragments: usize,
    solver: String,
    objective: String,
    margin: String,
}

#[derive(Serialize)]
struct Summary<'a> {
    fragments_tsv: &'a str,
    variant_json: &'a str,
    raw_fragments: usize,
    raw_connected_components: usize,
    components_after_chunking: usize,
    components_solved: usize,
    components_skipped: usize,
    exact_components: usize,
    greedy_components: usize,
    total_projected_fragments: usize,
    variant_rows: usize,
    solved_sites: usize,
    out_local_calls_tsv: &'a str,
    max_exact_sites: usize,
    max_component_sites: usize,
}

#[derive(Debug, Clone)]
struct HeapEdge {
    weight: f64,
    src: usize,
    dst: usize,
    sign: i8,
}

impl PartialEq for HeapEdge {
    fn eq(&self, other: &Self) -> bool {
        self.weight == other.weight
            && self.src == other.src
            && self.dst == other.dst
            && self.sign == other.sign
    }
}

impl Eq for HeapEdge {}

impl PartialOrd for HeapEdge {
    fn partial_cmp(&self, other: &Self) -> Option<Ordering> {
        Some(self.cmp(other))
    }
}

impl Ord for HeapEdge {
    fn cmp(&self, other: &Self) -> Ordering {
        self.weight
            .partial_cmp(&other.weight)
            .unwrap_or(Ordering::Equal)
            .then_with(|| other.src.cmp(&self.src))
            .then_with(|| other.dst.cmp(&self.dst))
            .then_with(|| other.sign.cmp(&self.sign))
    }
}

#[derive(Debug, Default)]
struct Dsu {
    parent: HashMap<u64, u64>,
    rank: HashMap<u64, u8>,
}

impl Dsu {
    fn find(&mut self, x: u64) -> u64 {
        if !self.parent.contains_key(&x) {
            self.parent.insert(x, x);
            self.rank.insert(x, 0);
            return x;
        }

        let mut root = x;
        while self.parent[&root] != root {
            root = self.parent[&root];
        }

        let mut cur = x;
        while self.parent[&cur] != cur {
            let next = self.parent[&cur];
            self.parent.insert(cur, root);
            cur = next;
        }

        root
    }

    fn union(&mut self, a: u64, b: u64) {
        let mut ra = self.find(a);
        let mut rb = self.find(b);
        if ra == rb {
            return;
        }

        let rank_a = self.rank[&ra];
        let rank_b = self.rank[&rb];
        if rank_a < rank_b {
            std::mem::swap(&mut ra, &mut rb);
        }

        self.parent.insert(rb, ra);
        if rank_a == rank_b {
            self.rank.insert(ra, rank_a + 1);
        }
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

fn read_variant_json(path: &str) -> Result<Vec<VariantRow>> {
    let text = fs::read_to_string(path)
        .with_context(|| format!("failed to read variant JSON: {path}"))?;
    let value: Value = serde_json::from_str(&text)
        .with_context(|| format!("failed to parse variant JSON: {path}"))?;

    let mut rows = Vec::new();
    for row in as_array_root(&value)? {
        if let Some(pos) = row.get("pos").and_then(Value::as_u64) {
            rows.push(VariantRow { pos });
        } else if let Some(pos_str) = row.get("pos").and_then(Value::as_str) {
            if let Ok(pos) = pos_str.parse::<u64>() {
                rows.push(VariantRow { pos });
            }
        }
    }
    rows.sort_by_key(|r| r.pos);
    Ok(rows)
}

fn parse_json_vec_u64(text: &str) -> Result<Vec<u64>> {
    let values: Vec<Value> = serde_json::from_str(text)?;
    let mut out = Vec::with_capacity(values.len());
    for value in values {
        if let Some(x) = value.as_u64() {
            out.push(x);
        } else if let Some(x) = value.as_i64() {
            if x >= 0 {
                out.push(x as u64);
            }
        } else if let Some(x) = value.as_str() {
            out.push(x.parse()?);
        }
    }
    Ok(out)
}

fn parse_json_vec_i8(text: &str) -> Result<Vec<i8>> {
    let values: Vec<Value> = serde_json::from_str(text)?;
    let mut out = Vec::with_capacity(values.len());
    for value in values {
        if let Some(x) = value.as_i64() {
            out.push(x as i8);
        } else if let Some(x) = value.as_str() {
            out.push(x.parse()?);
        }
    }
    Ok(out)
}

fn parse_json_vec_f64(text: &str) -> Result<Vec<f64>> {
    let values: Vec<Value> = serde_json::from_str(text)?;
    let mut out = Vec::with_capacity(values.len());
    for value in values {
        if let Some(x) = value.as_f64() {
            out.push(x);
        } else if let Some(x) = value.as_str() {
            out.push(x.parse()?);
        }
    }
    Ok(out)
}

fn load_fragments(path: &str) -> Result<Vec<RawFragment>> {
    let mut reader = ReaderBuilder::new()
        .delimiter(b'\t')
        .from_path(path)
        .with_context(|| format!("failed to read fragments TSV: {path}"))?;

    let headers = reader.headers()?.clone();
    let idx_fragment = headers
        .iter()
        .position(|h| h == "fragment_id")
        .context("fragments TSV missing fragment_id column")?;
    let idx_positions = headers
        .iter()
        .position(|h| h == "positions_json")
        .context("fragments TSV missing positions_json column")?;
    let idx_alleles = headers
        .iter()
        .position(|h| h == "alleles_json")
        .context("fragments TSV missing alleles_json column")?;
    let idx_weights = headers
        .iter()
        .position(|h| h == "weights_json")
        .context("fragments TSV missing weights_json column")?;

    let mut out = Vec::new();
    for record in reader.records() {
        let record = record?;
        let positions = parse_json_vec_u64(record.get(idx_positions).unwrap_or(""))?;
        let alleles = parse_json_vec_i8(record.get(idx_alleles).unwrap_or(""))?;
        let weights = parse_json_vec_f64(record.get(idx_weights).unwrap_or(""))?;

        if positions.len() != alleles.len() || positions.len() != weights.len() {
            continue;
        }

        let mut items = Vec::new();
        for ((pos, allele), weight) in positions.into_iter().zip(alleles).zip(weights) {
            if !matches!(allele, -1 | 1) || weight <= 0.0 {
                continue;
            }
            items.push((pos, allele, weight));
        }

        if items.len() >= 2 {
            items.sort_by_key(|(pos, _, _)| *pos);
            out.push(RawFragment {
                fragment_id: record.get(idx_fragment).unwrap_or("").to_string(),
                items,
            });
        }
    }

    Ok(out)
}

fn fragment_cost_local(x: &[i8], sites: &[usize], alleles: &[i8], weights: &[f64]) -> f64 {
    let mut plus = 0.0;
    let mut minus = 0.0;
    for ((idx, allele), weight) in sites.iter().zip(alleles).zip(weights) {
        if x[*idx] != *allele {
            plus += *weight;
        }
        if x[*idx] != -*allele {
            minus += *weight;
        }
    }
    plus.min(minus)
}

fn fragment_orientation(x: &[i8], sites: &[usize], alleles: &[i8], weights: &[f64]) -> i8 {
    let mut plus = 0.0;
    let mut minus = 0.0;
    for ((idx, allele), weight) in sites.iter().zip(alleles).zip(weights) {
        if x[*idx] != *allele {
            plus += *weight;
        }
        if x[*idx] != -*allele {
            minus += *weight;
        }
    }
    if plus <= minus {
        1
    } else {
        -1
    }
}

fn total_cost(x: &[i8], fragments: &[ProjectedFragment]) -> f64 {
    fragments
        .iter()
        .map(|f| fragment_cost_local(x, &f.sites, &f.alleles, &f.weights))
        .sum()
}

fn exact_solve(n: usize, fragments: &[ProjectedFragment]) -> (Vec<i8>, f64, f64) {
    if n == 0 {
        return (Vec::new(), 0.0, 0.0);
    }

    let mut best = f64::INFINITY;
    let mut second = f64::INFINITY;
    let mut best_x = vec![1; n];
    let masks = 1_u64 << n.saturating_sub(1);

    for mask in 0..masks {
        let mut x = vec![1_i8; n];
        for (i, item) in x.iter_mut().enumerate().skip(1) {
            *item = if ((mask >> (i - 1)) & 1) == 1 { 1 } else { -1 };
        }

        let cost = total_cost(&x, fragments);
        if cost < best {
            second = best;
            best = cost;
            best_x = x;
        } else if cost < second {
            second = cost;
        }
    }

    let margin = if second.is_finite() { second - best } else { 0.0 };
    (best_x, best, margin)
}

fn greedy_init(n: usize, fragments: &[ProjectedFragment]) -> Vec<i8> {
    let mut pair: HashMap<(usize, usize), f64> = HashMap::new();

    for fragment in fragments {
        let k = fragment.sites.len();
        if k > 20 {
            continue;
        }
        for i in 0..k {
            for j in (i + 1)..k {
                let mut a = fragment.sites[i];
                let mut b = fragment.sites[j];
                let sign = fragment.alleles[i] * fragment.alleles[j];
                let weight = fragment.weights[i].min(fragment.weights[j]);
                if a > b {
                    std::mem::swap(&mut a, &mut b);
                }
                *pair.entry((a, b)).or_insert(0.0) += f64::from(sign) * weight;
            }
        }
    }

    let mut adj: Vec<Vec<(usize, i8, f64)>> = vec![Vec::new(); n];
    let mut pairs: Vec<_> = pair.into_iter().collect();
    pairs.sort_by_key(|((a, b), _)| (*a, *b));

    for ((a, b), value) in pairs {
        if value == 0.0 {
            continue;
        }
        let sign = if value > 0.0 { 1 } else { -1 };
        let weight = value.abs();
        adj[a].push((b, sign, weight));
        adj[b].push((a, sign, weight));
    }

    for edges in &mut adj {
        edges.sort_by(|a, b| {
            b.2.partial_cmp(&a.2)
                .unwrap_or(Ordering::Equal)
                .then_with(|| a.0.cmp(&b.0))
                .then_with(|| a.1.cmp(&b.1))
        });
    }

    let mut x = vec![0_i8; n];

    for start in 0..n {
        if x[start] != 0 {
            continue;
        }

        x[start] = 1;
        let mut heap = BinaryHeap::new();
        for (nb, sign, weight) in &adj[start] {
            heap.push(HeapEdge {
                weight: *weight,
                src: start,
                dst: *nb,
                sign: *sign,
            });
        }

        while let Some(edge) = heap.pop() {
            if x[edge.dst] != 0 {
                continue;
            }
            x[edge.dst] = x[edge.src] * edge.sign;
            for (nb, sign, weight) in &adj[edge.dst] {
                if x[*nb] == 0 {
                    heap.push(HeapEdge {
                        weight: *weight,
                        src: edge.dst,
                        dst: *nb,
                        sign: *sign,
                    });
                }
            }
        }
    }

    for value in &mut x {
        if *value == 0 {
            *value = 1;
        }
    }

    x
}

fn local_refine(
    mut x: Vec<i8>,
    fragments: &[ProjectedFragment],
    max_iter: usize,
) -> (Vec<i8>, f64) {
    let mut by_site: HashMap<usize, Vec<usize>> = HashMap::new();
    let mut costs = Vec::with_capacity(fragments.len());

    for (idx, fragment) in fragments.iter().enumerate() {
        let cost = fragment_cost_local(&x, &fragment.sites, &fragment.alleles, &fragment.weights);
        costs.push(cost);
        for site in &fragment.sites {
            by_site.entry(*site).or_default().push(idx);
        }
    }

    let mut best_cost: f64 = costs.iter().sum();

    for _ in 0..max_iter {
        let mut improved = false;

        for site in 0..x.len() {
            let Some(affected) = by_site.get(&site) else {
                continue;
            };

            let old_value = x[site];
            let old_sum: f64 = affected.iter().map(|idx| costs[*idx]).sum();

            x[site] = -x[site];
            let mut new_costs = Vec::new();
            let mut new_sum = 0.0;

            for idx in affected {
                let fragment = &fragments[*idx];
                let cost =
                    fragment_cost_local(&x, &fragment.sites, &fragment.alleles, &fragment.weights);
                new_sum += cost;
                new_costs.push((*idx, cost));
            }

            let delta = new_sum - old_sum;
            if delta < -1e-9 {
                for (idx, cost) in new_costs {
                    costs[idx] = cost;
                }
                best_cost += delta;
                improved = true;
            } else {
                x[site] = old_value;
            }
        }

        if !improved {
            break;
        }
    }

    (x, best_cost)
}

fn solve_greedy(
    n: usize,
    fragments: &[ProjectedFragment],
    max_iter: usize,
) -> (Vec<i8>, f64, String) {
    let x = greedy_init(n, fragments);
    let (x, cost) = local_refine(x, fragments, max_iter);
    (x, cost, String::new())
}

fn compute_site_support(
    x: &[i8],
    fragments: &[ProjectedFragment],
) -> (HashMap<usize, f64>, HashMap<usize, f64>) {
    let mut support = HashMap::new();
    let mut conflict = HashMap::new();

    for fragment in fragments {
        let orientation =
            fragment_orientation(x, &fragment.sites, &fragment.alleles, &fragment.weights);
        for ((idx, allele), weight) in fragment
            .sites
            .iter()
            .zip(&fragment.alleles)
            .zip(&fragment.weights)
        {
            if orientation * *allele == x[*idx] {
                *support.entry(*idx).or_insert(0.0) += *weight;
            } else {
                *conflict.entry(*idx).or_insert(0.0) += *weight;
            }
        }
    }

    (support, conflict)
}

fn project_fragments_to_positions(
    fragments: &[RawFragment],
    positions: &[u64],
) -> Vec<ProjectedFragment> {
    let pos_to_idx: HashMap<u64, usize> = positions
        .iter()
        .enumerate()
        .map(|(idx, pos)| (*pos, idx))
        .collect();

    let mut out = Vec::new();
    for fragment in fragments {
        let mut sites = Vec::new();
        let mut alleles = Vec::new();
        let mut weights = Vec::new();

        for (pos, allele, weight) in &fragment.items {
            if let Some(idx) = pos_to_idx.get(pos) {
                sites.push(*idx);
                alleles.push(*allele);
                weights.push(*weight);
            }
        }

        if sites.len() >= 2 {
            let _ = &fragment.fragment_id;
            out.push(ProjectedFragment {
                sites,
                alleles,
                weights,
            });
        }
    }

    out
}

fn format_six(value: f64) -> String {
    format!("{value:.6}")
}

fn write_local_calls(path: &str, rows: &[SolvedSite]) -> Result<()> {
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
        "phase_state",
        "local_phase_state",
        "confidence",
        "phase6_solver",
        "phase6_component_id",
        "phase6_component_sites",
        "phase6_component_fragments",
        "phase6_objective",
        "phase6_margin",
        "phase6_site_support",
        "phase6_site_conflict",
    ])?;

    for row in rows {
        writer.write_record([
            row.pos.to_string(),
            row.block_id.clone(),
            row.phase_state.to_string(),
            row.local_phase_state.to_string(),
            row.confidence.clone(),
            row.phase6_solver.clone(),
            row.phase6_component_id.clone(),
            row.phase6_component_sites.to_string(),
            row.phase6_component_fragments.to_string(),
            row.phase6_objective.clone(),
            row.phase6_margin.clone(),
            row.phase6_site_support.clone(),
            row.phase6_site_conflict.clone(),
        ])?;
    }
    writer.flush()?;
    Ok(())
}

fn write_components(path: &str, rows: &[ComponentRow]) -> Result<()> {
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
        "component_id",
        "start_pos",
        "end_pos",
        "span_bp",
        "n_sites",
        "n_fragments",
        "solver",
        "objective",
        "margin",
    ])?;

    for row in rows {
        writer.write_record([
            row.component_id.clone(),
            row.start_pos.to_string(),
            row.end_pos.to_string(),
            row.span_bp.to_string(),
            row.n_sites.to_string(),
            row.n_fragments.to_string(),
            row.solver.clone(),
            row.objective.clone(),
            row.margin.clone(),
        ])?;
    }
    writer.flush()?;
    Ok(())
}

fn main() -> Result<()> {
    let args = Args::parse();

    let variant_rows = read_variant_json(&args.variant_json)?;
    let variant_positions: Vec<u64> = variant_rows.iter().map(|r| r.pos).collect();
    let variant_pos_set: HashSet<u64> = variant_positions.iter().copied().collect();

    let mut raw_fragments = load_fragments(&args.fragments_tsv)?;
    for fragment in &mut raw_fragments {
        fragment
            .items
            .retain(|(pos, _, _)| variant_pos_set.contains(pos));
    }
    raw_fragments.retain(|f| f.items.len() >= 2);

    let mut dsu = Dsu::default();
    for fragment in &raw_fragments {
        for (pos, _, _) in &fragment.items {
            dsu.find(*pos);
        }
        if let Some((first, _, _)) = fragment.items.first() {
            for (pos, _, _) in fragment.items.iter().skip(1) {
                dsu.union(*first, *pos);
            }
        }
    }

    let mut comp_pos: HashMap<u64, HashSet<u64>> = HashMap::new();
    for fragment in &raw_fragments {
        for (pos, _, _) in &fragment.items {
            let root = dsu.find(*pos);
            comp_pos.entry(root).or_default().insert(*pos);
        }
    }

    let mut raw_components: Vec<Vec<u64>> = comp_pos
        .into_values()
        .map(|set| {
            let mut positions: Vec<u64> = set.into_iter().collect();
            positions.sort_unstable();
            positions
        })
        .filter(|positions| positions.len() >= args.min_component_sites)
        .collect();
    raw_components.sort_by_key(|positions| (positions[0], *positions.last().unwrap()));

    let mut components = Vec::new();
    for component in &raw_components {
        if component.len() <= args.max_component_sites {
            components.push(component.clone());
        } else {
            for chunk in component.chunks(args.max_component_sites) {
                if chunk.len() >= args.min_component_sites {
                    components.push(chunk.to_vec());
                }
            }
        }
    }

    let mut solved: HashMap<u64, SolvedSite> = HashMap::new();
    let mut component_rows = Vec::new();

    let mut exact_components = 0;
    let mut greedy_components = 0;
    let mut skipped_components = 0;
    let mut total_projected_fragments = 0;

    for (ci, positions) in components.iter().enumerate() {
        let projected = project_fragments_to_positions(&raw_fragments, positions);

        if projected.len() < args.min_fragments_per_component {
            skipped_components += 1;
            continue;
        }

        let n = positions.len();
        total_projected_fragments += projected.len();

        let (x, objective, margin, solver) = if n <= args.max_exact_sites {
            let (x, objective, margin) = exact_solve(n, &projected);
            exact_components += 1;
            (x, objective, margin.to_string(), "exact".to_string())
        } else {
            let (x, objective, margin) = solve_greedy(n, &projected, args.local_refine_iters);
            greedy_components += 1;
            (x, objective, margin, "greedy_refine".to_string())
        };

        let (support, conflict) = compute_site_support(&x, &projected);
        let block_id = format!("phase6_comp_{}", ci + 1);

        for (idx, pos) in positions.iter().enumerate() {
            let sup = support.get(&idx).copied().unwrap_or(0.0);
            let con = conflict.get(&idx).copied().unwrap_or(0.0);
            let confidence = if sup + con > 0.0 {
                sup / (sup + con)
            } else {
                0.0
            };

            solved.insert(
                *pos,
                SolvedSite {
                    pos: *pos,
                    block_id: block_id.clone(),
                    phase_state: x[idx],
                    local_phase_state: x[idx],
                    confidence: format_six(confidence),
                    phase6_solver: solver.clone(),
                    phase6_component_id: block_id.clone(),
                    phase6_component_sites: n,
                    phase6_component_fragments: projected.len(),
                    phase6_objective: format_six(objective),
                    phase6_margin: margin.clone(),
                    phase6_site_support: format_six(sup),
                    phase6_site_conflict: format_six(con),
                },
            );
        }

        component_rows.push(ComponentRow {
            component_id: block_id,
            start_pos: positions[0],
            end_pos: *positions.last().unwrap(),
            span_bp: positions.last().unwrap() - positions[0] + 1,
            n_sites: n,
            n_fragments: projected.len(),
            solver,
            objective: format_six(objective),
            margin,
        });
    }

    let mut out_rows = Vec::with_capacity(variant_rows.len());
    for row in &variant_rows {
        if let Some(site) = solved.get(&row.pos) {
            out_rows.push(site.clone());
        } else {
            out_rows.push(SolvedSite {
                pos: row.pos,
                block_id: "unassigned".to_string(),
                phase_state: 0,
                local_phase_state: 0,
                confidence: "0.000000".to_string(),
                phase6_solver: "unassigned".to_string(),
                phase6_component_id: String::new(),
                phase6_component_sites: 0,
                phase6_component_fragments: 0,
                phase6_objective: String::new(),
                phase6_margin: String::new(),
                phase6_site_support: "0.000000".to_string(),
                phase6_site_conflict: "0.000000".to_string(),
            });
        }
    }

    write_local_calls(&args.out_local_calls_tsv, &out_rows)?;
    write_components(&args.out_components_tsv, &component_rows)?;

    let summary = Summary {
        fragments_tsv: &args.fragments_tsv,
        variant_json: &args.variant_json,
        raw_fragments: raw_fragments.len(),
        raw_connected_components: raw_components.len(),
        components_after_chunking: components.len(),
        components_solved: component_rows.len(),
        components_skipped: skipped_components,
        exact_components,
        greedy_components,
        total_projected_fragments,
        variant_rows: variant_rows.len(),
        solved_sites: solved.len(),
        out_local_calls_tsv: &args.out_local_calls_tsv,
        max_exact_sites: args.max_exact_sites,
        max_component_sites: args.max_component_sites,
    };

    let summary_text = serde_json::to_string_pretty(&summary)?;
    fs::write(&args.out_summary_json, format!("{summary_text}\n"))
        .with_context(|| format!("failed to write {}", args.out_summary_json))?;
    println!("{summary_text}");

    Ok(())
}
