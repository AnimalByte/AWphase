#!/usr/bin/env python3
import argparse, csv, json
from collections import defaultdict
from pathlib import Path


def read_tsv(path):
    with open(path) as fh:
        return list(csv.DictReader(fh, delimiter='\t'))


def write_tsv(path, rows, fields):
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    with open(path, 'w', newline='') as fh:
        w = csv.DictWriter(fh, fieldnames=fields, delimiter='\t')
        w.writeheader()
        w.writerows(rows)


def load_site_priors(path):
    if not path or not Path(path).exists():
        return {}
    obj = json.load(open(path))
    entries = obj.get('entries', obj if isinstance(obj, list) else [])
    out = {}
    for e in entries:
        try:
            out[int(e['pos'])] = float(e.get('bias', 0.0))
        except Exception:
            pass
    return out


def load_local(path, phase_col='local_phase_state'):
    rows = []
    for r in read_tsv(path):
        try:
            rows.append({
                'pos': int(r['pos']),
                'phase_state': int(r.get(phase_col, 0) or 0),
                'confidence': float(r.get('confidence', 0.0) or 0.0),
                'block_id': r.get('block_id', ''),
            })
        except Exception:
            pass
    return sorted(rows, key=lambda x: x['pos'])


def load_selected_ids(path):
    ids = set()
    if not path or not Path(path).exists():
        return ids
    for r in read_tsv(path):
        rid = str(r.get('read_id', '')).strip()
        if rid:
            ids.add(rid)
    return ids


def load_haplotag(path):
    tags = {}
    if not path or not Path(path).exists():
        return tags
    for r in read_tsv(path):
        rid = str(r.get('read_id', '')).strip()
        if not rid:
            continue
        hp = str(r.get('hp', r.get('HP', r.get('haplotype', 'UNK')))).upper()
        if hp in {'1', 'HP1'}:
            hp = 'HP1'
        elif hp in {'2', 'HP2'}:
            hp = 'HP2'
        else:
            hp = 'UNK'
        tags[rid] = hp
    return tags


def load_super(path):
    support = defaultdict(float)
    if not path or not Path(path).exists():
        return support
    for r in read_tsv(path):
        try:
            pos = int(r.get('site_pos', r.get('pos')))
        except Exception:
            continue
        n = 1.0
        for k in ('support_reads', 'n_reads', 'support', 'read_count', 'count'):
            if str(r.get(k, '')).strip():
                try:
                    n = max(1.0, float(r[k]))
                except Exception:
                    pass
                break
        support[pos] += n
    return support


def load_obs_support(obs_tsv, selected_ids, hap_tags):
    selected_support = defaultdict(int)
    tagged_support = defaultdict(int)
    hp1_reads = defaultdict(int)
    hp2_reads = defaultdict(int)

    if not obs_tsv or not Path(obs_tsv).exists():
        return selected_support, tagged_support, hp1_reads, hp2_reads

    with open(obs_tsv) as fh:
        reader = csv.DictReader(fh, delimiter='\t')
        for r in reader:
            rid = str(r.get('read_id', '')).strip()
            if not rid:
                continue
            try:
                pos = int(r.get('site_pos', r.get('pos')))
            except Exception:
                continue
            try:
                allele = int(float(r.get('allele', 0)))
            except Exception:
                allele = 0
            try:
                ambig = int(float(r.get('is_ambiguous', 0)))
            except Exception:
                ambig = 0
            if rid in selected_ids and ambig == 0 and allele in (-1, 1):
                selected_support[pos] += 1
            hp = hap_tags.get(rid, 'UNK')
            if hp in {'HP1', 'HP2'} and ambig == 0 and allele in (-1, 1):
                tagged_support[pos] += 1
                if hp == 'HP1':
                    hp1_reads[pos] += 1
                else:
                    hp2_reads[pos] += 1
    return selected_support, tagged_support, hp1_reads, hp2_reads


def load_site_debug(path):
    out = {}
    if not path or not Path(path).exists():
        return out
    for r in read_tsv(path):
        try:
            pos = int(r.get('pos', r.get('site_pos')))
        except Exception:
            continue
        def f(k, d=0.0):
            try:
                return float(r.get(k, d) or d)
            except Exception:
                return d
        def i(k, d=0):
            try:
                return int(float(r.get(k, d) or d))
            except Exception:
                return d
        emitted = i('emitted_obs', 0)
        ambiguous = i('ambiguous_obs', 0)
        out[pos] = {
            'ambiguous_frac': ambiguous / emitted if emitted else 0.0,
            'mean_allele_confidence': f('mean_allele_confidence', 0.0),
            'mean_abs_score_delta': f('mean_abs_score_delta', 0.0),
        }
    return out


def bridge_support(a, b, super_support, max_gap_bp):
    if a is None or b is None or b <= a or (b - a) > max_gap_bp:
        return 0.0
    return sum(v for p, v in super_support.items() if a < p < b)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--local-calls-tsv', required=True)
    ap.add_argument('--selected-reads-tsv', required=True)
    ap.add_argument('--haplotag-tsv')
    ap.add_argument('--superreads-tsv')
    ap.add_argument('--obs-tsv', required=True)
    ap.add_argument('--site-debug-tsv')
    ap.add_argument('--site-priors-json')
    ap.add_argument('--chrom', default='chr20')
    ap.add_argument('--min-confidence', type=float, default=0.45)
    ap.add_argument('--min-selected-support', type=int, default=1)
    ap.add_argument('--min-tagged-support', type=int, default=0)
    ap.add_argument('--max-ambiguous-frac', type=float, default=0.30)
    ap.add_argument('--max-gap-bp', type=int, default=50000)
    ap.add_argument('--out-calls-tsv', required=True)
    ap.add_argument('--out-positions-tsv', required=True)
    ap.add_argument('--out-summary-json', required=True)
    args = ap.parse_args()

    local = load_local(args.local_calls_tsv)
    selected_ids = load_selected_ids(args.selected_reads_tsv)
    hap_tags = load_haplotag(args.haplotag_tsv)
    sel, tag, hp1, hp2 = load_obs_support(args.obs_tsv, selected_ids, hap_tags)
    sup = load_super(args.superreads_tsv)
    dbg = load_site_debug(args.site_debug_tsv)
    priors = load_site_priors(args.site_priors_json)

    rows = []
    positions = []
    cur = 0
    prev_pos = None
    prev_state = None

    for lc in local:
        pos = lc['pos']
        d = dbg.get(pos, {})
        local_conf = lc['confidence']
        weighted_conf = local_conf
        weighted_conf += 0.08 * min(sel.get(pos, 0), 3)
        weighted_conf += 0.05 * min(tag.get(pos, 0), 2)
        weighted_conf += 0.03 * min(sup.get(pos, 0.0), 3.0)
        weighted_conf = min(1.0, weighted_conf)

        keep = True
        reasons = []
        if lc['phase_state'] == 0:
            keep = False
            reasons.append('abstain_local')
        if weighted_conf < args.min_confidence:
            keep = False
            reasons.append('low_confidence')
        if sel.get(pos, 0) < args.min_selected_support and (tag.get(pos, 0) + sup.get(pos, 0.0)) < 2.0:
            keep = False
            reasons.append('low_selected_support')
        if tag.get(pos, 0) < args.min_tagged_support:
            keep = False
            reasons.append('low_tagged_support')
        if d.get('ambiguous_frac', 0.0) > args.max_ambiguous_frac:
            keep = False
            reasons.append('too_ambiguous')

        bridge = bridge_support(prev_pos, pos, sup, args.max_gap_bp)
        if keep:
            if prev_pos is None or (pos - prev_pos) > args.max_gap_bp or (prev_state is not None and lc['phase_state'] != prev_state and bridge < 1.0):
                cur += 1
            bid = f'scaffold_{cur}'
            positions.append(pos)
            prev_pos = pos
            prev_state = lc['phase_state']
            decision = 'keep'
        else:
            bid = ''
            decision = '|'.join(reasons) if reasons else 'filtered'

        conf = 0.0
        if keep:
            conf = max(0.0, min(1.0,
                0.35 * weighted_conf +
                0.18 * min(sel.get(pos, 0) / 4.0, 1.0) +
                0.12 * min(tag.get(pos, 0) / 3.0, 1.0) +
                0.12 * min(sup.get(pos, 0.0) / 3.0, 1.0) +
                0.08 * min(d.get('mean_allele_confidence', 0.0) * 4.0, 1.0) +
                0.08 * min(d.get('mean_abs_score_delta', 0.0) / 2.0, 1.0) +
                0.07 * (1.0 - min(abs(priors.get(pos, 0.0)), 1.0))
            ))

        rows.append({
            'chrom': args.chrom,
            'pos': pos,
            'local_phase_state': lc['phase_state'],
            'local_confidence': f'{local_conf:.6f}',
            'weighted_local_confidence': f'{weighted_conf:.6f}',
            'scaffold_phase_state': lc['phase_state'] if keep else 0,
            'scaffold_confidence': f'{conf:.6f}',
            'scaffold_block_id': bid,
            'selected_support': sel.get(pos, 0),
            'tagged_support': tag.get(pos, 0),
            'superread_support': f'{sup.get(pos, 0.0):.3f}',
            'site_prior_bias': f'{priors.get(pos, 0.0):.6f}',
            'ambiguous_frac': f'{d.get("ambiguous_frac", 0.0):.6f}',
            'mean_allele_confidence': f'{d.get("mean_allele_confidence", 0.0):.6f}',
            'mean_abs_score_delta': f'{d.get("mean_abs_score_delta", 0.0):.6f}',
            'hp1_reads': hp1.get(pos, 0),
            'hp2_reads': hp2.get(pos, 0),
            'decision': decision,
        })

    fields = list(rows[0].keys()) if rows else ['chrom', 'pos']
    write_tsv(args.out_calls_tsv, rows, fields)
    write_tsv(args.out_positions_tsv, [{'chrom': args.chrom, 'pos': p} for p in positions], ['chrom', 'pos'])
    summary = {
        'chrom': args.chrom,
        'total_local_sites': len(local),
        'scaffold_sites': len(positions),
        'scaffold_blocks': len({r['scaffold_block_id'] for r in rows if r['scaffold_block_id']}),
    }
    Path(args.out_summary_json).parent.mkdir(parents=True, exist_ok=True)
    json.dump(summary, open(args.out_summary_json, 'w'), indent=2)
    print(json.dumps(summary, indent=2))


if __name__ == '__main__':
    main()
