#!/usr/bin/env python3
import argparse, csv, json
from pathlib import Path
import pysam


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


def load_variants(path):
    obj = json.load(open(path))
    obj = obj.get('variants', obj) if isinstance(obj, dict) else obj
    out = {}
    for r in obj:
        try:
            out[int(r['pos'])] = (str(r['ref_allele']), str(r['alt_allele']))
        except Exception:
            pass
    return out


def resolve_sample(vcf, requested=None):
    samples = list(vcf.header.samples)
    if not samples:
        raise SystemExit('Truth VCF has no samples')
    if requested and requested in samples:
        return requested
    if requested and requested not in samples:
        if len(samples) == 1:
            return samples[0]
        raise SystemExit(f"Requested truth sample '{requested}' not found. Available samples: {samples}")
    return samples[0]


def load_truth(path, sample=None):
    vcf = pysam.VariantFile(path)
    sample = resolve_sample(vcf, sample)
    truth = {}
    for rec in vcf.fetch():
        if rec.alts is None or len(rec.alts) != 1:
            continue
        sd = rec.samples[sample]
        gt = sd.get('GT')
        phased = getattr(sd, 'phased', False)
        if not phased or gt is None or len(gt) != 2 or None in gt:
            continue
        if tuple(gt) == (0, 1):
            state = 1
        elif tuple(gt) == (1, 0):
            state = -1
        else:
            continue
        truth[(int(rec.pos), str(rec.ref), str(rec.alts[0]))] = state
    return truth, sample


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--candidate-features-tsv', required=True)
    ap.add_argument('--variant-json', required=True)
    ap.add_argument('--truth-vcf', required=True)
    ap.add_argument('--sample')
    ap.add_argument('--out-tsv', required=True)
    ap.add_argument('--out-summary-json', required=True)
    args = ap.parse_args()

    feats = read_tsv(args.candidate_features_tsv)
    var = load_variants(args.variant_json)
    truth, sample = load_truth(args.truth_vcf, args.sample)

    rows = []
    labeled = 0
    pos_rows = 0
    for r in feats:
        pos = int(r['pos'])
        ref, alt = var.get(pos, ('', ''))
        t = truth.get((pos, ref, alt))
        row = dict(r)
        row['truth_state'] = '' if t is None else t
        if t is None:
            row['label'] = ''
        else:
            row['label'] = 1 if int(r['candidate_phase_state']) == t else 0
            labeled += 1
            pos_rows += int(row['label'] == 1)
        rows.append(row)

    fields = list(rows[0].keys()) if rows else []
    write_tsv(args.out_tsv, rows, fields)
    summary = {'rows': len(rows), 'labeled_rows': labeled, 'positive_rows': pos_rows, 'truth_sample': sample}
    Path(args.out_summary_json).parent.mkdir(parents=True, exist_ok=True)
    json.dump(summary, open(args.out_summary_json, 'w'), indent=2)
    print(json.dumps(summary, indent=2))


if __name__ == '__main__':
    main()
