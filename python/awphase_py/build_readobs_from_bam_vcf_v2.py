#!/usr/bin/env python3
import argparse
import csv
import json
from collections import Counter, defaultdict
from pathlib import Path

import pysam


def normalize_chrom_name(bam, chrom: str) -> str:
    refs = set(bam.references)
    if chrom in refs:
        return chrom
    if chrom.startswith("chr") and chrom[3:] in refs:
        return chrom[3:]
    if ("chr" + chrom) in refs:
        return "chr" + chrom
    raise ValueError(f"Chromosome {chrom!r} not found in BAM header")


def gt_is_diploid_het(gt):
    if gt is None or len(gt) != 2:
        return False
    a, b = gt
    if a is None or b is None:
        return False
    return {a, b} == {0, 1}


def softclip_lengths(aln):
    left = 0
    right = 0
    cigar = aln.cigartuples or []
    if cigar and cigar[0][0] == 4:
        left = cigar[0][1]
    if cigar and cigar[-1][0] == 4:
        right = cigar[-1][1]
    return left, right


def nearest_neighbor_dist(sorted_positions):
    out = {}
    n = len(sorted_positions)
    for i, pos in enumerate(sorted_positions):
        left = abs(pos - sorted_positions[i - 1]) if i > 0 else 10**9
        right = abs(sorted_positions[i + 1] - pos) if i + 1 < n else 10**9
        out[pos] = min(left, right)
    return out


def penalty_and_flags(
    baseq,
    mapq,
    nm,
    dist_left,
    dist_right,
    softclip_left,
    softclip_right,
    neighbor_dist,
    args,
):
    penalty = 1.0
    flags = []

    near_end = min(dist_left, dist_right) <= args.end_buffer
    if near_end:
        penalty *= args.end_penalty
        flags.append("near_read_end")

    clipped_near_site = (
        (softclip_left > 0 and dist_left <= args.clip_buffer) or
        (softclip_right > 0 and dist_right <= args.clip_buffer)
    )
    if clipped_near_site:
        penalty *= args.clip_penalty
        flags.append("softclip_near_site")

    if nm is not None and nm >= args.nm_threshold:
        penalty *= args.nm_penalty
        flags.append("high_nm")

    if neighbor_dist <= args.context_bp:
        penalty *= args.context_penalty
        flags.append("dense_variant_context")

    eff_baseq = max(1, int(round(baseq * penalty)))
    return eff_baseq, clipped_near_site, near_end, penalty, ",".join(flags)


def score_context_window(query_seq, query_qual, center_qpos, ref_ctx, alt_ctx, radius):
    if query_seq is None or center_qpos is None:
        return 0.0, 0.0, 0, False, 0.0
    left = center_qpos - radius
    right = center_qpos + radius + 1
    if left < 0 or right > len(query_seq):
        return 0.0, 0.0, 0, False, 0.0

    qctx = query_seq[left:right].upper()
    quals = []
    if query_qual is not None:
        quals = list(query_qual[left:right])
    if len(qctx) != len(ref_ctx) or len(qctx) != len(alt_ctx):
        return 0.0, 0.0, 0, False, 0.0

    ref_score = 0.0
    alt_score = 0.0
    compared = 0
    for i, qbase in enumerate(qctx):
        if qbase not in {"A", "C", "G", "T"}:
            continue
        q = quals[i] if i < len(quals) else 20
        w = max(0.1, min(1.0, q / 40.0))
        ref_score += w if qbase == ref_ctx[i] else -0.35 * w
        alt_score += w if qbase == alt_ctx[i] else -0.35 * w
        compared += 1

    delta = alt_score - ref_score
    denom = abs(ref_score) + abs(alt_score) + 1e-6
    conf = min(1.0, abs(delta) / denom) if compared else 0.0
    ambiguous = compared < 3 or abs(delta) < 0.20
    return ref_score, alt_score, compared, ambiguous, conf


def local_context_scores(fasta, chrom, pos, ref, alt, aln, qpos, radius):
    if radius <= 0:
        return 0.0, 0.0, 0, True, 0.0
    start0 = max(0, pos - 1 - radius)
    end0 = pos + radius
    ref_ctx = fasta.fetch(chrom, start0, end0).upper()
    if len(ref_ctx) != (2 * radius + 1):
        return 0.0, 0.0, 0, True, 0.0
    if ref_ctx[radius] != ref:
        ref_ctx = ref_ctx[:radius] + ref + ref_ctx[radius + 1:]
    alt_ctx = ref_ctx[:radius] + alt + ref_ctx[radius + 1:]
    return score_context_window(
        aln.query_sequence,
        aln.query_qualities,
        qpos,
        ref_ctx,
        alt_ctx,
        radius,
    )


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--bam", required=True)
    ap.add_argument("--vcf", required=True)
    ap.add_argument("--fasta", required=True)
    ap.add_argument("--region", required=True, help="e.g. chr20:1-5000000")
    ap.add_argument("--sample", default="HG002")
    ap.add_argument("--min-mapq", type=int, default=20)
    ap.add_argument("--min-baseq", type=int, default=20)
    ap.add_argument("--end-buffer", type=int, default=10)
    ap.add_argument("--clip-buffer", type=int, default=12)
    ap.add_argument("--context-bp", type=int, default=5)
    ap.add_argument("--end-penalty", type=float, default=0.65)
    ap.add_argument("--clip-penalty", type=float, default=0.50)
    ap.add_argument("--context-penalty", type=float, default=0.80)
    ap.add_argument("--nm-threshold", type=int, default=5)
    ap.add_argument("--nm-penalty", type=float, default=0.80)
    ap.add_argument("--out-json", required=True)
    ap.add_argument("--out-summary", required=True)
    ap.add_argument("--out-site-tsv", required=True)
    ap.add_argument("--out-read-tsv", required=True)
    ap.add_argument("--out-obs-tsv", required=True)
    args = ap.parse_args()

    chrom, coords = args.region.split(":")
    region_start, region_end = [int(x.replace(",", "")) for x in coords.split("-")]

    bam = pysam.AlignmentFile(args.bam, "rb")
    vcf = pysam.VariantFile(args.vcf)
    fasta = pysam.FastaFile(args.fasta)
    bam_chrom = normalize_chrom_name(bam, chrom)

    total_records = 0
    records = []
    skip_counter = Counter()

    for rec in vcf.fetch(chrom, region_start - 1, region_end):
        total_records += 1
        row = {
            "pos": rec.pos,
            "ref": rec.ref.upper(),
            "alts": [a.upper() for a in (rec.alts or [])],
            "sample_present": args.sample in rec.samples,
            "gt": None,
            "usable": False,
            "skip_reason": "",
        }

        if len(row["alts"]) != 1:
            row["skip_reason"] = "non_biallelic"
            skip_counter[row["skip_reason"]] += 1
            records.append(row)
            continue

        alt = row["alts"][0]
        if len(row["ref"]) != 1 or len(alt) != 1:
            row["skip_reason"] = "non_snv"
            skip_counter[row["skip_reason"]] += 1
            records.append(row)
            continue

        if not row["sample_present"]:
            row["skip_reason"] = "missing_sample"
            skip_counter[row["skip_reason"]] += 1
            records.append(row)
            continue

        gt = rec.samples[args.sample].get("GT")
        row["gt"] = gt
        if not gt_is_diploid_het(gt):
            row["skip_reason"] = "non_het"
            skip_counter[row["skip_reason"]] += 1
            records.append(row)
            continue

        fasta_ref = fasta.fetch(chrom, rec.pos - 1, rec.pos).upper()
        if fasta_ref != row["ref"]:
            row["skip_reason"] = "ref_mismatch_vs_fasta"
            skip_counter[row["skip_reason"]] += 1
            records.append(row)
            continue

        row["usable"] = True
        records.append(row)

    usable_positions = sorted(r["pos"] for r in records if r["usable"])
    neighbor_dist = nearest_neighbor_dist(usable_positions)

    observations = []
    obs_rows = []
    site_rows = []
    read_aggr = defaultdict(lambda: {
        "obs_count": 0,
        "ref_obs": 0,
        "alt_obs": 0,
        "sum_signed_weight": 0.0,
        "sum_weight": 0.0,
        "reverse_obs": 0,
        "forward_obs": 0,
        "sites": [],
    })

    for row in records:
        if not row["usable"]:
            site_rows.append({
                "pos": row["pos"],
                "ref": row["ref"],
                "alt": ",".join(row["alts"]),
                "usable": False,
                "skip_reason": row["skip_reason"],
                "neighbor_dist": "",
                "reads_seen": 0,
                "ref_obs": 0,
                "alt_obs": 0,
                "other_base": 0,
                "deletions": 0,
                "refskips": 0,
                "low_mapq": 0,
                "low_baseq": 0,
                "secondary_or_supplementary": 0,
                "duplicates": 0,
                "emitted_obs": 0,
                "ambiguous_obs": 0,
                "alt_like_context_obs": 0,
                "ref_like_context_obs": 0,
                "mean_allele_confidence": 0.0,
                "mean_signed_score_delta": 0.0,
                "mean_abs_score_delta": 0.0,
                "mean_local_ref_score": 0.0,
                "mean_local_alt_score": 0.0,
            })
            continue

        pos = row["pos"]
        ref = row["ref"]
        alt = row["alts"][0]

        counts = Counter()
        emitted_obs = 0
        sum_allele_conf = 0.0
        sum_signed_delta = 0.0
        sum_abs_delta = 0.0
        sum_local_ref_score = 0.0
        sum_local_alt_score = 0.0

        for col in bam.pileup(
            bam_chrom,
            pos - 1,
            pos,
            truncate=True,
            min_base_quality=0,
            stepper="samtools",
            fastafile=fasta,
            ignore_overlaps=False,
            ignore_orphans=False,
        ):
            if col.reference_pos != pos - 1:
                continue

            for pr in col.pileups:
                counts["reads_seen"] += 1

                if pr.is_del:
                    counts["deletions"] += 1
                    continue
                if pr.is_refskip:
                    counts["refskips"] += 1
                    continue

                aln = pr.alignment

                if aln.is_secondary or aln.is_supplementary:
                    counts["secondary_or_supplementary"] += 1
                    continue
                if aln.is_duplicate:
                    counts["duplicates"] += 1
                    continue
                if aln.mapping_quality < args.min_mapq:
                    counts["low_mapq"] += 1
                    continue

                qpos = pr.query_position
                if qpos is None:
                    counts["null_query_position"] += 1
                    continue

                baseq = int(aln.query_qualities[qpos]) if aln.query_qualities is not None else 0
                if baseq < args.min_baseq:
                    counts["low_baseq"] += 1
                    continue

                base = aln.query_sequence[qpos].upper()

                if base == ref:
                    allele = -1
                    counts["ref_obs"] += 1
                elif base == alt:
                    allele = 1
                    counts["alt_obs"] += 1
                else:
                    counts["other_base"] += 1
                    continue

                softclip_left, softclip_right = softclip_lengths(aln)
                read_len = aln.query_length or len(aln.query_sequence)
                dist_left = qpos
                dist_right = read_len - qpos - 1
                nm = aln.get_tag("NM") if aln.has_tag("NM") else None

                eff_baseq, clipped_near_site, near_end, penalty, flags = penalty_and_flags(
                    baseq=baseq,
                    mapq=int(aln.mapping_quality),
                    nm=nm,
                    dist_left=dist_left,
                    dist_right=dist_right,
                    softclip_left=softclip_left,
                    softclip_right=softclip_right,
                    neighbor_dist=neighbor_dist[pos],
                    args=args,
                )

                weight = (eff_baseq / 40.0) * (int(aln.mapping_quality) / 60.0)

                radius = min(args.context_bp, qpos, read_len - qpos - 1)
                local_ref_score, local_alt_score, context_bases_compared, is_ambiguous, allele_confidence = local_context_scores(
                    fasta=fasta,
                    chrom=chrom,
                    pos=pos,
                    ref=ref,
                    alt=alt,
                    aln=aln,
                    qpos=qpos,
                    radius=radius,
                )
                allele_score_delta = local_alt_score - local_ref_score

                weight_v3 = weight
                if is_ambiguous:
                    weight_v3 *= 0.5
                    counts["ambiguous_obs"] += 1
                else:
                    if allele_score_delta > 0:
                        counts["alt_like_context_obs"] += 1
                    elif allele_score_delta < 0:
                        counts["ref_like_context_obs"] += 1
                weight_v3 *= (0.75 + 0.25 * allele_confidence)
                weight_v3 = max(0.01, weight_v3)

                sum_allele_conf += allele_confidence
                sum_signed_delta += allele_score_delta
                sum_abs_delta += abs(allele_score_delta)
                sum_local_ref_score += local_ref_score
                sum_local_alt_score += local_alt_score

                obs = {
                    "read_id": aln.query_name,
                    "site_pos": pos,
                    "allele": allele,
                    "baseq": eff_baseq,
                    "mapq": int(aln.mapping_quality),
                    "raw_baseq": baseq,
                    "weight_v2": weight,
                    "weight_v3": weight_v3,
                    "local_ref_score": local_ref_score,
                    "local_alt_score": local_alt_score,
                    "allele_score_delta": allele_score_delta,
                    "allele_confidence": allele_confidence,
                    "is_ambiguous": is_ambiguous,
                    "context_bases_compared": context_bases_compared,
                    "variant_class": "snv",
                    "penalty_multiplier": penalty,
                }
                observations.append(obs)
                emitted_obs += 1

                obs_rows.append({
                    "read_id": aln.query_name,
                    "site_pos": pos,
                    "ref": ref,
                    "alt": alt,
                    "allele": allele,
                    "raw_baseq": baseq,
                    "effective_baseq": eff_baseq,
                    "mapq": int(aln.mapping_quality),
                    "nm": nm if nm is not None else "",
                    "is_reverse": int(aln.is_reverse),
                    "is_read1": int(aln.is_read1),
                    "is_read2": int(aln.is_read2),
                    "is_proper_pair": int(aln.is_paired and aln.is_proper_pair),
                    "softclip_left": softclip_left,
                    "softclip_right": softclip_right,
                    "distance_to_left_end": dist_left,
                    "distance_to_right_end": dist_right,
                    "near_end": int(near_end),
                    "clipped_near_site": int(clipped_near_site),
                    "neighbor_variant_dist": neighbor_dist[pos],
                    "local_ref_score": f"{local_ref_score:.6f}",
                    "local_alt_score": f"{local_alt_score:.6f}",
                    "allele_score_delta": f"{allele_score_delta:.6f}",
                    "allele_confidence": f"{allele_confidence:.6f}",
                    "is_ambiguous": int(is_ambiguous),
                    "context_bases_compared": context_bases_compared,
                    "variant_class": "snv",
                    "penalty_multiplier": f"{penalty:.4f}",
                    "flags": flags,
                    "weight_v2": f"{weight:.6f}",
                    "weight_v3": f"{weight_v3:.6f}",
                })

                ra = read_aggr[aln.query_name]
                ra["obs_count"] += 1
                ra["sum_signed_weight"] += allele * weight_v3
                ra["sum_weight"] += abs(weight_v3)
                ra["sites"].append(pos)
                if allele == -1:
                    ra["ref_obs"] += 1
                elif allele == 1:
                    ra["alt_obs"] += 1
                if aln.is_reverse:
                    ra["reverse_obs"] += 1
                else:
                    ra["forward_obs"] += 1

        site_rows.append({
            "pos": pos,
            "ref": ref,
            "alt": alt,
            "usable": True,
            "skip_reason": "",
            "neighbor_dist": neighbor_dist[pos],
            "reads_seen": counts["reads_seen"],
            "ref_obs": counts["ref_obs"],
            "alt_obs": counts["alt_obs"],
            "other_base": counts["other_base"],
            "deletions": counts["deletions"],
            "refskips": counts["refskips"],
            "low_mapq": counts["low_mapq"],
            "low_baseq": counts["low_baseq"],
            "secondary_or_supplementary": counts["secondary_or_supplementary"],
            "duplicates": counts["duplicates"],
            "emitted_obs": emitted_obs,
            "ambiguous_obs": counts["ambiguous_obs"],
            "alt_like_context_obs": counts["alt_like_context_obs"],
            "ref_like_context_obs": counts["ref_like_context_obs"],
            "mean_allele_confidence": f"{(sum_allele_conf / emitted_obs) if emitted_obs else 0.0:.6f}",
            "mean_signed_score_delta": f"{(sum_signed_delta / emitted_obs) if emitted_obs else 0.0:.6f}",
            "mean_abs_score_delta": f"{(sum_abs_delta / emitted_obs) if emitted_obs else 0.0:.6f}",
            "mean_local_ref_score": f"{(sum_local_ref_score / emitted_obs) if emitted_obs else 0.0:.6f}",
            "mean_local_alt_score": f"{(sum_local_alt_score / emitted_obs) if emitted_obs else 0.0:.6f}",
        })

    read_rows = []
    for read_id, ra in read_aggr.items():
        signed = ra["sum_signed_weight"]
        if signed > 0.2:
            tag = "hap_alt_like"
        elif signed < -0.2:
            tag = "hap_ref_like"
        else:
            tag = "ambiguous"

        read_rows.append({
            "read_id": read_id,
            "obs_count": ra["obs_count"],
            "ref_obs": ra["ref_obs"],
            "alt_obs": ra["alt_obs"],
            "forward_obs": ra["forward_obs"],
            "reverse_obs": ra["reverse_obs"],
            "sum_weight": f"{ra['sum_weight']:.6f}",
            "sum_signed_weight": f"{signed:.6f}",
            "pseudo_haplotag": tag,
            "sites": ",".join(map(str, sorted(set(ra["sites"]))[:100])),
        })

    summary = {
        "region": args.region,
        "sample": args.sample,
        "bam": args.bam,
        "vcf": args.vcf,
        "fasta": args.fasta,
        "min_mapq": args.min_mapq,
        "min_baseq": args.min_baseq,
        "total_vcf_records_seen": total_records,
        "usable_het_biallelic_snvs": len(usable_positions),
        "skip_counts": dict(skip_counter),
        "observations_emitted": len(observations),
        "informative_reads": len(read_aggr),
        "sites_with_observations": sum(1 for r in site_rows if r["usable"] and int(r["emitted_obs"]) > 0),
        "mean_observations_per_observed_site": (
            sum(int(r["emitted_obs"]) for r in site_rows if r["usable"] and int(r["emitted_obs"]) > 0)
            / max(1, sum(1 for r in site_rows if r["usable"] and int(r["emitted_obs"]) > 0))
        ),
    }

    Path(args.out_json).parent.mkdir(parents=True, exist_ok=True)
    with open(args.out_json, "w") as fh:
        json.dump(observations, fh)

    with open(args.out_summary, "w") as fh:
        json.dump(summary, fh, indent=2)

    with open(args.out_site_tsv, "w", newline="") as fh:
        writer = csv.DictWriter(fh, fieldnames=list(site_rows[0].keys()), delimiter="\t")
        writer.writeheader()
        writer.writerows(site_rows)

    with open(args.out_obs_tsv, "w", newline="") as fh:
        writer = csv.DictWriter(fh, fieldnames=list(obs_rows[0].keys()) if obs_rows else [
            "read_id","site_pos","ref","alt","allele","raw_baseq","effective_baseq","mapq","nm",
            "is_reverse","is_read1","is_read2","is_proper_pair","softclip_left","softclip_right",
            "distance_to_left_end","distance_to_right_end","near_end","clipped_near_site",
            "neighbor_variant_dist","local_ref_score","local_alt_score","allele_score_delta","allele_confidence",
            "is_ambiguous","context_bases_compared","variant_class","penalty_multiplier","flags","weight_v2","weight_v3"
        ], delimiter="\t")
        writer.writeheader()
        if obs_rows:
            writer.writerows(obs_rows)

    with open(args.out_read_tsv, "w", newline="") as fh:
        writer = csv.DictWriter(fh, fieldnames=list(read_rows[0].keys()) if read_rows else [
            "read_id","obs_count","ref_obs","alt_obs","forward_obs","reverse_obs",
            "sum_weight","sum_signed_weight","pseudo_haplotag","sites"
        ], delimiter="\t")
        writer.writeheader()
        if read_rows:
            writer.writerows(read_rows)

    print(json.dumps(summary, indent=2))


if __name__ == "__main__":
    main()
