#!/usr/bin/env python3
"""Process Mengbing VHH/nanobody NGS `*.collapsed.csv` files.

For every `<id>.collapsed.csv` (cols: nt,count) in a folder, this:
  1. translates each nt sequence (frame 0),
  2. segments it into FR1-4 / CDR1-3 against a fixed VHH scaffold (fuzzy anchor search),
  3. filters reads: drops N-containing, stop-codon, and low-quality (unfindable anchor
     or frameshift) reads,
  4. writes `<id>.translated_filtered.csv`,
  5. collapses by CDR amino-acid content -> `<id>.cdr_collapsed.csv`,
  6. writes `QC_summary.csv`.

USAGE
  /Users/jiezhou/miniconda3/envs/tfrc_env/bin/python process_vhh_ngs.py <folder>
  # e.g.
  python process_vhh_ngs.py "/…/NGS_data/20260715_NGS_Mengbing_BQ"

No third-party libraries required (genetic code is hand-coded). If a future library uses a
DIFFERENT scaffold, edit REF_NT / REF_AA / ANCHORS below (see the boundary table in the docstring).

Scaffold reference (constant framework; CDRs diversified):
  EVQLVESGGGLVQPGGSLRLSCAASGYTFTENTMHWFRQAPGKEREWVASIYSSSSYTYYADSVKGRFTISRDNSKNTAYLQMNSLRAEDTAVYYCARAYYGFDYWGQGTLVTVSS
  FR1 1-26 | CDR1 27-35 | FR2 36-49 | CDR2 50-59 | FR3 60-98 | CDR3 99-105 | FR4 106-116
"""
import argparse, csv, glob, os, sys, collections

REF_NT = ("gaggttcagctggtggagtctggcggtggcctggtgcagccagggggatccctccgtttgtcctgtgcagctt"
          "ctggctacacgtttacggaaaacacgatgcactggtttcgtcaggccccgggtaaggaacgggagtgggtcgc"
          "aTCgATTTATTCTTCTTCTAGCTATACTTaTtatgccgatagcgtcaagggccgtttcactatctcgagagac"
          "AACagtaaaaacacagcctacctaCAAATGAACAGCCTGAGAGCCgaggacactgccgtctattattgtgctcg"
          "cGCTTACTACGGTTTTgactactggggtcaaggaaccctgGTCACCGTCTCCTCA").upper()
REF_AA = ("EVQLVESGGGLVQPGGSLRLSCAASGYTFTENTMHWFRQAPGKEREWVASIYSSSSYTYYADSVKGRFTISRDN"
          "SKNTAYLQMNSLRAEDTAVYYCARAYYGFDYWGQGTLVTVSS")

# genetic code
_BASES = "TCAG"
_AAS = "FFLLSSSSYY**CC*WLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG"
CODON = {a+b+c: _AAS[i] for i, (a, b, c) in enumerate(
    (x, y, z) for x in _BASES for y in _BASES for z in _BASES)}

def translate(nt):
    return "".join(CODON.get(nt[j:j+3], "X") for j in range(0, len(nt) - 2, 3))

# conserved framework anchors: (name, motif, ref_nt_start)
ANCHORS = [
    ("A1", "TGTGCAGCTTCTGGC",      63),   # C A A S G (22-26) -> CDR1 after
    ("A2", "TGGTTTCGTCAGGCC",      105),  # W F R Q A (36-40) -> FR2
    ("A3", "AAGGAACGGGAGTGGGTC",   126),  # K E R E W V (43-48); CDR2 = A3_end+3
    ("A4", "TATGCCGATAGCGTCAAGGGC",177),  # Y A D S V K G (60-66) -> FR3
    ("A5", "TATTATTGTGCTCGC",      279),  # Y Y C A R (94-98) -> CDR3 after
    ("A6", "TGGGGTCAAGGAACCCTGG",  315),  # W G Q G T L (106-111) -> FR4
]
_A = {n: m for n, m, _ in ANCHORS}

def _fuzzy_find(read, motif, lo, hi, max_mm):
    L = len(motif); best_pos, best_mm = None, max_mm + 1
    end = min(hi, len(read) - L)
    for s in range(max(0, lo), end + 1):
        mm = 0
        for a, b in zip(read[s:s+L], motif):
            if a != b:
                mm += 1
                if mm >= best_mm:
                    break
        if mm < best_mm:
            best_mm, best_pos = mm, s
            if mm == 0:
                break
    return (best_pos, best_mm) if best_pos is not None else (None, None)

def segment(read):
    if "N" in read:
        return {"fail": "N_in_read"}
    pos = {}; search_from = 0; prev_ref = 0; slack = 30
    for name, motif, ref_off in ANCHORS:
        hi = search_from + (ref_off - prev_ref) + slack + len(motif)
        p, mm = _fuzzy_find(read, motif, search_from, hi, max(1, len(motif) // 5))
        if p is None:
            return {"fail": f"anchor_{name}_not_found"}
        pos[name] = p; search_from = p + len(motif); prev_ref = ref_off + len(motif)
    a1, a2, a3, a4, a5, a6 = (pos[k] for k in ("A1", "A2", "A3", "A4", "A5", "A6"))
    a1e = a1 + len(_A["A1"]); a2e_fr2_start = a2
    a3e = a3 + len(_A["A3"]); a5e = a5 + len(_A["A5"])
    frame = a1 % 3
    fr1, cdr1 = read[frame:a1e], read[a1e:a2]
    fr2, cdr2 = read[a2:a3e+3], read[a3e+3:a4]
    fr3, cdr3 = read[a4:a5e], read[a5e:a6]
    fr4 = read[a6:]
    if any(len(x) % 3 for x in (cdr1, fr2, cdr2, fr3, cdr3)):
        return {"fail": "frameshift"}
    regs = {"FR1": translate(fr1), "CDR1": translate(cdr1), "FR2": translate(fr2),
            "CDR2": translate(cdr2), "FR3": translate(fr3), "CDR3": translate(cdr3),
            "FR4": translate(fr4)}
    vhh = "".join(regs[k] for k in ("FR1", "CDR1", "FR2", "CDR2", "FR3", "CDR3", "FR4"))
    if "*" in vhh:
        return {"fail": "stop_codon"}
    regs["VHH_seq"] = vhh
    return regs

FILT_COLS = ["nt", "count", "VHH_seq", "FR1", "CDR1", "FR2", "CDR2", "FR3", "CDR3", "FR4"]

def process_folder(folder):
    files = sorted(glob.glob(os.path.join(folder, "*.collapsed.csv")))
    if not files:
        sys.exit(f"No *.collapsed.csv files in {folder}")
    summ = []
    for f in files:
        base = os.path.basename(f).replace(".collapsed.csv", "")
        filt_path = os.path.join(folder, base + ".translated_filtered.csv")
        coll_path = os.path.join(folder, base + ".cdr_collapsed.csv")
        n_in = reads_in = kept = kept_reads = 0
        drop = collections.Counter()
        cnt = collections.Counter(); nvar = collections.Counter(); rep = {}
        with open(f) as fh, open(filt_path, "w", newline="") as oh:
            rd = csv.reader(fh); next(rd, None)
            w = csv.writer(oh); w.writerow(FILT_COLS)
            for row in rd:
                if not row:
                    continue
                nt, c = row[0].upper(), int(row[1]); n_in += 1; reads_in += c
                r = segment(nt)
                if "fail" in r:
                    drop[r["fail"]] += 1; continue
                kept += 1; kept_reads += c
                w.writerow([nt, c] + [r[k] for k in FILT_COLS[2:]])
                key = (r["CDR1"], r["CDR2"], r["CDR3"])
                cnt[key] += c; nvar[key] += 1; rep.setdefault(key, r["VHH_seq"])
        with open(coll_path, "w", newline="") as oh:
            w = csv.writer(oh)
            w.writerow(["CDR1", "CDR2", "CDR3", "count", "n_nt_variants", "VHH_seq_example"])
            for key, c in cnt.most_common():
                w.writerow([key[0], key[1], key[2], c, nvar[key], rep[key]])
        lowq = drop.get("frameshift", 0) + sum(v for k, v in drop.items() if k.startswith("anchor"))
        summ.append(dict(sample=base, in_unique=n_in, in_reads=reads_in,
                         kept_unique=kept, kept_reads=kept_reads, cdr_unique=len(cnt),
                         drop_N=drop.get("N_in_read", 0), drop_stop=drop.get("stop_codon", 0),
                         drop_lowqual=lowq))
        print(f"  {base}: kept {kept}/{n_in} unique -> {len(cnt)} unique CDR combos")
    qc = os.path.join(folder, "QC_summary.csv")
    with open(qc, "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=list(summ[0].keys())); w.writeheader()
        w.writerows(summ)
    print(f"Wrote QC_summary.csv and per-sample *.translated_filtered.csv / *.cdr_collapsed.csv")

if __name__ == "__main__":
    ap = argparse.ArgumentParser(description="Translate + segment + filter + CDR-collapse VHH NGS collapsed.csv files.")
    ap.add_argument("folder", help="folder containing *.collapsed.csv files")
    args = ap.parse_args()
    # sanity: reference must translate to itself
    assert segment(REF_NT)["VHH_seq"] == REF_AA, "reference scaffold mis-specified"
    process_folder(args.folder)
