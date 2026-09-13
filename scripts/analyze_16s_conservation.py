#!/usr/bin/env python3
"""Measure per-position conservation of the 16S rRNA gene and plot the result.

Pipeline:
  1. Fetch the E. coli 16S reference (J01695.2 region 1268..2809 = exactly 1542 bp).
     This is the classical coordinate system used to locate the V1-V9 regions.
  2. For every sequence in the input FASTA files: fix the strand if the record was
     deposited in reverse orientation, align it to the reference with edlib
     (infix alignment), and project each base onto reference coordinates.
  3. Count A/C/G/T/ambiguous/gap at every position 1..1542. The denominator counts
     only sequences that actually cover a position, so the two ends of the gene are
     not distorted by the many partial records.
  4. Compute % identity (frequency of the most common base) and standardised
     Shannon entropy, smooth them with a sliding window, and draw a 4-panel figure.

Examples:
    python3 analyze_16s_conservation.py --limit 2000   # quick trial run
    python3 analyze_16s_conservation.py                # full dataset
"""

import argparse
import math
import os
import sys

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
DATA = os.environ.get("SIXTEENS_DATA_DIR", HERE)
OUT = os.environ.get("SIXTEENS_OUT_DIR", HERE)

FASTA_FILES = [os.path.join(DATA, "16s_sequences.fasta"),
               os.path.join(DATA, "16s_from_genomes.fasta")]
REF_FASTA = os.path.join(OUT, "ecoli_16s_reference.fasta")
OUT_TSV = os.path.join(OUT, "16s_position_stats.tsv")
OUT_EXCLUDED = os.path.join(OUT, "16s_excluded.tsv")
OUT_IDENTITY = os.path.join(OUT, "16s_alignment_identity.tsv")
OUT_PNG = os.path.join(OUT, "16s_conservation.png")
OUT_PDF = os.path.join(OUT, "16s_conservation.pdf")

REF_ACC = "J01695.2"
REF_START, REF_STOP = 1268, 2809   # NCBI's own annotation: "16S rRNA (rrnB)"
REF_LEN = REF_STOP - REF_START + 1  # 1542

BASES = "ACGT"
# Conserved motifs used to decide which strand a record was deposited on
ORIENT_MOTIFS = ("GTGCCAGCAGCCGCGGTAA", "GGATTAGATACCC", "AAACTCAAAGGAATTGACGG")
COMPLEMENT = str.maketrans("ACGTUNRYKMSWBDHVacgtunrykmswbdhv",
                           "TGCAANYRMKSWVHDBtgcaanyrmkswvhdb")

# Hypervariable regions in E. coli numbering (Chakravorty et al. 2007)
V_REGIONS = [("V1", 69, 99), ("V2", 137, 242), ("V3", 433, 497), ("V4", 576, 682),
             ("V5", 822, 879), ("V6", 986, 1043), ("V7", 1117, 1173),
             ("V8", 1243, 1294), ("V9", 1435, 1465)]
# Commonly sequenced amplicons (span between the usual primer pairs)
AMPLICONS = [("V1-V2", 27, 338), ("V1-V3", 27, 534), ("V3-V5", 341, 926),
             ("V4", 515, 806), ("V6-V9", 968, 1492), ("V1-V9", 27, 1492)]


def revcomp(seq):
    return seq.translate(COMPLEMENT)[::-1]


def read_fasta(path):
    name, chunks = None, []
    with open(path, encoding="utf-8") as fh:
        for line in fh:
            line = line.rstrip()
            if line.startswith(">"):
                if name:
                    yield name, "".join(chunks).upper()
                name, chunks = line[1:], []
            elif name is not None:
                chunks.append(line)
    if name:
        yield name, "".join(chunks).upper()


def get_reference():
    """Return the reference 16S sequence, downloading it once and caching it."""
    if os.path.exists(REF_FASTA):
        seq = next(read_fasta(REF_FASTA))[1]
        if len(seq) == REF_LEN:
            return seq
    from ncbi_common import eutils, make_session
    print("Downloading reference %s:%d-%d from NCBI..." % (REF_ACC, REF_START, REF_STOP),
          flush=True)
    text = eutils(make_session(), "efetch",
                  {"db": "nuccore", "id": REF_ACC, "rettype": "fasta", "retmode": "text",
                   "seq_start": REF_START, "seq_stop": REF_STOP})
    with open(REF_FASTA, "w", encoding="utf-8") as fh:
        fh.write(text if text.endswith("\n") else text + "\n")
    seq = next(read_fasta(REF_FASTA))[1]
    if len(seq) != REF_LEN:
        sys.exit("Reference is %d bp, expected %d bp" % (len(seq), REF_LEN))
    return seq


def fix_orientation(seq):
    """Reverse-complement the record if it sits on the opposite strand.

    Returns (sequence, was_flipped).
    """
    rc = revcomp(seq)
    fwd = sum(m in seq for m in ORIENT_MOTIFS)
    rev = sum(m in rc for m in ORIENT_MOTIFS)
    return (rc, True) if rev > fwd else (seq, False)


def map_to_reference(seq, ref, edlib):
    """Project every base of seq onto reference coordinates.

    Returns (int8 array of length REF_LEN, n_matching_bases, n_readable_bases).
    Each array element holds the ASCII code of the base at that reference position;
    0 means the sequence does not cover the position, and ord('-') means the
    sequence covers it but lacks a base there (a deletion relative to the reference).

    "n_readable_bases" counts only positions carrying a real A/C/G/T. Ambiguity
    codes (N, R, Y, ...) are excluded from both the numerator and the denominator of
    the identity, so an old low-quality record full of Ns is not wrongly judged to be
    something other than 16S.
    """
    # edlib's "HW" mode requires the query to fit inside the target, so the shorter
    # sequence has to be the query
    swapped = len(seq) > len(ref)
    query, target = (ref, seq) if swapped else (seq, ref)
    res = edlib.align(query, target, mode="HW", task="path")
    if res["editDistance"] < 0 or not res["locations"]:
        return None, 0, 0

    t_start = res["locations"][0][0]  # start of the alignment on the target
    out = np.zeros(REF_LEN, dtype=np.int8)
    qi, ti = 0, t_start
    n_match = 0

    for num, op in _cigar_ops(res["cigar"]):
        for _ in range(num):
            # which index refers to the reference depends on whether roles were swapped
            if swapped:
                ref_i, seq_i = qi, ti
            else:
                ref_i, seq_i = ti, qi
            if op in "=XM":
                if 0 <= ref_i < REF_LEN:
                    out[ref_i] = ord(seq[seq_i])
                    if seq[seq_i] == ref[ref_i]:
                        n_match += 1
                qi += 1
                ti += 1
            elif op == "I":      # present in the query, absent from the target
                if swapped and 0 <= qi < REF_LEN:
                    out[qi] = ord("-")   # reference has a base here, the sequence does not
                qi += 1
            elif op == "D":      # present in the target, absent from the query
                if not swapped and 0 <= ti < REF_LEN:
                    out[ti] = ord("-")
                ti += 1
    n_known = int(np.isin(out, [ord(b) for b in BASES]).sum())
    return out, n_match, n_known


def _cigar_ops(cigar):
    num = ""
    for ch in cigar or "":
        if ch.isdigit():
            num += ch
        else:
            yield int(num or 1), ch
            num = ""


def smooth(values, window):
    """Centred moving average that ignores NaN values."""
    if window <= 1:
        return values
    v = np.asarray(values, dtype=float)
    ok = np.isfinite(v)
    filled = np.where(ok, v, 0.0)
    kernel = np.ones(window)
    num = np.convolve(filled, kernel, mode="same")
    den = np.convolve(ok.astype(float), kernel, mode="same")
    with np.errstate(invalid="ignore", divide="ignore"):
        return np.where(den > 0, num / den, np.nan)


def build_matrix(ref, args):
    import edlib

    counts = np.zeros((REF_LEN, 6), dtype=np.int64)  # A C G T ambiguous gap
    idx = {b: i for i, b in enumerate(BASES)}
    excluded = []
    per_seq = []
    n_used = n_flipped = 0
    n_total = 0

    for path in FASTA_FILES:
        for header, seq in read_fasta(path):
            n_total += 1
            if args.limit and n_total > args.limit:
                break
            acc = header.split()[0]
            if len(seq) < 100:
                excluded.append((acc, "sequence too short (%d bp)" % len(seq)))
                continue
            seq, flipped = fix_orientation(seq)
            n_flipped += flipped

            mapped, n_match, n_cov = map_to_reference(seq, ref, edlib)
            if mapped is None or n_cov == 0:
                excluded.append((acc, "could not be aligned to the reference"))
                continue
            if n_cov < args.min_coverage:
                excluded.append((acc, "only %d readable bases on the reference frame (<%d)"
                                 % (n_cov, args.min_coverage)))
                continue
            identity = n_match / n_cov
            if identity < args.min_identity:
                excluded.append((acc, "identity only %.1f%% (<%.0f%%)"
                                 % (identity * 100, args.min_identity * 100)))
                continue

            per_seq.append((acc, identity * 100, n_cov, header))
            for b, col in idx.items():
                counts[mapped == ord(b), col] += 1
            covered = mapped != 0
            known = covered & np.isin(mapped, [ord(b) for b in BASES])
            counts[covered & (mapped == ord("-")), 5] += 1
            counts[covered & ~known & (mapped != ord("-")), 4] += 1
            n_used += 1
        if args.limit and n_total > args.limit:
            break

    return counts, n_used, n_flipped, excluded, per_seq


def compute_stats(counts):
    total = counts.sum(axis=1).astype(float)          # sequences covering the position
    acgt = counts[:, :4].astype(float)
    acgt_sum = acgt.sum(axis=1)

    with np.errstate(invalid="ignore", divide="ignore"):
        pct = np.where(total[:, None] > 0, counts / total[:, None] * 100, np.nan)
        freq = np.where(acgt_sum[:, None] > 0, acgt / acgt_sum[:, None], np.nan)
        identity = np.where(acgt_sum > 0, acgt.max(axis=1) / acgt_sum * 100, np.nan)
        entropy = -np.nansum(np.where(freq > 0, freq * np.log2(freq), 0.0), axis=1)
        entropy = np.where(acgt_sum > 0, entropy, np.nan)
    emax = np.nanmax(entropy)
    entropy_std = entropy / emax if emax > 0 else entropy
    return total, pct, identity, entropy, entropy_std


def write_tsv(path, total, pct, identity, entropy, entropy_std):
    with open(path, "w", encoding="utf-8") as fh:
        fh.write("position\tn_covered\tpct_A\tpct_C\tpct_G\tpct_T\tpct_ambiguous\t"
                 "pct_gap\tpct_identity\tentropy_bits\tentropy_std\n")
        for i in range(REF_LEN):
            fh.write("%d\t%d\t%s\t%s\n"
                     % (i + 1, int(total[i]),
                        "\t".join("%.4f" % v for v in pct[i]),
                        "\t".join("%.4f" % v for v in (identity[i], entropy[i], entropy_std[i]))))


def plot(total, pct, identity, entropy_std, n_used, window, png, pdf):
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    from matplotlib.colors import LinearSegmentedColormap

    x = np.arange(1, REF_LEN + 1)
    ident_s = smooth(identity, window)
    ent_s = smooth(entropy_std, window)

    fig, axes = plt.subplots(4, 1, figsize=(12, 12.5), sharex=True,
                             gridspec_kw={"height_ratios": [1.0, 1.15, 0.62, 0.34],
                                          "hspace": 0.28})

    def shade(ax, label_on_top=False):
        for name, a, b in V_REGIONS:
            ax.axvspan(a, b, color="0.88", zorder=0)
            if label_on_top:
                ax.text((a + b) / 2, 1.02, name, transform=ax.get_xaxis_transform(),
                        ha="center", va="bottom", fontsize=9)

    # --- Panel A: % sequence identity ---
    ax = axes[0]
    shade(ax)
    ax.plot(x, ident_s, color="#e8262f", lw=1.1)
    ax.set_ylim(60, 102)
    ax.set_yticks([60, 65, 70, 75, 80, 85, 90, 95, 100])
    ax.set_ylabel("% sequence identity")
    ax.set_title("16S rRNA gene conservation — %s LPSN sequences mapped to E. coli %s "
                 "(%d bp sliding window)" % (format(n_used, ","), REF_ACC, window),
                 fontsize=12, pad=14)
    valid = np.isfinite(ident_s)
    cons = int(np.nanargmax(np.where(valid, ident_s, -np.inf)))
    hyper = int(np.nanargmin(np.where(valid, ident_s, np.inf)))
    ax.annotate("Conserved region", xy=(x[cons], ident_s[cons]),
                xytext=(x[cons] + 130, 101), ha="center", fontsize=9,
                arrowprops=dict(arrowstyle="->", lw=0.9, color="0.25"))
    ax.annotate("Hypervariable region", xy=(x[hyper], ident_s[hyper]),
                xytext=(x[hyper] + 190, 63.5), ha="center", fontsize=9,
                arrowprops=dict(arrowstyle="->", lw=0.9, color="0.25"))

    # --- Panel B: standardised entropy plus amplicon spans ---
    ax = axes[1]
    shade(ax, label_on_top=True)
    ax.plot(x, ent_s, color="black", lw=1.1)
    ax.set_ylim(-0.62, 1.05)
    ax.set_yticks([0, 0.25, 0.5, 0.75, 1.0])
    ax.set_ylabel("Standardized entropy")
    ax.axhline(0, color="0.6", lw=0.6)
    for k, (name, a, b) in enumerate(AMPLICONS):
        y = -0.10 - 0.085 * k
        ax.plot([a, b], [y, y], color="#8b2c2c", lw=1.8, solid_capstyle="butt",
                clip_on=False)
        ax.text((a + b) / 2, y + 0.018, name, ha="center", va="bottom",
                fontsize=8, color="#8b2c2c")

    # --- Panel C: nucleotide composition, one colour tone per row ---
    ax = axes[2]
    rows = [("A", 0, "#2f7d4f"), ("C", 1, "#2b5fa8"), ("G", 2, "#d08a1e"),
            ("T", 3, "#b83b36"), ("gap", 5, "#5b5b5b")]
    for r, (label, col, color) in enumerate(rows):
        cmap = LinearSegmentedColormap.from_list(label, ["white", color])
        ax.imshow(np.nan_to_num(pct[:, col])[None, :], aspect="auto", origin="upper",
                  cmap=cmap, vmin=0, vmax=100, interpolation="nearest",
                  extent=(0.5, REF_LEN + 0.5, r + 0.5, r - 0.5))
    ax.set_ylim(len(rows) - 0.5, -0.5)
    ax.set_yticks(range(len(rows)))
    ax.set_yticklabels([label for label, _, _ in rows])
    ax.set_ylabel("Composition\nper position")
    for r in range(len(rows) - 1):
        ax.axhline(r + 0.5, color="white", lw=1.2)
    for _, a, b in V_REGIONS:
        ax.axvline(a, color="0.3", lw=0.5, alpha=0.6)
        ax.axvline(b, color="0.3", lw=0.5, alpha=0.6)
    ax.text(1.005, 0.5, "darker = higher\nshare (0→100%)", transform=ax.transAxes,
            fontsize=7.5, color="0.35", va="center")

    # --- Panel D: how many sequences cover each position ---
    ax = axes[3]
    shade(ax)
    ax.fill_between(x, 0, total, color="0.55", lw=0)
    ax.set_ylabel("Sequences\ncovering")
    ax.set_xlabel("Position along the 16S gene (E. coli numbering, bp)")
    ax.set_ylim(0, max(total) * 1.08)

    for a in axes:
        a.set_xlim(1, REF_LEN)
        a.spines[["top", "right"]].set_visible(False)
    axes[2].spines[["top", "right"]].set_visible(True)

    fig.savefig(png, dpi=300, bbox_inches="tight")
    fig.savefig(pdf, bbox_inches="tight")
    print("Wrote %s and %s" % (png, pdf), flush=True)


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--window", type=int, default=25, help="smoothing window in bp (default 25)")
    ap.add_argument("--min-identity", type=float, default=0.60,
                    help="minimum identity to the reference (default 0.60 - low enough to "
                         "keep Archaea, whose 16S is only ~60-75%% identical to E. coli)")
    ap.add_argument("--min-coverage", type=int, default=300,
                    help="minimum number of A/C/G/T bases a sequence must contribute "
                         "(default 300)")
    ap.add_argument("--limit", type=int, default=0,
                    help="process only the first N sequences (trial run)")
    args = ap.parse_args()

    ref = get_reference()
    print("Reference: %d bp | 515F primer lands at position %d"
          % (len(ref), ref.find("GTGCCAGCAGCCGCGGTAA") + 1), flush=True)

    counts, n_used, n_flipped, excluded, per_seq = build_matrix(ref, args)
    print("Used %d sequences (%d strand-flipped, %d excluded)"
          % (n_used, n_flipped, len(excluded)), flush=True)

    with open(OUT_IDENTITY, "w", encoding="utf-8") as fh:
        fh.write("accession\tpct_identity_vs_ecoli\tbases_used\tdescription\n")
        for acc, ident, cov, header in per_seq:
            fh.write("%s\t%.2f\t%d\t%s\n" % (acc, ident, cov, header))

    with open(OUT_EXCLUDED, "w", encoding="utf-8") as fh:
        fh.write("accession\treason\n")
        for acc, reason in excluded:
            fh.write("%s\t%s\n" % (acc, reason))

    total, pct, identity, entropy, entropy_std = compute_stats(counts)
    write_tsv(OUT_TSV, total, pct, identity, entropy, entropy_std)
    print("Wrote %s" % OUT_TSV, flush=True)

    mid = slice(400, 1100)
    print("Coverage: median %d sequences per position (core region %d)"
          % (int(np.median(total)), int(np.median(total[mid]))), flush=True)
    print("%% identity: mean %.1f | lowest %.1f (position %d) | highest %.1f (position %d)"
          % (np.nanmean(identity), np.nanmin(identity), int(np.nanargmin(identity)) + 1,
             np.nanmax(identity), int(np.nanargmax(identity)) + 1), flush=True)

    plot(total, pct, identity, entropy_std, n_used, args.window, OUT_PNG, OUT_PDF)


if __name__ == "__main__":
    main()
