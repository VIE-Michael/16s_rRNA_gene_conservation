#!/usr/bin/env python3
"""Export the 16S conservation results to an Excel workbook for per-position lookup.

Reads 16s_position_stats.tsv (produced by analyze_16s_conservation.py) and writes
16s_conservation_by_position.xlsx with two sheets:

  - "Per position": one row per position along the 16S gene (E. coli numbering),
    giving the most common nucleotide, the % conservation, the share of each
    nucleotide, the gap share, how many sequences cover the position, and which
    V region the position belongs to.
  - "Region summary": aggregated statistics for V1-V9, for the conserved stretches
    between them and for the whole gene, each with its consensus sequence.

Example:
    python3 export_conservation_excel.py
"""

import argparse
import os

HERE = os.path.dirname(os.path.abspath(__file__))
OUT = os.environ.get("SIXTEENS_OUT_DIR", HERE)
IN_TSV = os.path.join(OUT, "16s_position_stats.tsv")
OUT_XLSX = os.path.join(OUT, "16s_conservation_by_position.xlsx")

# Hypervariable regions in E. coli numbering (Chakravorty et al. 2007)
V_REGIONS = [("V1", 69, 99), ("V2", 137, 242), ("V3", 433, 497), ("V4", 576, 682),
             ("V5", 822, 879), ("V6", 986, 1043), ("V7", 1117, 1173),
             ("V8", 1243, 1294), ("V9", 1435, 1465)]

BASE_FILL = {"A": "D6EFDF", "C": "D6E2F5", "G": "FBECD2", "T": "F8DAD8", "-": "E8E8E8"}


def region_of(pos):
    for name, a, b in V_REGIONS:
        if a <= pos <= b:
            return name
    return ""


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--src", default=IN_TSV)
    ap.add_argument("--out", default=OUT_XLSX)
    ap.add_argument("--gap-dominant", type=float, default=50.0,
                    help="flag a position as indel-dominated when %% gap exceeds this")
    args = ap.parse_args()

    if not os.path.exists(args.src):
        raise SystemExit("%s not found - run analyze_16s_conservation.py first." % args.src)

    import openpyxl
    from openpyxl.formatting.rule import ColorScaleRule
    from openpyxl.styles import Alignment, Font, PatternFill
    from openpyxl.utils import get_column_letter

    rows = []
    with open(args.src, encoding="utf-8") as fh:
        header = fh.readline().rstrip("\n").split("\t")
        col = {name: i for i, name in enumerate(header)}
        for line in fh:
            p = line.rstrip("\n").split("\t")
            pos = int(p[col["position"]])
            pcts = {b: float(p[col["pct_" + b]]) for b in "ACGT"}
            consensus = max(pcts, key=pcts.get)
            gap = float(p[col["pct_gap"]])
            rows.append({
                "position": pos,
                "consensus": consensus,
                "conservation": float(p[col["pct_identity"]]),
                "entropy": float(p[col["entropy_std"]]),
                "A": pcts["A"], "C": pcts["C"], "G": pcts["G"], "T": pcts["T"],
                "gap": gap,
                "ambiguous": float(p[col["pct_ambiguous"]]),
                "n": int(p[col["n_covered"]]),
                "region": region_of(pos),
                "note": "indel-dominated" if gap >= args.gap_dominant else "",
            })

    wb = openpyxl.Workbook()

    # ------------------------------------------------------------------ sheet 1
    ws = wb.active
    ws.title = "Per position"
    headers = ["Position", "Nucleotide", "% conservation", "Standardized entropy",
               "% A", "% C", "% G", "% T", "% gap", "% ambiguous",
               "Sequences covering", "V region", "Note"]
    ws.append(headers)
    for r in rows:
        ws.append([r["position"], r["consensus"], r["conservation"], r["entropy"],
                   r["A"], r["C"], r["G"], r["T"], r["gap"], r["ambiguous"],
                   r["n"], r["region"], r["note"]])

    bold = Font(bold=True)
    for c in ws[1]:
        c.font = bold
        c.alignment = Alignment(horizontal="center", vertical="center", wrap_text=True)
    ws.freeze_panes = "C2"
    ws.auto_filter.ref = "A1:M%d" % (len(rows) + 1)

    for i, r in enumerate(rows, start=2):
        cell = ws.cell(row=i, column=2)
        cell.alignment = Alignment(horizontal="center")
        cell.font = Font(bold=True)
        cell.fill = PatternFill("solid", fgColor=BASE_FILL[r["consensus"]])
        for c in range(3, 11):
            ws.cell(row=i, column=c).number_format = "0.00"
        ws.cell(row=i, column=4).number_format = "0.000"

    # colour scale on % conservation: red (variable) -> yellow -> green (conserved)
    ws.conditional_formatting.add(
        "C2:C%d" % (len(rows) + 1),
        ColorScaleRule(start_type="num", start_value=25, start_color="F4776B",
                       mid_type="num", mid_value=70, mid_color="FFE08A",
                       end_type="num", end_value=100, end_color="63BE7B"))

    widths = [9, 12, 13, 13, 8, 8, 8, 8, 9, 12, 13, 10, 17]
    for i, w in enumerate(widths, start=1):
        ws.column_dimensions[get_column_letter(i)].width = w

    # ------------------------------------------------------------------ sheet 2
    ws2 = wb.create_sheet("Region summary")
    ws2.append(["Region", "From", "To", "Positions", "Mean % conservation",
                "Mean entropy", "Mean % gap", "Consensus sequence"])
    for c in ws2[1]:
        c.font = bold
        c.alignment = Alignment(horizontal="center", vertical="center", wrap_text=True)

    def block(name, a, b):
        sel = [r for r in rows if a <= r["position"] <= b]
        return [name, a, b, len(sel),
                sum(r["conservation"] for r in sel) / len(sel),
                sum(r["entropy"] for r in sel) / len(sel),
                sum(r["gap"] for r in sel) / len(sel),
                "".join(r["consensus"] for r in sel)]

    blocks = []
    prev_end = 0
    for name, a, b in V_REGIONS:
        if a > prev_end + 1:
            blocks.append(("conserved %d-%d" % (prev_end + 1, a - 1), prev_end + 1, a - 1))
        blocks.append((name, a, b))
        prev_end = b
    last = rows[-1]["position"]
    if prev_end < last:
        blocks.append(("conserved %d-%d" % (prev_end + 1, last), prev_end + 1, last))
    blocks.append(("WHOLE GENE", 1, last))

    for name, a, b in blocks:
        ws2.append(block(name, a, b))
    for i in range(2, ws2.max_row + 1):
        label = ws2.cell(row=i, column=1).value
        ws2.cell(row=i, column=1).font = Font(bold=not label.startswith("conserved"))
        for c in (5, 6, 7):
            ws2.cell(row=i, column=c).number_format = "0.00"
        ws2.cell(row=i, column=6).number_format = "0.000"
        ws2.cell(row=i, column=8).alignment = Alignment(vertical="top")
        ws2.cell(row=i, column=8).font = Font(name="Menlo", size=8)
    ws2.conditional_formatting.add(
        "E2:E%d" % ws2.max_row,
        ColorScaleRule(start_type="num", start_value=70, start_color="F4776B",
                       mid_type="num", mid_value=85, mid_color="FFE08A",
                       end_type="num", end_value=100, end_color="63BE7B"))
    for i, w in enumerate([18, 8, 8, 11, 14, 13, 12, 70], start=1):
        ws2.column_dimensions[get_column_letter(i)].width = w
    ws2.freeze_panes = "A2"

    wb.save(args.out)
    print("Wrote %s (%d positions, %d regions)" % (args.out, len(rows), ws2.max_row - 1))
    cons = "".join(r["consensus"] for r in rows)
    print("Consensus sequence is %d bp, starts: %s..." % (len(cons), cons[:40]))


if __name__ == "__main__":
    main()
