#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Figure 1 as a heatmap: antibiotic (rows) x species (columns).

    python ve_hinh1_heatmap.py

Why a heatmap instead of a forest plot:
  - The paper's main finding is that imipenem/meropenem and ertapenem/
    doripenem shift in OPPOSITE directions. In a heatmap those are adjacent
    rows, so a reader sees one blue band sitting directly above a red band.
  - It is far more compact: ~90 pairs fit on half a page.
  - Blank cells carry information too: where one system sets no breakpoint.

Outputs ket_qua/figure1.png and ket_qua/figure1.pdf
"""

from pathlib import Path
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.colors import TwoSlopeNorm, LinearSegmentedColormap
from matplotlib.patches import Rectangle

FILE_KQ = Path("./ket_qua/bang2_ket_qua.csv")
RA_PNG = Path("./ket_qua/figure1.png")
RA_PDF = Path("./ket_qua/figure1.pdf")

N_MIN = 100          # minimum isolates for a cell to be shown
HIEN_SO = True       # print the value inside each cell

# Column order: Enterobacterales first, then non-fermenters, then Gram-positive
THU_TU_LOAI = [
    "Escherichia coli", "Klebsiella pneumoniae", "Salmonella enterica",
    "Enterobacter cloacae", "Citrobacter freundii", "Serratia marcescens",
    "Pseudomonas aeruginosa", "Acinetobacter baumannii",
    "Staphylococcus aureus",
]
TEN_COT = {
    "Escherichia coli": "E. coli",
    "Klebsiella pneumoniae": "K. pneumoniae",
    "Salmonella enterica": "S. enterica",
    "Enterobacter cloacae": "E. cloacae",
    "Citrobacter freundii": "C. freundii",
    "Serratia marcescens": "S. marcescens",
    "Pseudomonas aeruginosa": "P. aeruginosa",
    "Acinetobacter baumannii": "A. baumannii",
    "Staphylococcus aureus": "S. aureus",
}

# Row order inside each class is fixed, so the carbapenem pattern lines up
NHOM_KS = [
    ("Penicillins", [
        "penicillin", "ampicillin", "amoxicillin/clavulanic acid",
        "ampicillin/sulbactam", "piperacillin", "piperacillin/tazobactam",
        "oxacillin"]),
    ("Cephalosporins", [
        "cefuroxime", "cefoxitin", "cefotaxime", "ceftriaxone", "ceftazidime",
        "cefepime", "cefazolin", "cefpodoxime", "ceftazidime/avibactam",
        "ceftolozane/tazobactam", "cefiderocol"]),
    ("Carbapenems", ["imipenem", "meropenem", "ertapenem", "doripenem"]),
    ("Monobactams", ["aztreonam"]),
    ("Aminoglycosides", ["amikacin", "gentamicin", "tobramycin"]),
    ("Fluoroquinolones", [
        "ciprofloxacin", "levofloxacin", "moxifloxacin", "nalidixic acid"]),
    ("Tetracyclines", ["tetracycline", "minocycline", "tigecycline"]),
    ("Other agents", [
        "trimethoprim/sulfamethoxazole", "chloramphenicol", "colistin",
        "nitrofurantoin", "fosfomycin", "vancomycin", "teicoplanin",
        "daptomycin", "linezolid", "erythromycin", "clindamycin",
        "rifampin", "fusidic acid", "quinupristin/dalfopristin",
        "mupirocin", "azithromycin"]),
]

# Diverging palette. Muted so it still reads when printed in greyscale.
CMAP = LinearSegmentedColormap.from_list("venh", [
    "#1a4a7a", "#4e82b4", "#a8c6de", "#f4f4f4",
    "#e8b4b4", "#c25c5c", "#8c2020",
])
MAU_TRONG = "#e9e9e9"     # no breakpoint in one of the two systems
MAU_IT_MAU = "#f7f7f7"    # below N_MIN


def nap():
    if not FILE_KQ.exists():
        raise SystemExit(f"{FILE_KQ} not found. Run phan_tich_nguong.py first.")
    d = pd.read_csv(FILE_KQ)
    d = d[d["n"] >= N_MIN].copy()
    if d.empty:
        raise SystemExit("No pair reaches the minimum sample size.")
    return d


def dung_ma_tran(d):
    hang = []
    nhan_nhom = []
    for ten_nhom, ds in NHOM_KS:
        co = [k for k in ds if k in set(d["khang_sinh"])]
        if co:
            nhan_nhom.append((ten_nhom, len(hang), len(hang) + len(co)))
            hang.extend(co)
    # anything not in the lists above
    con_lai = sorted(set(d["khang_sinh"]) - set(hang))
    if con_lai:
        nhan_nhom.append(("Unclassified", len(hang), len(hang) + len(con_lai)))
        hang.extend(con_lai)

    cot = [s for s in THU_TU_LOAI if s in set(d["loai"])]

    M = np.full((len(hang), len(cot)), np.nan)
    N = np.zeros_like(M)
    for r in d.itertuples():
        if r.khang_sinh in hang and r.loai in cot:
            M[hang.index(r.khang_sinh), cot.index(r.loai)] = r.chenh_diem_pt
            N[hang.index(r.khang_sinh), cot.index(r.loai)] = r.n
    return M, N, hang, cot, nhan_nhom


def ve(d):
    M, N, hang, cot, nhan_nhom = dung_ma_tran(d)
    nr, nc = M.shape

    gh = np.nanmax(np.abs(M))
    gh = float(np.ceil(gh / 10.0) * 10)
    norm = TwoSlopeNorm(vmin=-gh, vcenter=0.0, vmax=gh)

    fig, ax = plt.subplots(figsize=(1.02 * nc + 5.2, 0.34 * nr + 2.6))

    # grey backdrop so empty cells are visibly "no data", not white space
    ax.add_patch(Rectangle((-0.5, -0.5), nc, nr, facecolor=MAU_TRONG,
                           edgecolor="none", zorder=0))

    im = ax.imshow(M, cmap=CMAP, norm=norm, aspect="auto",
                   interpolation="nearest", zorder=1)

    # cell values
    if HIEN_SO:
        for i in range(nr):
            for j in range(nc):
                v = M[i, j]
                if np.isnan(v):
                    continue
                # pick text colour by how dark the cell is
                sang = abs(v) / gh < 0.55
                ax.text(j, i, f"{v:+.0f}" if abs(v) >= 0.5 else "0",
                        ha="center", va="center", fontsize=7.4,
                        color="#222222" if sang else "white", zorder=3)

    # thin white grid between cells
    ax.set_xticks(np.arange(-0.5, nc, 1), minor=True)
    ax.set_yticks(np.arange(-0.5, nr, 1), minor=True)
    ax.grid(which="minor", color="white", linewidth=1.1)
    ax.tick_params(which="minor", length=0)

    ax.set_xticks(range(nc))
    ax.set_xticklabels([TEN_COT.get(c, c) for c in cot], fontsize=8.8,
                       style="italic", rotation=42, ha="left")
    ax.xaxis.set_ticks_position("top")
    ax.set_yticks(range(nr))
    ax.set_yticklabels(hang, fontsize=8.2)
    ax.tick_params(axis="both", length=0)

    # separators and labels for antibiotic classes
    for ten, bd, kt in nhan_nhom:
        if bd > 0:
            ax.axhline(bd - 0.5, color="#333333", linewidth=1.3, zorder=4)
        giua = (bd + kt - 1) / 2
        ax.text(-0.40, giua, ten, transform=ax.get_yaxis_transform(),
                fontsize=8.2, fontweight="bold", va="center", ha="right",
                color="#222222")
        # thin bracket tying the label to its rows
        ax.plot([-0.375, -0.375], [bd - 0.38, kt - 1 + 0.38],
                transform=ax.get_yaxis_transform(), color="#999999",
                linewidth=0.9, clip_on=False, zorder=5)

    for s in ax.spines.values():
        s.set_visible(False)

    cb = fig.colorbar(im, ax=ax, fraction=0.028, pad=0.06,
                      ticks=np.linspace(-gh, gh, 5))
    cb.set_label("Difference in resistance rate\n(EUCAST minus CLSI, percentage points)",
                 fontsize=8.4, labelpad=10)
    cb.ax.tick_params(labelsize=7.8)
    cb.outline.set_visible(False)
    cb.ax.text(0.5, 1.05, "EUCAST higher", transform=cb.ax.transAxes,
               fontsize=7.4, ha="center", va="bottom", color="#8c2020")
    cb.ax.text(0.5, -0.06, "CLSI higher", transform=cb.ax.transAxes,
               fontsize=7.4, ha="center", va="top", color="#1a4a7a")

    ax.text(0.0, -0.035,
            f"Grey cells: no comparable breakpoint in one system, or fewer "
            f"than {N_MIN:,} isolates.",
            transform=ax.transAxes, fontsize=7.6, color="#666666",
            va="top", ha="left")

    fig.subplots_adjust(left=0.40)
    fig.savefig(RA_PNG, dpi=400, bbox_inches="tight", facecolor="white")
    fig.savefig(RA_PDF, bbox_inches="tight", facecolor="white")
    plt.close(fig)
    print(f"{nr} antibiotics x {nc} species -> {RA_PNG}")


if __name__ == "__main__":
    ve(nap())
