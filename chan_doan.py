#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Chan doan nguon goc cac gia tri MIC trong du lieu BV-BRC.

    python chan_doan.py

Voi moi file trong du_lieu/, dem so dong co MIC theo:
  - Evidence (Laboratory Method / Computational Method / trong)
  - Computational Method (vd "MIC XGBoost Model" = MIC DU DOAN, khong phai do that)
  - Laboratory Typing Method
  - co / khong co dau bat dang thuc
Ket qua ghi ra ket_qua/chan_doan.txt - gui file do la du.
"""

from pathlib import Path
from collections import Counter
import pandas as pd

THU_MUC = Path("./du_lieu")
RA = Path("./ket_qua/chan_doan.txt")
COT = ["Genome ID", "Antibiotic", "Measurement Sign", "Measurement Value",
       "Laboratory Typing Method", "Evidence", "Computational Method", "PubMed"]


def mot_file(f, ghi):
    dem = {k: Counter() for k in ("ev", "cm", "pp", "dau")}
    tong = co_mic = lab = co_pubmed_lab = 0
    chung_lab, chung_all = set(), set()
    for c in pd.read_csv(f, usecols=lambda x: x in COT, low_memory=False,
                         chunksize=500_000):
        tong += len(c)
        m = c[c["Measurement Value"].notna()]
        co_mic += len(m)
        chung_all.update(m["Genome ID"].unique())
        ev = m["Evidence"].fillna("(trong)") if "Evidence" in m else pd.Series("(khong co cot)", index=m.index)
        cm = m["Computational Method"].fillna("(trong)") if "Computational Method" in m else pd.Series("(khong co cot)", index=m.index)
        dem["ev"].update(ev)
        dem["cm"].update(cm)
        dem["pp"].update(m["Laboratory Typing Method"].fillna("(trong)"))
        la_lab = ev.str.lower().str.contains("laboratory") & (cm == "(trong)")
        ml = m[la_lab]
        lab += len(ml)
        chung_lab.update(ml["Genome ID"].unique())
        co_pubmed_lab += ml["PubMed"].notna().sum()
        dem["dau"].update(ml["Measurement Sign"].fillna("(trong)"))

    ghi(f"\n=== {f.name} ===")
    ghi(f"tong dong: {tong:,} | co MIC: {co_mic:,} | so chung co MIC: {len(chung_all):,}")
    if tong == 2_500_000:
        ghi("  !! DUNG 2.500.000 DONG -> file bi BV-BRC cat, tai chua du")
    pt = 100 * lab / co_mic if co_mic else 0
    ghi(f"MIC DO THAT (Evidence=Laboratory, khong co Computational Method): "
        f"{lab:,} dong ({pt:.1f}%), {len(chung_lab):,} chung")
    ghi(f"  trong do co ma PubMed: {co_pubmed_lab:,}")
    ghi("Evidence:            " + ", ".join(f"{k}={v:,}" for k, v in dem["ev"].most_common(5)))
    ghi("Computational Method: " + ", ".join(f"{k}={v:,}" for k, v in dem["cm"].most_common(5)))
    ghi("Phuong phap:         " + ", ".join(f"{k}={v:,}" for k, v in dem["pp"].most_common(6)))
    ghi("Dau (chi MIC do that): " + ", ".join(f"{k}={v:,}" for k, v in dem["dau"].most_common(6)))


def main():
    RA.parent.mkdir(exist_ok=True)
    dong = []

    def ghi(x):
        print(x)
        dong.append(x)

    for f in sorted(THU_MUC.glob("BVBRC_genome_amr_*.csv")):
        mot_file(f, ghi)
    RA.write_text("\n".join(dong), encoding="utf-8")
    print(f"\nDa ghi {RA}")


if __name__ == "__main__":
    main()
