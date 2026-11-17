#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Ghep nuoc phan lap + nam thu thap tu bang bo gen BV-BRC, roi phan tich
muc venh CLSI - EUCAST tach theo khu vuc.

    python phan_tich_khu_vuc.py

Can:
  du_lieu/                     8-9 file BVBRC_genome_amr_*.csv (nhu cu)
  du_lieu/bo_gen/              file tai tu tab GENOMES cua tung taxon
                               (ten file tuy y, mien la .csv)
  nguong.csv, ten_thuoc.csv    nhu cu
  phan_tich_nguong.py          cung thu muc (script nay dung lai ham cua no)

Dau ra trong ket_qua/:
  khu_vuc_mo_ta.csv       so chung va % co thong tin nuoc, theo loai x khu vuc
  bang3_theo_khu_vuc.csv  ket qua chinh, tinh rieng cho tung khu vuc
  so_sanh_WPR_EUR.csv     cac cap du mau o ca hai khu vuc, dat canh nhau
  figure2.png / .pdf      WPR (truc doc) so voi chau Au (truc ngang)
  nuoc_khong_xep_duoc.csv ten nuoc chua co trong bang quy doi - kiem tay
"""

from pathlib import Path
import re
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

import phan_tich_nguong as P

THU_MUC_BO_GEN = Path("./du_lieu/bo_gen")
RA = Path("./ket_qua")
N_MIN_KHU_VUC = 50     # nguong mau trong MOI khu vuc; thap hon toan bo vi chia nho

# Dai Loan khong co trong danh sach thanh vien WHO WPRO. Mac dinh xep vao
# "Tay Thai Binh Duong theo dia ly" va bao cao rieng trong bai. Doi thanh
# False neu muon dung dung danh sach thanh vien WHO.
TINH_DAI_LOAN_VAO_WPR = True

# --------------------------------------------------------------------------
# Bang quy doi nuoc -> khu vuc. Chu thuong, da bo dau va phan sau dau ':'.
# --------------------------------------------------------------------------
WPR = {
    "australia", "brunei", "brunei darussalam", "cambodia", "china",
    "cook islands", "fiji", "japan", "kiribati", "laos", "lao pdr",
    "lao people's democratic republic", "malaysia", "marshall islands",
    "micronesia", "federated states of micronesia", "mongolia", "nauru",
    "new zealand", "niue", "palau", "papua new guinea", "philippines",
    "south korea", "korea", "republic of korea", "korea, republic of",
    "samoa", "singapore", "solomon islands", "tonga", "tuvalu", "vanuatu",
    "viet nam", "vietnam", "american samoa", "french polynesia", "guam",
    "hong kong", "macao", "macau", "new caledonia",
    "northern mariana islands", "pitcairn", "tokelau", "wallis and futuna",
}
DAI_LOAN = {"taiwan"}
EUR = {
    "albania", "andorra", "armenia", "austria", "azerbaijan", "belarus",
    "belgium", "bosnia and herzegovina", "bulgaria", "croatia", "cyprus",
    "czech republic", "czechia", "denmark", "estonia", "finland", "france",
    "georgia", "germany", "greece", "hungary", "iceland", "ireland", "israel",
    "italy", "kazakhstan", "kosovo", "kyrgyzstan", "latvia", "lithuania",
    "luxembourg", "malta", "moldova", "monaco", "montenegro", "netherlands",
    "north macedonia", "macedonia", "norway", "poland", "portugal", "romania",
    "russia", "russian federation", "san marino", "serbia", "slovakia",
    "slovenia", "spain", "sweden", "switzerland", "tajikistan", "turkey",
    "turkiye", "turkmenistan", "ukraine", "united kingdom", "uk", "england",
    "scotland", "wales", "northern ireland", "uzbekistan",
    "england and wales", "the netherlands", "great britain",
}
AMR = {  # chau My
    "usa", "united states", "united states of america", "canada", "mexico",
    "brazil", "argentina", "chile", "colombia", "peru", "venezuela",
    "ecuador", "bolivia", "paraguay", "uruguay", "guatemala", "honduras",
    "nicaragua", "costa rica", "panama", "cuba", "dominican republic",
    "haiti", "jamaica", "trinidad and tobago", "puerto rico", "el salvador",
    "guyana", "suriname", "belize", "barbados", "bahamas",
    "turks and caicos islands", "cayman islands", "bermuda",
}
EMR = {  # Dong Dia Trung Hai cua WHO
    "afghanistan", "bahrain", "djibouti", "egypt", "iran", "iraq", "jordan",
    "kuwait", "lebanon", "libya", "morocco", "oman", "pakistan", "qatar",
    "saudi arabia", "somalia", "sudan", "syria", "tunisia",
    "united arab emirates", "uae", "yemen", "palestinian territory",
    "palestine", "west bank", "gaza",
}
AFR = {  # chau Phi cua WHO
    "algeria", "angola", "benin", "botswana", "burkina faso", "burundi",
    "cameroon", "cape verde", "central african republic", "chad", "comoros",
    "congo", "democratic republic of the congo", "cote d'ivoire",
    "ivory coast", "equatorial guinea", "eritrea", "eswatini", "swaziland",
    "ethiopia", "gabon", "gambia", "the gambia", "ghana", "guinea",
    "guinea-bissau", "kenya", "lesotho", "liberia", "madagascar", "malawi",
    "mali", "mauritania", "mauritius", "mozambique", "namibia", "niger",
    "nigeria", "rwanda", "senegal", "seychelles", "sierra leone",
    "south africa", "south sudan", "tanzania", "togo", "uganda", "zambia",
    "zimbabwe",
}
SEAR = {  # Dong Nam A cua WHO (khac WPR)
    "india", "bangladesh", "nepal", "sri lanka", "thailand", "myanmar",
    "indonesia", "bhutan", "maldives", "timor-leste", "east timor",
    "north korea", "democratic people's republic of korea",
}


TAT_CA = WPR | EUR | AMR | SEAR | EMR | AFR | DAI_LOAN


def chuan_nuoc(x):
    if pd.isna(x):
        return None
    s = str(x).strip().lower()
    if s not in TAT_CA:
        s = re.split(r"[:,(]", s)[0]      # "china: zhejiang", "united kingdom (england)"
    s = re.sub(r"\s+", " ", s).strip(" .")
    return s or None


def khu_vuc(nuoc):
    if nuoc is None or (not isinstance(nuoc, str) and pd.isna(nuoc)) or nuoc in {
            "no travel", "no travel reported", "no travel information",
            "not stated", "n", "y", "control", "europe", "asian continent",
            "european cruise - multiple destinations", "china/nigeria",
            "", "missing", "not collected", "unknown",
                               "not applicable", "na", "n/a", "-"}:
        return "Khong ro"
    if nuoc in WPR:
        return "WPR"
    if nuoc in DAI_LOAN:
        return "WPR" if TINH_DAI_LOAN_VAO_WPR else "Khac"
    if nuoc in EUR:
        return "EUR"
    if nuoc in AMR:
        return "AMR"
    if nuoc in SEAR:
        return "SEAR"
    if nuoc in EMR:
        return "EMR"
    if nuoc in AFR:
        return "AFR"
    return "Khac"


# --------------------------------------------------------------------------
# Doc bang bo gen
# --------------------------------------------------------------------------
def tim_cot(cols, *tu_khoa):
    for c in cols:
        cl = c.lower()
        if all(k in cl for k in tu_khoa):
            return c
    return None


def nap_bo_gen():
    files = sorted(THU_MUC_BO_GEN.glob("*.csv")) + sorted(THU_MUC_BO_GEN.glob("*.tsv"))
    if not files:
        raise SystemExit(f"Khong thay file nao trong {THU_MUC_BO_GEN}. "
                         "Tai tab GENOMES cua tung taxon bo vao day.")
    khung = []
    for f in files:
        sep = "\t" if f.suffix == ".tsv" else ","
        dau = pd.read_csv(f, sep=sep, nrows=0)
        c_id = tim_cot(dau.columns, "genome", "id")
        c_nuoc = tim_cot(dau.columns, "isolation", "country")
        c_nam = tim_cot(dau.columns, "collection", "year")
        c_ngay = tim_cot(dau.columns, "collection", "date")
        if c_id is None or c_nuoc is None:
            print(f"  [bo qua] {f.name}: khong thay cot Genome ID / Isolation Country")
            continue
        cot = [c for c in (c_id, c_nuoc, c_nam, c_ngay) if c]
        d = pd.read_csv(f, sep=sep, usecols=cot, dtype={c_id: str}, low_memory=False)
        d = d.rename(columns={c_id: "Genome ID", c_nuoc: "nuoc_goc"})
        nam = pd.to_numeric(d[c_nam], errors="coerce") if c_nam else pd.Series(np.nan, index=d.index)
        if c_ngay:
            nam2 = pd.to_numeric(d[c_ngay].astype(str).str.extract(r"((?:19|20)\d\d)")[0], errors="coerce")
            nam = nam.fillna(nam2)
        d["nam_thu_thap"] = nam
        khung.append(d[["Genome ID", "nuoc_goc", "nam_thu_thap"]])
        print(f"  {f.name}: {len(d):,} bo gen")
    g = pd.concat(khung, ignore_index=True).drop_duplicates("Genome ID")
    g["Genome ID"] = g["Genome ID"].str.strip()
    g["nuoc"] = g["nuoc_goc"].apply(chuan_nuoc)
    g["khu_vuc"] = g["nuoc"].apply(khu_vuc)
    return g


# --------------------------------------------------------------------------
# Hinh 2
# --------------------------------------------------------------------------
def ve_hinh2(ss, duong_dan):
    if ss.empty:
        print("  khong co cap nao du mau o ca hai khu vuc, bo qua hinh 2")
        return
    fig, ax = plt.subplots(figsize=(6.2, 6.0))
    gh = max(5, np.ceil(max(ss[["chenh_EUR", "chenh_WPR"]].abs().max()) / 5) * 5)
    ax.plot([-gh, gh], [-gh, gh], color="#999999", lw=0.9, ls="--", zorder=1)
    ax.axhline(0, color="#cccccc", lw=0.8, zorder=0)
    ax.axvline(0, color="#cccccc", lw=0.8, zorder=0)
    mau = np.where(ss["khac_ro"], "#8c2d2d", "#2d4a6e")
    ax.scatter(ss["chenh_EUR"], ss["chenh_WPR"], c=mau, s=30, zorder=4)
    ax.errorbar(ss["chenh_EUR"], ss["chenh_WPR"],
                xerr=[ss["chenh_EUR"] - ss["eur_duoi"], ss["eur_tren"] - ss["chenh_EUR"]],
                yerr=[ss["chenh_WPR"] - ss["wpr_duoi"], ss["wpr_tren"] - ss["chenh_WPR"]],
                fmt="o", ms=5, color="#2d4a6e", ecolor="#9fb3c8", elinewidth=0.9,
                capsize=0, zorder=3)
    # chi ghi nhan toi da 6 cap khac biet ro nhat, tranh chong chu
    nhan = ss[ss["khac_ro"]].reindex(
        ss[ss["khac_ro"]]["hieu_WPR_tru_EUR"].abs().sort_values(ascending=False).index).head(6)
    for r in nhan.itertuples():
        ax.annotate(f"{r.loai.split()[0][0]}. {r.loai.split()[1]} – {r.khang_sinh}",
                    (r.chenh_EUR, r.chenh_WPR), fontsize=7, xytext=(4, 3),
                    textcoords="offset points", color="#333333", style="italic")
    ax.set_xlim(-gh, gh)
    ax.set_ylim(-gh, gh)
    ax.set_xlabel("European Region: EUCAST minus CLSI\n(percentage points)", fontsize=9)
    ax.set_ylabel("Western Pacific Region: EUCAST minus CLSI\n(percentage points)", fontsize=9)
    ax.tick_params(labelsize=8)
    for s in ("top", "right"):
        ax.spines[s].set_visible(False)
    ax.text(0.02, 0.98, "Dashed line: equal difference in both regions\n"
            "Red: confidence intervals do not overlap",
            transform=ax.transAxes, fontsize=7.4, color="#666666", va="top")
    fig.tight_layout()
    fig.savefig(duong_dan.with_suffix(".png"), dpi=400, facecolor="white")
    fig.savefig(duong_dan.with_suffix(".pdf"), facecolor="white")
    plt.close(fig)


# --------------------------------------------------------------------------
def main():
    RA.mkdir(exist_ok=True)

    print("\n[1/5] Doc bang bo gen")
    g = nap_bo_gen()

    print("\n[2/5] Doc du lieu MIC (dung lai phan_tich_nguong.py)")
    d, _ = P.doc_tat_ca()
    d["Genome ID"] = d["Genome ID"].astype(str).str.strip()

    print("\n[3/5] Ghep nuoc phan lap")
    d = d.merge(g[["Genome ID", "nuoc", "khu_vuc", "nam_thu_thap"]],
                on="Genome ID", how="left")
    d["khu_vuc"] = d["khu_vuc"].fillna("Khong ro")
    ghep = d["nuoc"].notna().mean() * 100
    print(f"  {ghep:.1f}% ket qua MIC co ten nuoc")

    khong_xep = (d[d["khu_vuc"] == "Khac"]["nuoc"].value_counts()
                 .rename_axis("nuoc").reset_index(name="so_ket_qua"))
    khong_xep.to_csv(RA / "nuoc_khong_xep_duoc.csv", index=False)
    if len(khong_xep):
        print(f"  {len(khong_xep)} ten nuoc xep vao 'Khac' - kiem file nuoc_khong_xep_duoc.csv")

    print("\n[4/5] Ap nguong")
    bp = P.nap_nguong()
    d = P.ap_hai_he(d, bp)

    # Mo ta theo loai x khu vuc (dem chung, khong dem ket qua)
    mt = (d.groupby(["loai", "khu_vuc"])["Genome ID"].nunique()
          .unstack(fill_value=0))
    mt["tong"] = mt.sum(axis=1)
    for k in ("WPR", "EUR", "AMR", "SEAR", "EMR", "AFR", "Khac", "Khong ro"):
        if k not in mt:
            mt[k] = 0
    mt["pct_co_nuoc"] = 100 * (1 - mt["Khong ro"] / mt["tong"])
    mt["pct_WPR"] = 100 * mt["WPR"] / mt["tong"]
    nam = d.groupby("loai")["nam_thu_thap"].agg(["min", "max"])
    mt = mt.join(nam.rename(columns={"min": "nam_thu_som", "max": "nam_thu_muon"}))
    mt = mt[["tong", "WPR", "EUR", "AMR", "SEAR", "EMR", "AFR", "Khac", "Khong ro",
             "pct_co_nuoc", "pct_WPR", "nam_thu_som", "nam_thu_muon"]]
    mt.reset_index().to_csv(RA / "khu_vuc_mo_ta.csv", index=False)

    # Nuoc dong gop nhieu nhat trong tung khu vuc (dem chung)
    tn = (d.dropna(subset=["nuoc"]).groupby(["khu_vuc", "nuoc"])["Genome ID"]
          .nunique().rename("so_chung").reset_index()
          .sort_values(["khu_vuc", "so_chung"], ascending=[True, False]))
    tn.to_csv(RA / "nuoc_theo_khu_vuc.csv", index=False)
    print("  5 nuoc nhieu nhat o WPR:",
          ", ".join(f"{r.nuoc} {r.so_chung}" for r in tn[tn.khu_vuc == "WPR"].head(5).itertuples()))
    thieu = mt.index[mt["pct_co_nuoc"] < 5].tolist()
    if thieu:
        print(f"  !! Loai gan nhu khong ghep duoc nuoc: {thieu} - kiem file bo gen")
    print(mt.round(1).to_string())

    print("\n[5/5] Ket qua theo khu vuc")
    giu = P.N_TOI_THIEU
    P.N_TOI_THIEU = N_MIN_KHU_VUC
    kq = []
    for kv in ("WPR", "EUR", "AMR", "SEAR", "EMR", "AFR"):
        sub = d[d["khu_vuc"] == kv]
        if sub.empty:
            continue
        r, _ = P.tinh_ket_qua(sub)
        if len(r):
            r.insert(0, "khu_vuc", kv)
            kq.append(r)
    P.N_TOI_THIEU = giu
    if not kq:
        raise SystemExit("Khong khu vuc nao co cap du mau.")
    kq = pd.concat(kq, ignore_index=True)
    kq.to_csv(RA / "bang3_theo_khu_vuc.csv", index=False)
    print(f"  {len(kq)} dong (cap x khu vuc), N_MIN = {N_MIN_KHU_VUC} moi khu vuc")

    w = kq[kq.khu_vuc == "WPR"].set_index(["loai", "khang_sinh"])
    e = kq[kq.khu_vuc == "EUR"].set_index(["loai", "khang_sinh"])
    chung = w.index.intersection(e.index)
    ss = pd.DataFrame({
        "n_WPR": w.loc[chung, "n"], "n_EUR": e.loc[chung, "n"],
        "khang_clsi_WPR": w.loc[chung, "khang_clsi_pct"],
        "khang_clsi_EUR": e.loc[chung, "khang_clsi_pct"],
        "chenh_WPR": w.loc[chung, "chenh_diem_pt"],
        "wpr_duoi": w.loc[chung, "chenh_ktc_duoi"], "wpr_tren": w.loc[chung, "chenh_ktc_tren"],
        "chenh_EUR": e.loc[chung, "chenh_diem_pt"],
        "eur_duoi": e.loc[chung, "chenh_ktc_duoi"], "eur_tren": e.loc[chung, "chenh_ktc_tren"],
        "do_khop_WPR": w.loc[chung, "do_khop_pct"], "do_khop_EUR": e.loc[chung, "do_khop_pct"],
    }).reset_index()
    ss["hieu_WPR_tru_EUR"] = ss["chenh_WPR"] - ss["chenh_EUR"]
    # hai KTC khong chong nhau -> khac biet ro giua hai khu vuc
    ss["khac_ro"] = (ss["wpr_duoi"] > ss["eur_tren"]) | (ss["wpr_tren"] < ss["eur_duoi"])
    ss = ss.sort_values("hieu_WPR_tru_EUR", key=abs, ascending=False)
    ss.to_csv(RA / "so_sanh_WPR_EUR.csv", index=False)
    print(f"  {len(ss)} cap du mau o ca WPR va chau Au; "
          f"{int(ss['khac_ro'].sum())} cap co muc venh khac ro giua hai khu vuc")
    if len(ss):
        print(ss.head(10)[["loai", "khang_sinh", "n_WPR", "n_EUR", "chenh_WPR",
                           "chenh_EUR", "khac_ro"]].round(1).to_string(index=False))

    ve_hinh2(ss, RA / "figure2")
    print(f"\nXong. Xem {RA.resolve()}")


if __name__ == "__main__":
    main()
