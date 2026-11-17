#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
So sanh nguong CLSI va EUCAST tren du lieu MIC cong khai cua BV-BRC.

Cach dung:
    python phan_tich_nguong.py

Truoc khi chay, sua phan CAU HINH ben duoi cho dung duong dan may ban.

Dau ra (trong thu muc ket_qua/):
    bang1_mo_ta.csv       - mo ta du lieu theo loai
    bang2_ket_qua.csv     - ket qua chinh theo cap loai-khang sinh
    hinh1_forest.png      - hinh chinh cua bai
    nhat_ky_loc.csv       - so dong bi loai o moi buoc (phu luc 2)
"""

import os
import math
import warnings
from pathlib import Path

import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

warnings.filterwarnings("ignore", category=FutureWarning)


# ============================================================
# CAU HINH - sua phan nay cho dung may ban
# ============================================================

THU_MUC_DU_LIEU = Path("./du_lieu")      # noi chua cac file BVBRC_genome_amr_*.csv
FILE_NGUONG = Path("./nguong.csv")       # file nguong ban tu lap
FILE_TEN_THUOC = Path("./ten_thuoc.csv") # bang quy doi ten khang sinh
THU_MUC_RA = Path("./ket_qua")

# Taxon ID -> ten loai. Bo bot dong neu chua tai du.
CAC_LOAI = {
    573:  "Klebsiella pneumoniae",
    590:  "Salmonella enterica",
    1280: "Staphylococcus aureus",
    287:  "Pseudomonas aeruginosa",
    550:  "Enterobacter cloacae",
    546:  "Citrobacter freundii",
    615:  "Serratia marcescens",
    562:  "Escherichia coli",
    470:  "Acinetobacter baumannii",
}

# Nhom loai dung chung mot bang nguong.
# Trong nguong.csv, cot 'loai' dien ten nhom -> ap cho moi loai trong nhom.
# Dien ten mot loai cu the -> chi ap cho loai do, va THANG ngưỡng cua nhom.
NHOM = {
    "Enterobacterales": [
        "Escherichia coli", "Klebsiella pneumoniae", "Salmonella enterica",
        "Enterobacter cloacae", "Citrobacter freundii", "Serratia marcescens",
    ],
}

# Nguong co mau toi thieu de dua vao bang ket qua
N_TOI_THIEU = 100

# Doc file lon theo khoi de khong het RAM. Dat None de doc mot lan.
KICH_THUOC_KHOI = 500_000


# ============================================================
# BUOC 1 - DOC VA LOC DU LIEU
# ============================================================

COT_CAN = [
    "Genome ID", "Genome Name", "Antibiotic",
    "Resistant Phenotype", "Measurement Sign", "Measurement Value",
    "Measurement Unit", "Laboratory Typing Method",
    "Testing Standard", "Testing Standard Year", "PubMed",
    "Evidence", "Computational Method",
]

# BV-BRC luu ca MIC DU DOAN bang may hoc (Evidence = "Computational Method",
# Computational Method = "MIC XGBoost Model") chung voi MIC do that.
# Nghien cuu nay chi duoc dung MIC do trong phong xet nghiem.
CHI_GIU_MIC_DO_THAT = True


def chuan_hoa_ten(s):
    """Dua ten khang sinh ve dang so sanh duoc: chu thuong, bo khoang trang thua."""
    if pd.isna(s):
        return s
    return str(s).strip().lower()


def nap_ten_thuoc():
    """
    Doc bang quy doi ten khang sinh.

    Cot 'do_tin_cay' quyet dinh cach xu ly:
      chac            -> gop ten_goc thanh ten_chuan
      LOAI            -> bo han khoi phan tich
      CAN_QUYET_DINH  -> neu ten_chuan de trong thi bo, va IN CANH BAO
                         de ban biet con viec phai quyet dinh
    """
    if not FILE_TEN_THUOC.exists():
        print(f"  [canh bao] khong thay {FILE_TEN_THUOC}, bo qua buoc gop ten.")
        print("             Du lieu co nhieu ten trung nhau, ket qua se bi chia nho sai.")
        return {}, set(), []

    t = pd.read_csv(FILE_TEN_THUOC)
    t["ten_goc"] = t["ten_goc"].apply(chuan_hoa_ten)
    t["ten_chuan"] = t["ten_chuan"].apply(
        lambda x: chuan_hoa_ten(x) if pd.notna(x) and str(x).strip() else None
    )
    t["do_tin_cay"] = t["do_tin_cay"].fillna("").astype(str).str.strip()

    doi_ten, bo_han, chua_quyet = {}, set(), []

    for _, r in t.iterrows():
        goc, chuan, muc = r["ten_goc"], r["ten_chuan"], r["do_tin_cay"]
        if muc == "LOAI":
            bo_han.add(goc)
        elif muc == "CAN_QUYET_DINH":
            if chuan:
                doi_ten[goc] = chuan
            else:
                bo_han.add(goc)
                chua_quyet.append((goc, r.get("ghi_chu", "")))
        elif chuan:
            doi_ten[goc] = chuan

    return doi_ten, bo_han, chua_quyet


def gop_ten(d, doi_ten, bo_han):
    """Ap bang quy doi. Phai chay TRUOC drop_duplicates."""
    d["khang_sinh"] = d["khang_sinh"].replace(doi_ten)
    d = d[~d["khang_sinh"].isin(bo_han)]
    return d


def doc_mot_loai(duong_dan, taxon_id, ten_loai, nhat_ky, doi_ten=None, bo_han=None):
    """Doc mot file BV-BRC, loc, tra ve DataFrame sach."""
    khung = []
    doc = pd.read_csv(
        duong_dan, usecols=lambda c: c in COT_CAN,
        dtype={"Genome ID": str},
        low_memory=False, chunksize=KICH_THUOC_KHOI,
    )
    for khoi in doc:
        khung.append(khoi)
    d = pd.concat(khung, ignore_index=True)

    n0 = len(d)

    # Chi giu dong co gia tri MIC
    d = d[d["Measurement Value"].notna()].copy()
    n1 = len(d)

    # Loai MIC du doan bang may hoc. Chi giu ban ghi co bang chung phong xet nghiem.
    n_truoc_ev = len(d)
    if CHI_GIU_MIC_DO_THAT:
        ev = d["Evidence"].astype(str).str.lower() if "Evidence" in d else pd.Series("", index=d.index)
        cm = d["Computational Method"] if "Computational Method" in d else pd.Series(pd.NA, index=d.index)
        d = d[ev.str.contains("laboratory") & cm.isna()]
    n_bo_du_doan = n_truoc_ev - len(d)

    # Don vi phai la mg/L (hoac ug/mL, tuong duong)
    don_vi = d["Measurement Unit"].astype(str).str.lower()
    d = d[don_vi.str.contains("mg/l") | don_vi.str.contains("ug/ml")
          | don_vi.str.contains("µg/ml")]
    n2 = len(d)

    # Gia tri phai doi duoc sang so
    d["mic"] = pd.to_numeric(d["Measurement Value"], errors="coerce")
    d = d[d["mic"].notna() & (d["mic"] > 0)]
    n3 = len(d)

    # Chi giu phuong phap dinh luong. Loai khuech tan dia.
    pp = d["Laboratory Typing Method"].astype(str).str.lower()
    d = d[~pp.str.contains("disk", na=False)]
    n4 = len(d)

    # Chuan hoa
    d["khang_sinh"] = d["Antibiotic"].apply(chuan_hoa_ten)
    d = d[d["khang_sinh"].notna()]
    d["dau"] = d["Measurement Sign"].fillna("=").astype(str).str.strip()
    d["loai"] = ten_loai
    d["taxon_id"] = taxon_id

    # Gop ten khang sinh. PHAI chay truoc drop_duplicates, neu khong
    # mot chung co ca "geamycin" va "gentamicin" se bi dem thanh hai.
    truoc_gop = len(d)
    if doi_ten or bo_han:
        d = gop_ten(d, doi_ten or {}, bo_han or set())
    n_bo_theo_ten = truoc_gop - len(d)
    n5 = len(d)

    # Mot chung co the co nhieu ket qua cho cung mot khang sinh.
    # Giu ban ghi dau tien, ghi lai so luong trung.
    truoc_khi_bo = len(d)
    d = d.drop_duplicates(subset=["Genome ID", "khang_sinh"], keep="first")
    n6 = len(d)

    nhat_ky.append({
        "loai": ten_loai,
        "tong_dong_goc": n0,
        "co_MIC": n1,
        "bo_MIC_du_doan_hoac_khong_ro_nguon": n_bo_du_doan,
        "dung_don_vi": n2,
        "MIC_la_so": n3,
        "loai_khuech_tan_dia": n4,
        "sau_khi_gop_ten": n5,
        "bo_do_ten_thuoc": n_bo_theo_ten,
        "sau_khi_bo_trung": n6,
        "so_ban_ghi_trung_da_bo": truoc_khi_bo - n6,
    })
    return d


def doc_tat_ca():
    doi_ten, bo_han, chua_quyet = nap_ten_thuoc()
    if doi_ten:
        print(f"  bang ten: gop {len(doi_ten)} ten, loai {len(bo_han)} ten")
    if chua_quyet:
        print("\n  !! CON MUC CHUA QUYET DINH trong ten_thuoc.csv, "
              "dang tam loai khoi phan tich:")
        for goc, ghi_chu in chua_quyet:
            print(f"     - {goc}")
            if ghi_chu and pd.notna(ghi_chu):
                print(f"       {str(ghi_chu)[:110]}")
        print()

    nhat_ky = []
    cac_khung = []
    for taxon_id, ten in CAC_LOAI.items():
        dd = THU_MUC_DU_LIEU / f"BVBRC_genome_amr_{taxon_id}.csv"
        if not dd.exists():
            print(f"  [bo qua] khong thay {dd.name}")
            continue
        print(f"  doc {ten} ...", end=" ", flush=True)
        d = doc_mot_loai(dd, taxon_id, ten, nhat_ky, doi_ten, bo_han)
        print(f"{len(d):,} dong dung duoc")
        cac_khung.append(d)
    if not cac_khung:
        raise SystemExit("Khong doc duoc file nao. Kiem tra THU_MUC_DU_LIEU.")
    return pd.concat(cac_khung, ignore_index=True), pd.DataFrame(nhat_ky)


# ============================================================
# BUOC 2 - AP NGUONG
# ============================================================

def phan_loai(mic, dau, s_bp, r_bp):
    """
    Tra ve 'S', 'I', 'R', hoac None neu khong ket luan duoc chac chan.

    Quy uoc: S neu MIC <= s_bp; R neu MIC > r_bp; con lai la I.

    Xu ly dau bat dang thuc - day la phan quan trong nhat ve mat khoa hoc:
      "<=x" : MIC that nam dau do duoi x. Chi ket luan S khi x <= s_bp.
      "<x"  : tuong tu, chi ket luan S khi x <= s_bp.
      ">x"  : MIC that lon hon x. Chi ket luan R khi x >= r_bp.
      ">=x" : MIC that tu x tro len. Chi ket luan R khi x > r_bp.
      "="   : biet chinh xac, phan loai binh thuong.

    Khong ket luan duoc thi tra None, KHONG doan.
    """
    if pd.isna(s_bp) or pd.isna(r_bp) or pd.isna(mic):
        return None

    # Nguong nhay 0.001 = he do KHONG co hang nhay (EUCAST "susceptible,
    # increased exposure"; CLSI colistin chi co I va R). Khi do moi MIC
    # khong vuot nguong khang deu la I, ke ca gia tri bi cat "<=x".
    khong_co_S = s_bp <= 0.0011

    if dau in ("<=", "≤", "<"):
        if mic <= s_bp:
            return "S"
        if khong_co_S and mic <= r_bp:
            return "I"
        return None
    if dau == ">":
        return "R" if mic >= r_bp else None
    if dau in (">=", "≥"):
        return "R" if mic > r_bp else None

    # dau "=" hoac trong
    if mic <= s_bp:
        return "S"
    if mic > r_bp:
        return "R"
    return "I"


def nap_nguong():
    """
    Doc file nguong.csv. Cau truc bat buoc:
        loai, khang_sinh, he, phien_ban, nguong_nhay, nguong_khang
    - 'he' chi nhan hai gia tri: CLSI hoac EUCAST
    - 'loai' de trong nghia la ap cho moi loai
    """
    if not FILE_NGUONG.exists():
        raise SystemExit(
            f"Khong thay {FILE_NGUONG}. Xem file nguong_mau.csv de biet cau truc."
        )
    bp = pd.read_csv(FILE_NGUONG)
    can = {"loai", "khang_sinh", "he", "nguong_nhay", "nguong_khang"}
    thieu = can - set(bp.columns)
    if thieu:
        raise SystemExit(f"nguong.csv thieu cot: {thieu}")

    bp["khang_sinh"] = bp["khang_sinh"].apply(chuan_hoa_ten)
    bp["loai"] = bp["loai"].fillna("").astype(str).str.strip()
    bp["he"] = bp["he"].astype(str).str.strip().str.upper()

    # Bung ten nhom thanh tung loai. Dong ghi ten nhom tro thanh nhieu dong.
    mo_rong = []
    for _, r in bp.iterrows():
        if r["loai"] in NHOM:
            for loai in NHOM[r["loai"]]:
                moi = r.copy()
                moi["loai"] = loai
                moi["_tu_nhom"] = True
                mo_rong.append(moi)
        else:
            r = r.copy()
            r["_tu_nhom"] = False
            mo_rong.append(r)
    bp = pd.DataFrame(mo_rong).reset_index(drop=True)

    # Dong ghi ro ten loai THANG dong bung ra tu nhom.
    truoc = len(bp)
    bp = bp.sort_values("_tu_nhom")  # False truoc -> keep='first' giu dong rieng
    bp = bp.drop_duplicates(subset=["loai", "khang_sinh", "he"], keep="first")
    n_bi_thay = truoc - len(bp)
    if n_bi_thay:
        print(f"  {n_bi_thay} nguong nhom bi nguong rieng cua loai ghi de")

    return bp.drop(columns=["_tu_nhom"])


def ghep_nguong(d, bp, he):
    """Ghep nguong cua mot he vao du lieu. Uu tien nguong rieng cho loai."""
    b = bp[bp["he"] == he]

    rieng = b[b["loai"] != ""][["loai", "khang_sinh", "nguong_nhay", "nguong_khang"]]
    chung = b[b["loai"] == ""][["khang_sinh", "nguong_nhay", "nguong_khang"]]

    d = d.merge(rieng, on=["loai", "khang_sinh"], how="left")
    d = d.merge(chung, on="khang_sinh", how="left", suffixes=("", "_chung"))

    # Nguong rieng thang nguong chung
    d["s_bp"] = d["nguong_nhay"].fillna(d.get("nguong_nhay_chung"))
    d["r_bp"] = d["nguong_khang"].fillna(d.get("nguong_khang_chung"))

    bo = ["nguong_nhay", "nguong_khang", "nguong_nhay_chung", "nguong_khang_chung"]
    d = d.drop(columns=[c for c in bo if c in d.columns])
    return d


def ap_hai_he(d, bp):
    print("  ap nguong CLSI ...", end=" ", flush=True)
    t = ghep_nguong(d.copy(), bp, "CLSI")
    d["co_bp_clsi"] = t["s_bp"].notna().values & t["r_bp"].notna().values
    d["nhan_clsi"] = [
        phan_loai(m, dau, s, r)
        for m, dau, s, r in zip(t["mic"], t["dau"], t["s_bp"], t["r_bp"])
    ]
    print("xong")

    print("  ap nguong EUCAST ...", end=" ", flush=True)
    t = ghep_nguong(d.copy(), bp, "EUCAST")
    d["co_bp_eucast"] = t["s_bp"].notna().values & t["r_bp"].notna().values
    d["nhan_eucast"] = [
        phan_loai(m, dau, s, r)
        for m, dau, s, r in zip(t["mic"], t["dau"], t["s_bp"], t["r_bp"])
    ]
    print("xong")
    return d


# ============================================================
# BUOC 3 - THONG KE
# ============================================================

def wilson(k, n, z=1.96):
    """Khoang tin cay Wilson cho mot ty le. Tra ve (duoi, tren)."""
    if n == 0:
        return (np.nan, np.nan)
    p = k / n
    giua = p + z * z / (2 * n)
    nua = z * math.sqrt((p * (1 - p) + z * z / (4 * n)) / n)
    mau = 1 + z * z / n
    return ((giua - nua) / mau, (giua + nua) / mau)


def ktc_chenh_cap(b, c, n, z=1.96):
    """
    KTC cho chenh lech hai ty le do tren CUNG mot chung (du lieu cap doi).
    b = so chung CLSI goi R ma EUCAST khong goi R
    c = so chung EUCAST goi R ma CLSI khong goi R
    Dung sai so chuan cua hieu hai ty le cap doi.
    """
    if n == 0:
        return (np.nan, np.nan)
    d = (c - b) / n   # eucast - clsi
    var = (b + c - (c - b) ** 2 / n) / (n * n)
    if var < 0:
        var = 0.0
    nua = z * math.sqrt(var)
    return (d - nua, d + nua)


def tinh_ket_qua(d):
    """Tinh ba chi so chinh cho moi cap loai - khang sinh."""
    ok = d[d["nhan_clsi"].notna() & d["nhan_eucast"].notna()].copy()

    ok["doi"] = ok["nhan_clsi"] != ok["nhan_eucast"]
    ok["doi_lon"] = (
        ((ok["nhan_clsi"] == "S") & (ok["nhan_eucast"] == "R"))
        | ((ok["nhan_clsi"] == "R") & (ok["nhan_eucast"] == "S"))
    )
    ok["R_clsi"] = ok["nhan_clsi"] == "R"
    ok["R_eucast"] = ok["nhan_eucast"] == "R"

    dong = []
    for (loai, ks), g in ok.groupby(["loai", "khang_sinh"]):
        n = len(g)
        if n < N_TOI_THIEU:
            continue

        k_doi = int(g["doi"].sum())
        k_lon = int(g["doi_lon"].sum())
        k_rc = int(g["R_clsi"].sum())
        k_re = int(g["R_eucast"].sum())

        # o lech cho du lieu cap doi
        b = int((g["R_clsi"] & ~g["R_eucast"]).sum())
        c = int((~g["R_clsi"] & g["R_eucast"]).sum())

        d_duoi, d_tren = ktc_chenh_cap(b, c, n)
        doi_duoi, doi_tren = wilson(k_doi, n)

        # Tach S / I / R cho tung he, theo khuon bang cua Cusack 2019 CMI
        sir = {}
        for he, cot in (("clsi", "nhan_clsi"), ("eucast", "nhan_eucast")):
            for hang in ("S", "I", "R"):
                sir[f"{he}_{hang}_pct"] = 100 * (g[cot] == hang).sum() / n

        dong.append({
            "loai": loai,
            "khang_sinh": ks,
            "n": n,
            # Do khop phan loai: chi so chinh trong tai lieu linh vuc nay.
            # Bang 100 - ty le doi. Bao cao chi so nay de so duoc voi cac bai truoc.
            "do_khop_pct": 100 * (n - k_doi) / n,
            **sir,
            "khang_clsi_pct": 100 * k_rc / n,
            "khang_eucast_pct": 100 * k_re / n,
            "chenh_diem_pt": 100 * (k_re - k_rc) / n,
            "chenh_ktc_duoi": 100 * d_duoi,
            "chenh_ktc_tren": 100 * d_tren,
            "ty_le_doi_pct": 100 * k_doi / n,
            "doi_ktc_duoi": 100 * doi_duoi,
            "doi_ktc_tren": 100 * doi_tren,
            "ty_le_doi_lon_pct": 100 * k_lon / n,
            "clsi_R_eucast_khongR": b,
            "eucast_R_clsi_khongR": c,
        })

    kq = pd.DataFrame(dong)
    if len(kq):
        kq = kq.sort_values("do_khop_pct").reset_index(drop=True)
    return kq, ok


def bang_mo_ta(d, ok):
    """Bang 1: mo ta du lieu theo loai."""
    dong = []
    for loai, g in d.groupby("loai"):
        g_ok = ok[ok["loai"] == loai]
        nam = pd.to_numeric(g["Testing Standard Year"], errors="coerce").dropna()
        dong.append({
            "loai": loai,
            "so_chung": g["Genome ID"].nunique(),
            "so_cap_co_MIC": len(g),
            "so_khang_sinh": g["khang_sinh"].nunique(),
            "so_cap_phan_loai_duoc_ca_hai_he": len(g_ok),
            "ty_le_phan_loai_duoc_pct": 100 * len(g_ok) / len(g) if len(g) else 0,
            "so_nghien_cuu_nguon": g["PubMed"].nunique(dropna=True),
            "nam_som_nhat": int(nam.min()) if len(nam) else None,
            "nam_muon_nhat": int(nam.max()) if len(nam) else None,
        })
    return pd.DataFrame(dong).sort_values("so_cap_co_MIC", ascending=False)


# ============================================================
# BUOC 4 - HINH 1
# ============================================================

def ve_hinh1(kq, duong_dan, so_dong=30):
    """
    Forest plot: chenh lech ty le khang (EUCAST - CLSI) theo diem phan tram.
    Moi cap loai-khang sinh mot diem, co thanh khoang tin cay, duong doc o 0.
    """
    if kq.empty:
        print("  khong co du lieu de ve")
        return

    v = kq.reindex(kq["chenh_diem_pt"].abs().sort_values(ascending=False).index)
    v = v.head(so_dong).sort_values("chenh_diem_pt")

    nhan = [f"{r.loai.split()[0][0]}. {r.loai.split()[1]} — {r.khang_sinh}  (n={r.n:,})"
            for r in v.itertuples()]
    y = np.arange(len(v))

    cao = max(4.0, 0.32 * len(v) + 1.5)
    fig, ax = plt.subplots(figsize=(9.5, cao))

    mau = ["#b2182b" if x > 0 else "#2166ac" for x in v["chenh_diem_pt"]]
    ax.hlines(y, v["chenh_ktc_duoi"], v["chenh_ktc_tren"],
              color=mau, linewidth=2.0, alpha=0.75)
    ax.scatter(v["chenh_diem_pt"], y, color=mau, s=34, zorder=3)
    ax.axvline(0, color="#444444", linewidth=1.0, linestyle="-")

    ax.set_yticks(y)
    ax.set_yticklabels(nhan, fontsize=8.5)
    ax.set_xlabel("Chênh lệch tỷ lệ kháng, EUCAST trừ CLSI (điểm phần trăm)",
                  fontsize=10)
    ax.set_title("Mức vênh giữa hai hệ ngưỡng theo cặp loài – kháng sinh",
                 fontsize=11.5, pad=12)
    ax.grid(axis="x", alpha=0.25, linestyle=":")
    ax.spines[["top", "right"]].set_visible(False)

    ax.text(0.99, -0.055,
            "Điểm bên phải: EUCAST báo cáo kháng nhiều hơn CLSI",
            transform=ax.transAxes, ha="right", va="top",
            fontsize=8, color="#666666")

    fig.tight_layout()
    fig.savefig(duong_dan, dpi=300, bbox_inches="tight")
    plt.close(fig)


# ============================================================
# CHAY
# ============================================================

def main():
    THU_MUC_RA.mkdir(exist_ok=True)

    print("\n[1/5] Doc va loc du lieu")
    d, nhat_ky = doc_tat_ca()
    print(f"  tong cong {len(d):,} cap chung - khang sinh co MIC")

    print("\n[2/5] Nap bang nguong")
    bp = nap_nguong()
    print(f"  {len(bp):,} dong nguong, "
          f"{bp[bp.he=='CLSI'].shape[0]} CLSI / {bp[bp.he=='EUCAST'].shape[0]} EUCAST")

    print("\n[3/5] Ap nguong")
    d = ap_hai_he(d, bp)

    co_ca_hai = (d["nhan_clsi"].notna() & d["nhan_eucast"].notna()).sum()
    print(f"  phan loai duoc bang ca hai he: {co_ca_hai:,} / {len(d):,} "
          f"({100*co_ca_hai/len(d):.1f}%)")
    chi_mot = ((d["nhan_clsi"].notna()) ^ (d["nhan_eucast"].notna())).sum()
    print(f"  chi mot he phan loai duoc: {chi_mot:,}  <- bao cao rieng trong bai")

    print("\n[4/5] Tinh ket qua")
    kq, ok = tinh_ket_qua(d)
    b1 = bang_mo_ta(d, ok)
    print(f"  {len(kq)} cap loai - khang sinh dat n >= {N_TOI_THIEU}")

    b1.to_csv(THU_MUC_RA / "bang1_mo_ta.csv", index=False)
    kq.to_csv(THU_MUC_RA / "bang2_ket_qua.csv", index=False)
    nhat_ky.to_csv(THU_MUC_RA / "nhat_ky_loc.csv", index=False)

    # Cap loai - khang sinh chi MOT he co nguong (du N_TOI_THIEU MIC)
    g = d.groupby(["loai", "khang_sinh"]).agg(
        n_MIC=("mic", "size"), co_clsi=("co_bp_clsi", "any"),
        co_eucast=("co_bp_eucast", "any")).reset_index()
    mot = g[(g["co_clsi"] != g["co_eucast"]) & (g["n_MIC"] >= N_TOI_THIEU)].copy()
    mot["chi_co"] = np.where(mot["co_clsi"], "CLSI", "EUCAST")
    mot = mot.sort_values(["chi_co", "loai", "khang_sinh"])
    mot[["loai", "khang_sinh", "n_MIC", "chi_co"]].to_csv(
        THU_MUC_RA / "chi_mot_he.csv", index=False)
    print(f"  cap chi mot he co nguong: {len(mot)} "
          f"(chi CLSI {(mot.chi_co=='CLSI').sum()}, chi EUCAST {(mot.chi_co=='EUCAST').sum()})")

    # Bang dinh dang san theo khuon Cusack 2019 CMI, de dan thang vao ban thao
    if len(kq):
        bang = pd.DataFrame({
            "Organism (no. isolates)": kq["loai"] + " (" + kq["n"].map("{:,}".format) + ")",
            "Antimicrobial agent": kq["khang_sinh"],
            "CLSI S": kq["clsi_S_pct"].round(1),
            "CLSI I": kq["clsi_I_pct"].round(1),
            "CLSI R": kq["clsi_R_pct"].round(1),
            "EUCAST S": kq["eucast_S_pct"].round(1),
            "EUCAST I": kq["eucast_I_pct"].round(1),
            "EUCAST R": kq["eucast_R_pct"].round(1),
            "Category agreement (%)": kq["do_khop_pct"].round(1),
        })
        bang.to_csv(THU_MUC_RA / "bang2_dinh_dang_bai_bao.csv", index=False)

    print("\n[5/5] Ve hinh 1")
    ve_hinh1(kq, THU_MUC_RA / "hinh1_forest.png")

    print(f"\nXong. Ket qua trong {THU_MUC_RA.resolve()}\n")

    if len(kq):
        print("10 cap venh nhieu nhat:")
        cot = ["loai", "khang_sinh", "n", "khang_clsi_pct",
               "khang_eucast_pct", "chenh_diem_pt", "ty_le_doi_pct"]
        print(kq[cot].head(10).to_string(index=False,
              float_format=lambda x: f"{x:.1f}"))


if __name__ == "__main__":
    main()
