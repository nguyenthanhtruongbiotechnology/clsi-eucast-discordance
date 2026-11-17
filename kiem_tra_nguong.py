#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Kiem tra file nguong.csv truoc khi chay phan tich.

Cach dung:
    python kiem_tra_nguong.py

Bat cac loi chep tay thuong gap. KHONG kiem tra duoc gia tri co dung
voi M100 hay EUCAST khong - chuyen do chi ban doi chieu lai duoc.
"""

from pathlib import Path
import pandas as pd
import numpy as np

FILE_NGUONG = Path("./nguong.csv")
THU_MUC_DU_LIEU = Path("./du_lieu")

# Thang pha loang doi, don vi mg/L
THANG = [2.0 ** k for k in range(-11, 13)]   # 0.000488 ... 4096

CAC_LOAI = {
    573: "Klebsiella pneumoniae", 590: "Salmonella enterica",
    1280: "Staphylococcus aureus", 287: "Pseudomonas aeruginosa",
    550: "Enterobacter cloacae",
    546: "Citrobacter freundii", 615: "Serratia marcescens",
    562: "Escherichia coli", 470: "Acinetobacter baumannii",
}

# Phai khop y het NHOM trong phan_tich_nguong.py
NHOM = {
    "Enterobacterales": [
        "Escherichia coli", "Klebsiella pneumoniae", "Salmonella enterica",
        "Enterobacter cloacae", "Citrobacter freundii", "Serratia marcescens",
    ],
}

# EUCAST dung gia tri nay khi khong co hang "nhay o lieu thuong":
# moi gia tri duoi nguong khang deu la "nhay khi tang lieu".
# Day la gia tri hop le, khong phai loi chep.
EUCAST_KHONG_CO_HANG_S = 0.001

loi, canh_bao = [], []


def tren_thang(x, dung_sai=0.08):
    """
    Gia tri co nam tren thang pha loang doi khong.

    Dung sai 8% vi tai lieu hay lam tron: 0.06 thay cho 0.0625,
    0.12 thay cho 0.125, 0.03 thay cho 0.03125. Day deu la gia tri hop le.
    """
    if pd.isna(x):
        return True
    return any(abs(x - m) / m < dung_sai for m in THANG)


def nap():
    if not FILE_NGUONG.exists():
        raise SystemExit(f"Khong thay {FILE_NGUONG}")
    bp = pd.read_csv(FILE_NGUONG, comment="#")
    bp.columns = [c.strip() for c in bp.columns]
    return bp


def kiem_cau_truc(bp):
    can = ["loai", "khang_sinh", "he", "nguong_nhay", "nguong_khang"]
    thieu = [c for c in can if c not in bp.columns]
    if thieu:
        loi.append(f"Thieu cot bat buoc: {thieu}")
        raise SystemExit("\n".join(loi))

    if "phien_ban" not in bp.columns:
        canh_bao.append(
            "Khong co cot 'phien_ban'. Bai bao phai ghi ro phien ban nguong, "
            "them cot nay di."
        )

    he_la = set(bp["he"].astype(str).str.strip().str.upper()) - {"CLSI", "EUCAST"}
    if he_la:
        loi.append(f"Cot 'he' co gia tri la: {he_la}. Chi duoc CLSI hoac EUCAST.")


def kiem_gia_tri(bp):
    for i, r in bp.iterrows():
        d = i + 2   # dong trong file, tinh ca dong tieu de
        s, k = r["nguong_nhay"], r["nguong_khang"]
        ks = r["khang_sinh"]

        if pd.isna(s) or pd.isna(k):
            loi.append(f"Dong {d} ({ks}): thieu nguong_nhay hoac nguong_khang")
            continue

        if s <= 0 or k <= 0:
            loi.append(f"Dong {d} ({ks}): nguong phai duong, dang co S={s} R={k}")
            continue

        # Loi nang nhat: nguong nhay lon hon nguong khang
        if s > k:
            loi.append(
                f"Dong {d} ({ks}): nguong_nhay ({s}) LON HON nguong_khang ({k}). "
                f"Gan nhu chac chan la dao nguoc hai cot."
            )

        # 0.001 la cach EUCAST ghi "khong co hang nhay o lieu thuong".
        # Bo qua moi kiem tra hinh thuc cho dong nay.
        if abs(s - EUCAST_KHONG_CO_HANG_S) < 1e-9:
            # Hop le cho EUCAST, va cho CLSI voi colistin/polymyxin (chi co I va R)
            if str(r.get("he", "")).strip().upper() != "EUCAST" and \
                    str(ks).strip().lower() not in ("colistin", "polymyxin b"):
                canh_bao.append(
                    f"Dong {d} ({ks}): nguong_nhay=0.001 (khong co hang nhay) "
                    f"o he {r.get('he')}. Kiem lai trong M100."
                )
            continue

        # Khoang cach qua lon thuong la nham bang
        if k / s > 64:
            canh_bao.append(
                f"Dong {d} ({ks}): khoang S={s} den R={k} cach nhau {k/s:.0f} lan. "
                f"Kiem lai co nham bang khong."
            )

        for ten, v in (("nguong_nhay", s), ("nguong_khang", k)):
            if not tren_thang(v):
                canh_bao.append(
                    f"Dong {d} ({ks}): {ten}={v} khong nam tren thang pha loang doi. "
                    f"Nguong lam sang gan nhu luon la luy thua cua 2."
                )


def kiem_cap_doi(bp):
    """
    Moi cap loai-khang sinh phai co du ca hai he moi so sanh duoc.

    Quy tac giong ma phan tich: dong ghi ro ten loai THANG dong de trong.
    Dong de trong ap cho MOI loai chua co dong rieng.
    """
    bp["_loai"] = bp["loai"].fillna("").astype(str).str.strip()
    bp["_ks"] = bp["khang_sinh"].astype(str).str.strip().str.lower()
    bp["_he"] = bp["he"].astype(str).str.strip().str.upper()

    moi_loai = set(CAC_LOAI.values())

    def bung(ten):
        """Ten nhom -> tap loai. Ten loai -> chinh no. Trong -> moi loai."""
        if ten == "":
            return moi_loai
        if ten in NHOM:
            return set(NHOM[ten])
        return {ten}

    # Tap loai thuc su duoc phu, cho tung (khang sinh, he)
    phu = {}
    for (ks, he), g in bp.groupby(["_ks", "_he"]):
        tap = set()
        for ten in g["_loai"]:
            tap |= bung(ten)
        phu[(ks, he)] = tap

    cac_ks = sorted(bp["_ks"].unique())
    for ks in cac_ks:
        c = phu.get((ks, "CLSI"), set())
        e = phu.get((ks, "EUCAST"), set())
        ca_hai = c & e
        chi_c = c - e
        chi_e = e - c

        if not c and e:
            canh_bao.append(f"'{ks}': chi co EUCAST, khong co dong CLSI nao "
                            f"-> bi loai khoi phan tich")
        elif not e and c:
            canh_bao.append(f"'{ks}': chi co CLSI, khong co dong EUCAST nao "
                            f"-> bi loai khoi phan tich")
        elif chi_c or chi_e:
            if chi_c:
                canh_bao.append(f"'{ks}': {len(chi_c)} loai co CLSI nhung thieu "
                                f"EUCAST -> {sorted(chi_c)}")
            if chi_e:
                canh_bao.append(f"'{ks}': {len(chi_e)} loai co EUCAST nhung thieu "
                                f"CLSI -> {sorted(chi_e)}")

        if ca_hai:
            print(f"  ok  {ks}: {len(ca_hai)}/{len(moi_loai)} loai co du ca hai he")

    # Dong trung lap that su: cung loai, cung khang sinh, cung he
    trung = bp.groupby(["_loai", "_ks", "_he"]).size()
    for (lo, ks, he), n in trung.items():
        if n > 1:
            ten = lo if lo else "(moi loai)"
            loi.append(f"'{ten} / {ks} / {he}' co {n} dong. Phai duy nhat mot dong.")


def kiem_khop_du_lieu(bp):
    """Ten khang sinh trong nguong.csv co khop voi du lieu BV-BRC khong."""
    if not THU_MUC_DU_LIEU.exists():
        canh_bao.append(f"Khong thay {THU_MUC_DU_LIEU}, bo qua buoc doi chieu ten.")
        return

    co_trong_dl = set()
    for tid in CAC_LOAI:
        f = THU_MUC_DU_LIEU / f"BVBRC_genome_amr_{tid}.csv"
        if not f.exists():
            continue
        try:
            d = pd.read_csv(f, usecols=["Antibiotic"], low_memory=False)
            co_trong_dl |= set(d["Antibiotic"].dropna().str.strip().str.lower())
        except Exception as e:
            canh_bao.append(f"Khong doc duoc {f.name}: {e}")

    if not co_trong_dl:
        return

    trong_bp = set(bp["khang_sinh"].astype(str).str.strip().str.lower())
    khong_khop = sorted(trong_bp - co_trong_dl)
    for ks in khong_khop:
        gan = [a for a in co_trong_dl if ks.split("/")[0][:6] in a]
        goi_y = f"  Co phai y la: {gan[:3]}" if gan else ""
        loi.append(f"'{ks}' khong co trong du lieu BV-BRC.{goi_y}")


def main():
    bp = nap()
    print(f"Doc {len(bp)} dong nguong tu {FILE_NGUONG}\n")

    kiem_cau_truc(bp)
    kiem_gia_tri(bp)
    kiem_cap_doi(bp)
    kiem_khop_du_lieu(bp)

    if loi:
        print(f"LOI ({len(loi)}) - phai sua truoc khi chay phan tich:")
        for x in loi:
            print(f"  x {x}")
        print()

    if canh_bao:
        print(f"CANH BAO ({len(canh_bao)}) - kiem lai cho chac:")
        for x in canh_bao:
            print(f"  ! {x}")
        print()

    if not loi and not canh_bao:
        print("Khong thay loi cau truc nao.")

    # Tom tat
    n_ks = bp["_ks"].nunique()
    n_dong = len(bp)
    print(f"Tom tat: {n_dong} dong, {n_ks} khang sinh.")
    print("\nLuu y: script nay KHONG kiem tra duoc gia tri co dung voi M100 "
          "hay EUCAST khong. Phan do ban phai tu doi chieu lai voi tai lieu goc.")

    if loi:
        raise SystemExit(1)


if __name__ == "__main__":
    main()
