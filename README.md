# CLSI vs EUCAST breakpoint discordance in public MIC data

Code for the study *Discordance between CLSI and EUCAST breakpoints across nine bacterial species: a reinterpretation of publicly available minimum inhibitory concentration data* (submitted to Western Pacific Surveillance and Response).

The analysis takes laboratory-measured minimum inhibitory concentrations (MICs) from the [Bacterial and Viral Bioinformatics Resource Center (BV-BRC)](https://www.bv-brc.org), interprets each MIC with both CLSI M100 (36th edition, 2026) and EUCAST (version 16.1, 2026) breakpoints, and compares the two interpretations within the same isolates, overall and by WHO region of isolation.

## Contents

| File | Purpose |
|---|---|
| `chan_doan.py` | Diagnostic: counts laboratory-measured vs model-predicted MICs in each download |
| `kiem_tra_nguong.py` | Checks the breakpoint table for structural errors before analysis |
| `phan_tich_nguong.py` | Main analysis: filtering, interpretation, agreement and paired differences |
| `phan_tich_khu_vuc.py` | Links country of isolation and repeats the analysis by WHO region |
| `ve_hinh1_heatmap.py` | Draws Figure 1 |
| `ten_thuoc.csv` | Antimicrobial name standardization (merged and excluded names) |
| `nguong_EUCAST.csv` | EUCAST v16.1 breakpoints used |
| `nguong_CLSI_template.csv` | CLSI rows used, **values left blank** (see below) |
| `requirements.txt` | Python packages |

Script comments are in Vietnamese without diacritics.

## CLSI breakpoints

CLSI M100 is a copyrighted document, so its breakpoint values are not distributed here. `nguong_CLSI_template.csv` lists every species–agent row used, with empty `nguong_nhay` (susceptible) and `nguong_khang` (resistant) columns. Fill them from your own copy of M100 (36th edition), then combine with the EUCAST file:

```bash
python -c "import pandas as pd; pd.concat([pd.read_csv('nguong_EUCAST.csv'), pd.read_csv('nguong_CLSI_template.csv')]).to_csv('nguong.csv', index=False)"
```

Conventions in the breakpoint table:

- An isolate is susceptible if MIC ≤ `nguong_nhay` and resistant if MIC > `nguong_khang`; anything between is intermediate. A CLSI resistant breakpoint published as "≥ x" is entered as the doubling dilution below x.
- `nguong_nhay = 0.001` means the system has no susceptible category (EUCAST "susceptible, increased exposure" agents; CLSI colistin).
- `loai = Enterobacterales` applies the row to all six Enterobacterales species; a row naming a species overrides it.

## Data

Data are not included. To reproduce:

1. For each taxon (562 *E. coli*, 573 *K. pneumoniae*, 590 *S. enterica*, 550 *E. cloacae*, 546 *C. freundii*, 615 *S. marcescens*, 287 *P. aeruginosa*, 470 *A. baumannii*, 1280 *S. aureus*), open `https://www.bv-brc.org/view/Taxonomy/<taxon>#view_tab=amr`, apply the filter **Evidence = Laboratory Method** and download as CSV to `du_lieu/BVBRC_genome_amr_<taxon>.csv`. The filter is required for *E. coli*, which otherwise exceeds the 2 500 000-record download limit; most MIC records in BV-BRC are machine-learning predictions and are excluded by the analysis in any case.
2. For the regional analysis, download the Genomes table of each taxon (with the Isolation Country and Collection Year columns) to `du_lieu/bo_gen/`.

## Running

```bash
python -m venv .venv
.venv\Scripts\activate          # Windows; on macOS/Linux: source .venv/bin/activate
pip install -r requirements.txt

python chan_doan.py
python kiem_tra_nguong.py
python phan_tich_nguong.py
python ve_hinh1_heatmap.py
python phan_tich_khu_vuc.py
```

Outputs are written to `ket_qua/`.

## Citation

Nguyễn Thành Trương. Discordance between CLSI and EUCAST breakpoints across nine bacterial species: a reinterpretation of publicly available minimum inhibitory concentration data. Submitted, 2026.

## License

Code: MIT License. EUCAST breakpoint values are reproduced from the freely available EUCAST breakpoint tables (https://www.eucast.org).
