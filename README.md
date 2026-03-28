# DELkit — DNA-Encoded Library Analysis Pipeline

> **Complete cheminformatics pipeline from raw FASTQ to amino acid hits.**  
> Validated on PRJNA1215981 — a TfR1 peptide DEL published in PNAS February 2026.

---

## What this is

DELkit is an end-to-end DEL analysis pipeline built in a single Jupyter notebook. It takes raw NGS sequencing data (FASTQ) through every step of a real drug discovery workflow:

```
Raw FASTQ  →  Barcode decoding  →  Count matrix  →  Enrichment stats
     →  SAR mining  →  Chemical space  →  ML hit prediction  →  Amino acid hits
```

The pipeline was validated by independently reproducing the key findings of:

> He Q, Wang Y et al. *Enhanced Screening via a Pure DNA-Encoded Peptide Library Enabled by a Novel Fmoc Modification.* PNAS 2026. [DOI: 10.1073/pnas.2524999123](https://www.pnas.org/doi/10.1073/pnas.2524999123)

---

## Results — PRJNA1215981 (TfR1 Peptide DEL)

Starting from 49 million raw reads with no barcode reference:

| Metric | Value |
|---|---|
| Sequencing files processed | 3 (2 selection rounds + 1 input) |
| Total reads processed | ~49M |
| Decode rate | **89%** |
| Unique peptide combinations | 8,171,664 |
| Hits identified (FE ≥ 5×) | 26,662 |
| Top enrichment | **160.8×** |
| Primary pharmacophore identified | **D-Chg at position 3** |

**Top hit decoded:** `D-Chg — D-Chg — D-Chg — D-Chg` (160.8× enrichment)  
**Paper's best binder TR17:** `D-BrF — D-Chg — D-Chg — D-Chg — D-Lue` (Kd = 110 nM)

19 out of 20 top hits contain D-cyclohexylglycine (D-Chg) at position 3 — independently reproducing the paper's SAR conclusion that position 3 is the primary pharmacophore for TfR1 binding.

---

## Quick Start

### Option A — Run on bundled example data (no downloads needed)

```bash
git clone https://github.com/YOUR_USERNAME/delkit.git
cd delkit
pip install -r requirements.txt
jupyter lab DELkit_Analysis.ipynb
```

Run cells top to bottom. The example dataset (1,440 synthetic compounds, realistic DEL structure) is already in `data/`.

### Option B — Run on real PRJNA1215981 data

```bash
# 1. Install SRA toolkit (requires conda)
conda install -c bioconda sra-tools

# 2. Download raw sequencing data (~22GB)
prefetch PRJNA1215981
fasterq-dump --split-files --threads 4 --progress \
    SRR32132802 SRR32132803 SRR32132804 \
    SRR32132805 SRR32132806 SRR32132807

# 3. Move FASTQ files to project directory
mv SRR321328*.fastq /path/to/delkit/

# 4. Open notebook and set RUN_FASTQ = True in Cell 3b
jupyter lab DELkit_Analysis.ipynb
```

The barcode positions for PRJNA1215981 are pre-configured in the notebook:
```python
BARCODE_POSITIONS = [(19,28), (30,39), (41,50), (53,62)]
```
These were decoded from the read structure using positional conservation analysis — no reference required.

---

## Repository Structure

```
delkit/
├── DELkit_Analysis.ipynb        ← main notebook (12 cells, fully documented)
├── README.md
├── requirements.txt
├── LICENSE
│
├── data/
│   ├── count_file.csv           ← bundled example (1,440 compounds)
│   ├── cycle1_barcode_map.csv   ← barcode → building block, cycle 1
│   ├── cycle2_barcode_map.csv   ← barcode → building block, cycle 2
│   ├── cycle3_barcode_map.csv   ← barcode → building block, cycle 3
│   ├── library_design.csv       ← combinatorial topology
│   ├── ground_truth_labels.csv  ← for validation benchmarking
│   ├── selection_rep1.fastq     ← example NGS reads
│   ├── selection_rep2.fastq
│   └── bead_only_control.fastq
│
├── figures/                     ← generated plots
│   ├── volcano.png
│   ├── bb_enrichment.png
│   ├── chemical_space.png
│   └── feature_importance.png
│
└── outputs/                     ← generated analysis files
    ├── enrichment.csv
    ├── hits.csv
    ├── bb_enrichment_bb*.csv
    └── resynthesis_candidates.csv
```

---

## Notebook Walkthrough

| Cell | Description |
|---|---|
| **Cell 1** | Install dependencies |
| **Cell 2** | Imports and setup |
| **Cell 3** | Load data — change `COUNT_CSV` here for your own dataset |
| **Cell 3b** | FASTQ parser — set `RUN_FASTQ=True` for raw NGS input |
| **Cell 4** | Noise filter (NegBin background model) + CPM normalization |
| **Cell 5** | Enrichment calculation — fold enrichment, replicate consistency, composite score |
| **Cell 6** | Volcano plot |
| **Cell 7** | Hit selection and export |
| **Cell 8** | SAR: building block enrichment per cycle (Mann-Whitney significance) |
| **Cell 9** | Chemical space — UMAP/PCA with hit overlay |
| **Cell 10** | Random Forest hit classifier + active learning resynthesis prioritization |
| **Cell 11** | Summary |

---

## Supported Input Formats

| Format | Description | Key columns |
|---|---|---|
| **PB-DELASA** | Standard DEL count matrix | `BB1, BB2, BB3, sel_count_rep1, sel_count_rep2, ctrl_count` |
| **KinDEL** | insitro kinase benchmark | `seq_target_*, seq_matrix_*, molecule_hash` |
| **PRJNA1215981** | This repo's real data format | `bc1, bc2, bc3, bc4, sel_count_rep1, sel_count_rep2, ctrl_count` |
| **Generic** | Any CSV | `compound_id, sel_count, ctrl_count` |
| **FASTQ** | Raw NGS reads | Configured via `BARCODE_POSITIONS` in Cell 3b |

---

## Key Technical Decisions

### Why the NegBin model failed on PRJNA1215981
The input library control (SRR32132804) was sequenced at 10× lower depth than the selection files. Average control reads per compound = 1.4. Standard NegBin GLM flagged 99.6% of compounds as significant — a complete breakdown. The pipeline switches automatically to a composite enrichment score based on fold enrichment rank, absolute selection count rank, and replicate consistency rank.

### How barcode positions were decoded without a reference
Positional conservation analysis on 5,000 reads: positions with >80% base frequency = constant spacers, positions with ~25% each base = barcode. This revealed the structure `[19bp prefix][9bp BC1][AG][9bp BC2][GT][9bp BC3][GA+C][9bp BC4]` entirely from data.

### Why this is sequential selection, not replicates
SRR32132803 (round 1) and SRR32132802 (round 2) show 3× depth imbalance with near-zero replicate correlation. Round 2 represents surviving binders after a second stringent wash — not a replicate of round 1. The correct enrichment model compares round 2 directly to the input library.

---

## Dependencies

```
pandas >= 1.5
numpy >= 1.23
scipy >= 1.9
scikit-learn >= 1.1
matplotlib >= 3.6
seaborn >= 0.12
rdkit-pypi >= 2022.9
umap-learn >= 0.5
lightgbm >= 3.3
shap >= 0.41
statsmodels >= 0.13
```

Install all:
```bash
pip install -r requirements.txt
```

---

## Reproducing the PRJNA1215981 Results

The full analysis requires downloading ~22GB from NCBI SRA. Expected runtime on a modern laptop:

| Step | Time |
|---|---|
| FASTQ download | ~10 min (depends on connection) |
| fasterq-dump conversion | ~20 min |
| FASTQ parsing (3 files, 49M reads) | ~35 min |
| Enrichment + SAR + ML | ~15 min |
| **Total** | **~80 min** |

Key results to reproduce:
- Decode rate: **89%** on all 3 files
- Unique compounds: **8,171,664**
- Top hit: **D-Chg — D-Chg — D-Chg — D-Chg** at 160.8×
- SAR: **D-Chg at position 3** in 19/20 top hits

---

## Citation

If you use this pipeline, please cite the original paper whose data was used for validation:

```bibtex
@article{he2026enhanced,
  title={Enhanced Screening via a Pure DNA-Encoded Peptide Library 
         Enabled by a Novel Fmoc Modification},
  author={He, Qiujin and Wang, Yanhui and others},
  journal={Proceedings of the National Academy of Sciences},
  year={2026},
  doi={10.1073/pnas.2524999123}
}
```

---

## License

MIT License. See [LICENSE](LICENSE) for details.

---

## Contact

Built by Faraz Shaikh as a DEL cheminformatics portfolio project.  
Questions, issues, or collaboration: open a GitHub issue or connect on LinkedIn.
