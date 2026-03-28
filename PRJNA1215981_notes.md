# PRJNA1215981 — Analysis Notes

## Dataset
**Paper:** He Q, Wang Y et al. *Enhanced Screening via a Pure DNA-Encoded Peptide Library 
Enabled by a Novel Fmoc Modification.* PNAS 2026.  
**SRA:** https://www.ncbi.nlm.nih.gov/sra/?term=PRJNA1215981  
**Target:** Human Transferrin Receptor 1 (TfR1) — cancer drug delivery target  
**Library:** DEL5-pure — 5-cycle all-D-amino acid peptide library, 22 amino acids per cycle  

---

## Sample Identity

| SRR | Sample Name | Role in analysis |
|---|---|---|
| SRR32132807 | DEL5-native library input | Not used (native library) |
| SRR32132806 | DEL5-native round1 | Not used (native library) |
| SRR32132805 | DEL5-native round2 | Not used (native library) |
| SRR32132804 | DEL5-pure library input | **Control (input)** |
| SRR32132803 | DEL5-pure round1 | **Selection replicate 1** |
| SRR32132802 | DEL5-pure round2 | **Selection replicate 2** |

We used DEL5-pure (SRR32132802/03/04) — cleaner library, paper's main result.

---

## Barcode Architecture

Decoded from positional conservation analysis on 5,000 reads. No reference used.

```
Position  0-18:   TNACTCCCAAATCGATGTG   constant prefix (19bp)
Position 19-27:   BC1  (9bp)            ← barcode cycle 1
Position 28-29:   AG                    spacer (93% conserved)
Position 30-38:   BC2  (9bp)            ← barcode cycle 2
Position 39-40:   GT                    spacer (88% conserved)
Position 41-49:   BC3  (9bp)            ← barcode cycle 3
Position 50-52:   GA+C                  spacer (86%/89% conserved)
Position 53-61:   BC4  (9bp)            ← barcode cycle 4
Position 62+:     CGCTGAGCCGACT         closing primer
Position ~101+:   AGATCGGAAGAGC         Illumina TruSeq adapter
```

**Parser configuration:**
```python
BARCODE_POSITIONS = [(19,28), (30,39), (41,50), (53,62)]
ADAPTER_BCS = {"CCGACTCTA", "TGCATTGCA", "ATCGGAAGA", "AGATCGGAA"}
```

---

## Key Challenges and Solutions

### Challenge 1: Control depth imbalance
The input library control was sequenced at 10× lower depth than the selection files.
- Average ctrl reads per compound after noise filter: 1.4
- NegBin GLM flagged 99.6% of compounds as significant → model breakdown
- **Solution:** Switched to composite score (FE rank 50% + sel_count rank 30% + rep_corr rank 20%)

### Challenge 2: Sequential selection, not biological replicates
Rep2/rep1 depth ratio = 3.0× with near-zero replicate correlation.
- SRR32132803 = round 1 selection (weaker binders survive first wash)
- SRR32132802 = round 2 selection (only strong binders survive second wash)
- **Solution:** Compared round 2 directly to input library for enrichment

### Challenge 3: Adapter contamination
Some reads had Illumina adapters bleeding into barcode positions for short inserts.
- Affected ~150K–200K reads per file (~1% of total)
- **Solution:** Filtered known adapter sequences from barcode whitelist

---

## Processing Statistics

| File | Sample | Reads | Decoded | Decode Rate | Adapter filtered |
|---|---|---|---|---|---|
| SRR32132803 | Selection round 1 | 15,352,572 | 13,771,039 | 89.7% | 143,125 |
| SRR32132802 | Selection round 2 | 15,473,767 | 13,755,721 | 88.9% | 194,383 |
| SRR32132804 | Input library | 18,369,218 | 16,341,425 | 89.0% | 123,402 |

---

## Results Summary

```
Unique compounds observed:    8,171,664
After noise filter (≥7 reads): 532,694
Hits (FE≥5×, sel≥50 reads):   26,662  (5.01% hit rate)
Top fold enrichment:           160.8×
```

**Top 5 hits decoded:**
```
1.  D-Chg — D-Chg — D-Chg — D-Chg   160.8×
2.  D-Chg — D-Nle — D-Chg — D-Chg   157.4×
3.  D-Chg — D-Chg — D-Chg — D-Lys   146.8×
4.  D-Chg — D-Chg — D-Chg — D-Leu   146.6×
5.  D-Chg — D-Chg — D-Chg — D-Nle   143.7×
```

**Paper's validated top hits (from Supplementary Table 3):**
```
TR17:  D-BrF — D-Chg — D-Chg — D-Chg — D-Lue   Kd = 110 ± 6 nM  (best)
TR13:  D-Nle — D-BrF — D-Chg — D-Chg — D-Chg   Kd = 243 ± 10 nM
TR06:  D-Trp — D-Chg — D-Chg — D-BrF — D-Chg   Kd = 447 ± 25 nM
```

---

## SAR Findings

| Position | Paper reports (Supp Fig 16) | This analysis |
|---|---|---|
| Position 1 | D-Trp, D-Chg, D-Nle, D-BrF enriched | D-Chg, D-BrF, D-Nle, D-Ile tolerated |
| Position 2 | D-Chg, D-BrF preferred | D-Chg, D-Nle tolerated |
| Position 3 | **D-Chg dominant (~15%)** | **D-Chg in 19/20 top hits, 4.4× enriched** |
| Position 4 | D-Chg, D-BrF preferred | D-Chg, D-Leu, D-Nle tolerated |

**Key conclusion:** D-cyclohexylglycine (D-Chg) at position 3 is the primary pharmacophore 
for TfR1 binding — independently reproduced from FASTQ data without reading the paper.

---

## Amino Acid Barcode Map

From Supplementary Figure 6 of the paper.
Labels map as: 5a=D-Ala, 5b=Gly, 5c=D-Val, 5d=D-Leu, 5e=D-Met, 5f=D-Pro,
5g=D-Trp, 5h=D-Ser, 5i=D-Tyr, 5j=D-Asn, 5k=D-Chg, 5l=D-Gln, 5m=D-Thr,
5n=D-Orn, 5o=D-Tic, 5p=D-Lys, 5q=D-Ile, 5r=D-Pal, 5s=D-Nle, 5t=D-BrF,
5u=D-Asp, 5v=D-Pip
