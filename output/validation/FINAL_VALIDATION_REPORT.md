# TP53-SVE Final Validation Report

## Model Configuration

- **Features:** 10 (reduced from 34 to prevent overfitting)
- **Feature Set:** In_DBD, Hydrophobic_Exposure, DDE_Contact, ARES, DBD_TM, DNA_Contact_Score, SS_Change_Pct, Total_SS_Changes, Charge_Change, Helix_to_Coil
- **Classifier:** Fisher's Linear Discriminant Analysis
- **Training Set:** 20 pathogenic + 5 benign = 25 labeled

---

## Validation Scorecard

| Test | Metric | Result | Pass? |
|------|--------|--------|-------|
| LOOCV (34 features) | AUC / Accuracy | 0.9300 / 92.0% | Baseline |
| LOOCV (10 features) | AUC / Accuracy | 0.9900 / 92.0% | ✅ |
| Permutation (train/test) | p-value | 0.0620 | ❌ |
| Permutation (LOOCV) | p-value | 0.0000 | ✅ |
| Bootstrap 95% CI | AUC range | [1.0000, 1.0000] | ✅ |
| External (pseudo) | Agreement | 11/18 (61%) | ❌ |
| Delta Delta G (ARES) | Pearson r | 0.2786 | ⚠️ |
| Delta Delta G (RMSD) | Pearson r | 0.6009 | ✅ |

## Detailed LOOCV Results (10-Feature Model)

| Mutation | True | Predicted | Score | Correct |
|----------|------|-----------|-------|---------|
| R213Q | Pathogenic | Benign | 40.4 | ❌ |
| P278S | Pathogenic | Pathogenic | 78.7 | ✅ |
| R248W | Pathogenic | Pathogenic | 76.0 | ✅ |
| R249S | Pathogenic | Pathogenic | 74.8 | ✅ |
| R175H | Pathogenic | Pathogenic | 72.1 | ✅ |
| R158L | Pathogenic | Pathogenic | 83.6 | ✅ |
| R158H | Pathogenic | Pathogenic | 83.6 | ✅ |
| H193R | Pathogenic | Pathogenic | 82.5 | ✅ |
| H179R | Pathogenic | Pathogenic | 50.9 | ✅ |
| G245S | Pathogenic | Pathogenic | 100.0 | ✅ |
| C176F | Pathogenic | Pathogenic | 53.6 | ✅ |
| C135Y | Pathogenic | Pathogenic | 58.1 | ✅ |
| M237I | Pathogenic | Pathogenic | 63.1 | ✅ |
| R248Q | Pathogenic | Pathogenic | 88.8 | ✅ |
| Y220C | Pathogenic | Pathogenic | 66.8 | ✅ |
| R273L | Pathogenic | Pathogenic | 68.0 | ✅ |
| R273H | Pathogenic | Pathogenic | 100.0 | ✅ |
| R282W | Pathogenic | Pathogenic | 71.6 | ✅ |
| R273C | Pathogenic | Pathogenic | 80.5 | ✅ |
| V157F | Pathogenic | Pathogenic | 66.3 | ✅ |
| R337H | Benign | Benign | 17.2 | ✅ |
| P47S | Benign | Benign | 9.1 | ✅ |
| A189V | Benign | Pathogenic | 50.3 | ❌ |
| K132R | Benign | Benign | 29.9 | ✅ |
| P72R | Benign | Benign | 0.0 | ✅ |

## Pseudo-External Validation (25 Unlabeled Variants)

| Mutation | Score | Prediction | Known Classification |
|----------|-------|------------|---------------------|
| E221K | 100.0 | Pathogenic | Unknown |
| K164E | 98.4 | Pathogenic | Unknown |
| G245D | 91.2 | Pathogenic | Unknown |
| H179Q | 85.0 | Pathogenic | Unknown |
| G244C | 82.5 | Pathogenic | Unknown |
| L194R | 82.4 | Pathogenic | Likely Oncogenic (COSMIC) |
| K292R | 81.0 | Pathogenic | Unknown |
| R175G | 78.9 | Pathogenic | Unknown |
| T253I | 78.7 | Pathogenic | Unknown |
| D281G | 76.3 | Pathogenic | Likely Oncogenic (COSMIC) |
| R249T | 75.4 | Pathogenic | Unknown |
| E285K | 75.3 | Pathogenic | Likely Oncogenic (COSMIC) |
| P278R | 74.8 | Pathogenic | Unknown |
| W146C | 73.7 | Pathogenic | Unknown |
| C242F | 73.3 | Pathogenic | Unknown |
| Y220N | 72.5 | Pathogenic | Unknown |
| C238Y | 71.6 | Pathogenic | Unknown |
| V173L | 70.0 | Pathogenic | Unknown |
| C176R | 69.7 | Pathogenic | Unknown |
| R156H | 69.6 | Pathogenic | Unknown |
| N247D | 69.5 | Pathogenic | Likely Oncogenic (COSMIC) |
| I195T | 69.4 | Pathogenic | Likely Oncogenic (COSMIC) |
| L130F | 69.0 | Pathogenic | Unknown |
| C141Y | 68.7 | Pathogenic | Unknown |
| R282Q | 68.1 | Pathogenic | Unknown |
| P152L | 68.1 | Pathogenic | Unknown |
| R267W | 68.0 | Pathogenic | Unknown |
| G187T | 67.8 | Pathogenic | Unknown |
| Y126D | 66.8 | Pathogenic | Unknown |
| C242S | 66.6 | Pathogenic | Unknown |
| V143A | 65.6 | Pathogenic | Temperature-Sensitive |
| D49R | 65.4 | Pathogenic | Unknown |
| R249G | 65.3 | Pathogenic | Unknown |
| R306P | 63.5 | Pathogenic | Unknown |
| Y220S | 62.8 | Pathogenic | Unknown |
| F113N | 62.3 | Pathogenic | Unknown |
| N239D | 61.8 | Pathogenic | Likely Oncogenic (COSMIC) |
| N210S | 61.2 | Pathogenic | Unknown |
| S166L | 61.1 | Pathogenic | Unknown |
| Y234C | 60.5 | Pathogenic | Unknown |
| R196Q | 59.7 | Pathogenic | Unknown |
| A159V | 59.6 | Pathogenic | Unknown |
| P250L | 58.9 | Pathogenic | Unknown |
| R110P | 58.7 | Pathogenic | Unknown |
| E286K | 58.7 | Pathogenic | Unknown |
| V31I | 58.2 | Pathogenic | Unknown |
| G266V | 56.1 | Pathogenic | Unknown |
| R248G | 55.7 | Pathogenic | Unknown |
| R280K | 54.8 | Pathogenic | Likely Oncogenic (COSMIC) |
| T387V | 54.7 | Pathogenic | Unknown |
| F270L | 54.5 | Pathogenic | Unknown |
| I251F | 53.2 | Pathogenic | Unknown |
| A138V | 52.5 | Pathogenic | Likely Benign (literature) |
| K319T | 50.9 | Pathogenic | Unknown |
| R280T | 49.8 | Pathogenic | Unknown |
| V272M | 49.5 | Pathogenic | Gain-of-Function |
| R175C | 49.3 | Pathogenic | Unknown |
| E294G | 48.6 | Pathogenic | Unknown |
| T155I | 48.1 | Pathogenic | Unknown |
| V217M | 47.6 | Pathogenic | Unknown |
| S240N | 47.3 | Pathogenic | Unknown |
| S366A | 47.2 | Pathogenic | Unknown |
| G108L | 46.7 | Pathogenic | Unknown |
| K382R | 46.6 | Pathogenic | Likely Benign (literature) |
| R248L | 46.5 | Pathogenic | Unknown |
| Y163C | 46.3 | Pathogenic | Unknown |
| Q5T | 46.2 | Pathogenic | Unknown |
| C176S | 45.1 | Pathogenic | Unknown |
| S241F | 44.7 | Pathogenic | Likely Oncogenic (COSMIC) |
| D281E | 41.7 | Pathogenic | Unknown |
| L14V | 41.1 | Pathogenic | Unknown |
| G302R | 40.8 | Pathogenic | Unknown |
| S313R | 40.2 | Pathogenic | Unknown |
| E258K | 40.1 | Pathogenic | Unknown |
| A276V | 39.9 | Pathogenic | Unknown |
| H380R | 37.7 | Pathogenic | Unknown |
| S378Q | 37.5 | Pathogenic | Unknown |
| P36G | 37.2 | Pathogenic | Unknown |
| E56D | 37.2 | Pathogenic | Unknown |
| Q331R | 35.2 | Benign | Unknown |
| R342P | 34.9 | Benign | Likely Oncogenic (COSMIC) |
| S46E | 34.2 | Benign | Unknown |
| E339D | 34.0 | Benign | Unknown |
| D393E | 32.1 | Benign | Unknown |
| E352D | 31.9 | Benign | Unknown |
| K372S | 31.8 | Benign | Unknown |
| Q375H | 31.8 | Benign | Unknown |
| F341K | 28.2 | Benign | Unknown |
| M40N | 28.1 | Benign | Unknown |
| T125M | 27.1 | Benign | Likely Oncogenic (COSMIC) |
| L25V | 26.8 | Benign | Unknown |
| W23R | 26.5 | Benign | Likely Oncogenic (COSMIC) |
| P34L | 25.6 | Benign | Unknown |
| D21N | 24.4 | Benign | Unknown |
| N345S | 22.6 | Benign | Likely Benign (literature) |
| E326N | 21.7 | Benign | Unknown |
| W91L | 17.4 | Benign | Unknown |
| F54S | 12.8 | Benign | Unknown |
| L344R | 11.5 | Benign | Likely Oncogenic (COSMIC) |
| P85S | 11.4 | Benign | Unknown |
| N29S | 6.5 | Benign | Unknown |
| P71V | 3.9 | Benign | Unknown |
| L22F | 0.7 | Benign | Likely Oncogenic (COSMIC) |

**Agreement with known classifications: 11/18 (61%)**

---

## What Is Still Needed for Full Publication Readiness

### Current Dataset Composition

| Category | Count | Details |
|----------|------|---------|
| Pathogenic (labeled) | 20 | IARC TP53 Database confirmed hotspots |
| Benign (labeled) | 5 | ClinVar verified polymorphisms |
| Unlabeled (used as pseudo-external) | 103 | Mixed classifications |
| **Total with AlphaFold structures** | **128** | |

### Minimum Requirements for Tier-1 Publication

| Requirement | Current | Needed | Gap |
|------------|---------|--------|-----|
| Benign controls | 5 | 20-30 | +15-25 new benign variants |
| Total labeled samples | 25 | 60-80 | +35-55 new labeled variants |
| Independent test set | 0 (pseudo only) | 20-30 | Completely new variants not in training |
| AlphaFold structures needed | 0 new | 50-80 new | Must run AlphaFold for each new variant |
| Experimental Delta Delta G comparison | 16 mutations | 20-30 | +4-14 more published Delta Delta G values |

### Specific Benign Variants to Add (from ClinVar/gnomAD)

- E11Q - Benign (rs28934873)
- D21N - Benign polymorphism
- P34L - Benign (rs77697176)
- V31I - Benign polymorphism
- A63T - Likely Benign (rs764060847)
- P85S - Likely Benign
- R110L - Benign (rs78378222 region)
- S127T - Likely Benign
- Q144R - Likely Benign
- N210S - Likely Benign
- R267W - Likely Benign
- E294G - Likely Benign
- G302R - Likely Benign
- K319T - Likely Benign
- R342L - Likely Benign
- S366A - Likely Benign
- Q375H - Likely Benign
- Q331R - Likely Benign
- E339D - Likely Benign
- A347T - Likely Benign

### Specific Pathogenic Variants to Add for External Test

- R110P - Pathogenic (not in training set)
- Y163C - Pathogenic hotspot
- V173L - Pathogenic
- C242S - Zinc ligand disruption
- G244D - L3 loop
- R196* - Likely pathogenic
- P250L - Pathogenic
- E258K - Pathogenic
- R267P - Pathogenic
- S215G - Likely pathogenic

### Steps to Complete

1. **Generate AlphaFold structures** for 30-50 new variants (requires AlphaFold server access)
2. **Run the full pipeline** (Phases 1-3) on new structures
3. **Split into train/test**: Original 50 = training, New 30-50 = independent test
4. **Report independent test AUC** — this is the number a reviewer cares about
5. **Run LOOCV-permutation test** on the combined larger dataset
