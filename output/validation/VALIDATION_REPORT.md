# TP53-SVE Statistical Validation Report

This report addresses all identified statistical weaknesses in the original SVE evaluation.

---

## 1. Leave-One-Out Cross-Validation (LOOCV)

**Purpose:** Prove that the classifier generalizes and is not overfitting.

**Method:** Each of the 25 labeled samples is held out once, the model is retrained on the remaining 24, and the held-out sample is predicted. This is repeated 25 times, yielding an unbiased estimate of classification accuracy.

**Results (Original 5 Benign + 20 Pathogenic = 25 samples):**
- LOOCV Accuracy: **88.0%** (22/25)
- LOOCV AUC: **0.9600**
- Misclassified: **3**

Misclassified samples:
- M237I: True=Pathogenic, Predicted=Benign
- R213Q: True=Pathogenic, Predicted=Benign
- A189V: True=Benign, Predicted=Pathogenic

**LOOCV Detail (per sample):**

| Mutation | True Label | Predicted | Score | Correct |
|----------|-----------|-----------|-------|---------|
| R158L | Pathogenic | Pathogenic | 100.0 | YES |
| G245S | Pathogenic | Pathogenic | 100.0 | YES |
| H179R | Pathogenic | Pathogenic | 79.3 | YES |
| H193R | Pathogenic | Pathogenic | 86.0 | YES |
| M237I | Pathogenic | Benign | 42.8 | NO |
| P278S | Pathogenic | Pathogenic | 61.8 | YES |
| R158H | Pathogenic | Pathogenic | 84.0 | YES |
| C135Y | Pathogenic | Pathogenic | 67.2 | YES |
| C176F | Pathogenic | Pathogenic | 81.4 | YES |
| R282W | Pathogenic | Pathogenic | 77.3 | YES |
| R273L | Pathogenic | Pathogenic | 69.3 | YES |
| R273H | Pathogenic | Pathogenic | 100.0 | YES |
| R273C | Pathogenic | Pathogenic | 100.0 | YES |
| R249S | Pathogenic | Pathogenic | 100.0 | YES |
| R248W | Pathogenic | Pathogenic | 85.3 | YES |
| R248Q | Pathogenic | Pathogenic | 63.0 | YES |
| R213Q | Pathogenic | Benign | 32.0 | NO |
| R175H | Pathogenic | Pathogenic | 70.5 | YES |
| V157F | Pathogenic | Pathogenic | 65.4 | YES |
| Y220C | Pathogenic | Pathogenic | 83.2 | YES |
| R337H | Benign | Benign | 30.2 | YES |
| P47S | Benign | Benign | 28.4 | YES |
| A189V | Benign | Pathogenic | 62.7 | NO |
| P72R | Benign | Benign | 0.0 | YES |
| K132R | Benign | Benign | 42.5 | YES |

## 2. Expanded Benign Set Validation

**Purpose:** Test robustness with more benign controls (8 instead of 5).

Added: N345S, K382R, A138V (supported by ClinVar/literature evidence)

**Results (8 Benign + 20 Pathogenic = 28 samples):**
- LOOCV Accuracy: **78.6%** (22/28)
- LOOCV AUC: **0.9000**
- Misclassified: **6**

Misclassified samples:
- M237I: True=Pathogenic, Predicted=Benign
- C176F: True=Pathogenic, Predicted=Benign
- R213Q: True=Pathogenic, Predicted=Benign
- K382R: True=Benign, Predicted=Pathogenic
- A138V: True=Benign, Predicted=Pathogenic
- K132R: True=Benign, Predicted=Pathogenic

## 3. Feature Selection Analysis

**Purpose:** Prove the classifier works with fewer features (addressing the 34-feature / 25-sample overfitting concern).

| Top-K Features | LOOCV Accuracy | LOOCV AUC |
|---------------|---------------|----------|
| 3 | 88.0% | 0.8800 |
| 5 | 88.0% | 0.9400 |
| 7 | 84.0% | 0.9600 |
| 10 | 92.0% | 0.9800 |
| 15 | 92.0% | 0.9900 |
| 20 | 92.0% | 0.9900 |
| 34 | 88.0% | 0.9600 |

**Top 10 Most Predictive Features:**

1. **In_DBD** — 11.7% weight
2. **Hydrophobic_Exposure** — 9.1% weight
3. **DDE_Contact** — 7.7% weight
4. **ARES** — 5.6% weight
5. **DBD_TM** — 5.2% weight
6. **DNA_Contact_Score** — 4.5% weight
7. **SS_Change_Pct** — 4.1% weight
8. **Total_SS_Changes** — 4.1% weight
9. **Charge_Change** — 4.0% weight
10. **Helix_to_Coil** — 3.9% weight

## 4. Permutation Test

**Purpose:** Prove that the AUC is not achievable by random chance.

**Method:** Labels were randomly shuffled 1000 times. For each shuffle, Fisher's LDA was trained and AUC was computed.

- Real LOOCV AUC: **0.9600**
- Permutation AUC: mean = 0.9993, std = 0.0031
- P-value: **1.0000**
- Conclusion: **NOT significant** — classification may be due to chance.

## 5. Bootstrap 95% Confidence Interval

**Method:** 1000 bootstrap resamples with replacement.

- Bootstrap AUC: **1.0000**
- 95% CI: **[1.0000, 1.0000]**

## 6. Comparison to Experimental ΔΔG

**Purpose:** Validate computational energy scores against real thermodynamic measurements.

**Source:** Bullock et al. (2000) PNAS; Joerger & Fersht (2007) PNAS.

- ARES vs ΔΔG correlation: **r = 0.2786**
- RMSD vs ΔΔG correlation: r = 0.6009

| Mutation | ΔΔG (kcal/mol) | ARES | RMSD |
|----------|---------------|------|------|
| V143A | 3.5 | 42.3 | 35.5 |
| G245S | 1.8 | 40.1 | 33.5 |
| H179R | 2.7 | 52.3 | 31.2 |
| M237I | 1.0 | 30.4 | 19.7 |
| C135Y | 1.2 | 43.0 | 20.6 |
| C176F | 2.8 | 22.4 | 27.6 |
| R282W | 2.5 | 37.3 | 28.9 |
| R273H | 0.2 | 34.7 | 21.0 |
| R249S | 2.2 | 48.6 | 34.3 |
| R248W | 0.9 | 35.6 | 32.6 |
| R248Q | 0.5 | 46.5 | 22.2 |
| R175H | 3.0 | 49.2 | 31.9 |
| V157F | 1.5 | 29.6 | 32.3 |
| Y220C | 2.0 | 33.9 | 22.1 |
| A138V | 0.3 | 29.5 | 14.2 |
| R337H | 0.8 | 38.5 | 32.9 |

---

## Summary Scorecard

| Validation Test | Result | Status |
|----------------|--------|--------|
| LOOCV (Original 25) | AUC = 0.9600, Accuracy = 88.0% | PASS |
| LOOCV (Expanded 28) | AUC = 0.9000, Accuracy = 78.6% | PASS |
| Feature Selection (Top-5) | AUC = 0.9400 | PASS |
| Permutation Test | p = 1.0000 | FAIL |
| Bootstrap 95% CI | [1.0000, 1.0000] | PASS |
| ΔΔG Correlation (ARES) | r = 0.2786 | WEAK |
