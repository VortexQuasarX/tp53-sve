# KAN vs LDA Comparison Report (Full Dataset)
    
## Dataset
- Total samples: 128
- Labeled samples for training: 67 (56 Pathogenic, 11 Benign)
- Unlabeled/Unknown samples: 61
- Features: 34

## Results
| Metric | Fisher LDA | Kolmogorov-Arnold Network (KAN) |
|--------|------------|---------------------------------|
| Train Acc | ~1.000 | 1.000 |
| LOOCV Acc | 0.881 | 0.866 |
| LOOCV AUC | 0.890 | 0.854 |
| LOOCV MCC | 0.598 | 0.529 |

## Overfitting Assessment
KAN exhibits severe overfitting (n=67 samples vs ~432 parameters).
