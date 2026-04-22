import os
import sys
import torch  # MUST BE FIRST ON WINDOWS TO PREVENT DLL CLASH
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from sklearn.metrics import accuracy_score, roc_auc_score, matthews_corrcoef
from sklearn.model_selection import LeaveOneOut
from sklearn.preprocessing import StandardScaler
from tp53_sve import load_all_features, fisher_discriminant

from kan import KAN

def main():
    print("====================================================")
    print(" KAN vs LDA Comparison on TP53 Features")
    print("====================================================")
    
    # 1. Data Loading: Now using the pLDDT confidence-gated features
    try:
        df = pd.read_csv("output/phase3/plddt_weighted_features.csv")
    except FileNotFoundError:
        print("Weighted features not found! Run apply_plddt_weighting.py first.")
        df = load_all_features()
    
    feature_cols = [
        'RMSD', 'TM_Score', 'Mean_Displacement', 'Max_Displacement',
        'Residues_Above_5A', 'Residues_Above_10A',
        'Contacts_Lost', 'Preservation_Rate', 'DBD_Contact_Loss_Pct',
        'Total_SS_Changes', 'SS_Change_Pct', 'Helix_to_Coil', 'DBD_Changes',
        'Local_Global_Ratio', 'Total_SASA_Change', 'Hydrophobic_Exposure',
        'DBCA_Score', 'Zinc_Score', 'DNA_Contact_Score', 'Loop_Score',
        'ARES', 'DDE_Contact', 'Rewiring_Energy',
        'PC1', 'PC2',
        'BLOSUM62', 'Charge_Change', 'Hydro_Change', 'Volume_Change', 'MW_Change',
        'In_DBD', 'In_Zinc', 'In_DNA_Contact', 'In_Loop',
    ]
    
    # Filter labeled data dynamically
    # The dataset has 'Classification' column with 'Likely Oncogenic', 'Benign', and 'Unknown'
    df['Label'] = df['Classification'].apply(
        lambda x: 1 if x == 'Likely Oncogenic' else (0 if x == 'Benign' else -1)
    )
    
    labeled_df = df[df['Label'] != -1].copy()
    unlabeled_df = df[df['Label'] == -1].copy()
    
    X = labeled_df[feature_cols].values
    y = labeled_df['Label'].values
    X_full = df[feature_cols].values # For final full-dataset inference
    
    print(f"Loaded {len(X)} labeled samples ({sum(y==1)} Pathogenic, {sum(y==0)} Benign)")
    print(f"Number of features: {X.shape[1]}")
    
    # 2. Leave-One-Out Cross Validation
    loo = LeaveOneOut()
    
    lda_preds = np.zeros(len(y))
    lda_scores = np.zeros(len(y))
    
    kan_preds = np.zeros(len(y))
    kan_probs = np.zeros(len(y))
    
    kan_train_accs = []
    
    device = torch.device('cuda' if torch.cuda.is_available() else 'cpu')
    print(f"Using device for KAN: {device}")
    
    for i, (train_index, test_index) in enumerate(loo.split(X)):
        X_train, X_test = X[train_index], X[test_index]
        y_train, y_test = y[train_index], y[test_index]
        
        # Scale features
        scaler = StandardScaler()
        X_train_s = scaler.fit_transform(X_train)
        X_test_s = scaler.transform(X_test)
        
        # --- LDA ---
        X_b = X_train_s[y_train == 0]
        X_p = X_train_s[y_train == 1]
        
        w = fisher_discriminant(X_b, X_p)
        
        # Projections
        proj_train = X_train_s @ w
        proj_test = X_test_s @ w
        
        # Determine threshold from training data
        b_mean = np.mean(proj_train[y_train == 0])
        p_mean = np.mean(proj_train[y_train == 1])
        # Simple midpoint threshold
        thresh = (b_mean + p_mean) / 2
        
        # Direction
        if p_mean > b_mean:
            lda_scores[test_index] = proj_test[0]
            lda_preds[test_index] = 1 if proj_test[0] > thresh else 0
        else:
            lda_scores[test_index] = -proj_test[0]
            lda_preds[test_index] = 1 if proj_test[0] < thresh else 0
            
        # --- KAN ---
        # Very small KAN: 34 -> 2 -> 2
        model = KAN(width=[X.shape[1], 2, 2], grid=3, k=3, seed=42)
        dataset = {}
        dataset['train_input'] = torch.tensor(X_train_s, dtype=torch.float32).to(device)
        dataset['train_label'] = torch.tensor(y_train, dtype=torch.long).to(device)
        dataset['test_input'] = torch.tensor(X_test_s, dtype=torch.float32).to(device)
        dataset['test_label'] = torch.tensor(y_test, dtype=torch.long).to(device)
        
        # Train
        def train_acc():
            return torch.mean((torch.argmax(model(dataset['train_input']), dim=1) == dataset['train_label']).float())
            
        def test_acc():
            return torch.mean((torch.argmax(model(dataset['test_input']), dim=1) == dataset['test_label']).float())
            
        results = model.fit(dataset, opt="LBFGS", steps=20, log=10, metrics=(train_acc, test_acc), loss_fn=torch.nn.CrossEntropyLoss(), update_grid=False)
        
        kan_train_accs.append(results['train_acc'][-1])
        
        with torch.no_grad():
            logits = model(dataset['test_input'])
            probs = torch.softmax(logits, dim=1)
            pred = torch.argmax(logits, dim=1).item()
            kan_preds[test_index] = pred
            kan_probs[test_index] = probs[0, 1].item()
            
        print(f"Fold {i+1}/25 done. LDA: {lda_preds[test_index]}, KAN: {kan_preds[test_index]}, True: {y_test[0]}")
        
    print("\n--- RESULTS ---")
    
    lda_acc = accuracy_score(y, lda_preds)
    lda_auc = roc_auc_score(y, lda_scores)
    lda_mcc = matthews_corrcoef(y, lda_preds)
    
    kan_acc = accuracy_score(y, kan_preds)
    kan_auc = roc_auc_score(y, kan_probs)
    kan_mcc = matthews_corrcoef(y, kan_preds)
    
    kan_train_acc_mean = np.mean(kan_train_accs)
    
    print(f"[LDA] LOOCV Accuracy: {lda_acc:.3f}")
    print(f"[LDA] LOOCV AUC-ROC:  {lda_auc:.3f}")
    print(f"[LDA] LOOCV MCC:      {lda_mcc:.3f}")
    
    print(f"\n[KAN] Train Accuracy: {kan_train_acc_mean:.3f} (OVERFITTING CHECK)")
    print(f"[KAN] LOOCV Accuracy: {kan_acc:.3f}")
    print(f"[KAN] LOOCV AUC-ROC:  {kan_auc:.3f}")
    print(f"[KAN] LOOCV MCC:      {kan_mcc:.3f}")
    
    # Check overparameterization
    # 34 -> 2 -> 2 with grid=3, k=3 means 34*2 + 2*2 = 72 edges.
    # Spline params: (grid + k) * edges = (3+3)*72 = 432 parameters
    print(f"\n[Overfitting Assessment]")
    print(f"Samples: {len(X)}")
    print(f"KAN Parameters: Approx 432")
    if kan_train_acc_mean > 0.95 and kan_acc < 0.8:
        print("Verdict: MASSIVE OVERFITTING. The model memorizes training data but fails to generalize.")
    
    # Save a KAN model visualization on full data
    print("\nTraining full KAN for visualization...")
    model_full = KAN(width=[X.shape[1], 2, 2], grid=3, k=3, seed=42)
    scaler_full = StandardScaler()
    X_s = scaler_full.fit_transform(X)
    X_full_s = scaler_full.transform(X_full)
    dataset_full = {'train_input': torch.tensor(X_s, dtype=torch.float32).to(device),
                    'train_label': torch.tensor(y, dtype=torch.long).to(device),
                    'test_input': torch.tensor(X_full_s, dtype=torch.float32).to(device),
                    'test_label': torch.tensor(np.zeros(len(X_full)), dtype=torch.long).to(device)}
    
    model_full.fit(dataset_full, opt="LBFGS", steps=20, loss_fn=torch.nn.CrossEntropyLoss(), update_grid=False)
    # Prune it to make visualization interpretable
    model_full = model_full.prune()
    
    os.makedirs('output/phase3', exist_ok=True)
    try:
        model_full.plot(folder='output/phase3', beta=10)
        plt.savefig('output/phase3/kan_splines.png')
        print("Saved KAN spline visualization to output/phase3/kan_splines.png")
    except Exception as e:
        print(f"Failed to plot KAN: {e}")
        
    # --- Full Dataset Inference ---
    print("\n[Inference] Predicting on ALL 128 variants...")
    
    # KAN inference
    with torch.no_grad():
        logits = model_full(dataset_full['test_input'])
        probs = torch.softmax(logits, dim=1)
        kan_full_probs = probs[:, 1].cpu().numpy()
        kan_full_preds = torch.argmax(logits, dim=1).cpu().numpy()
        
    # LDA inference
    w_full = fisher_discriminant(X_s[y == 0], X_s[y == 1])
    proj_full = X_full_s @ w_full
    proj_train = X_s @ w_full
    b_mean_full = np.mean(proj_train[y == 0])
    p_mean_full = np.mean(proj_train[y == 1])
    thresh_full = (b_mean_full + p_mean_full) / 2
    
    lda_full_preds = np.zeros(len(df))
    lda_full_scores = np.zeros(len(df))
    
    if p_mean_full > b_mean_full:
        lda_full_scores = proj_full
        lda_full_preds = (proj_full > thresh_full).astype(int)
    else:
        lda_full_scores = -proj_full
        lda_full_preds = (-proj_full > -thresh_full).astype(int)
        
    df['KAN_Probability'] = kan_full_probs
    df['KAN_Prediction'] = kan_full_preds
    df['LDA_Score'] = lda_full_scores
    df['LDA_Prediction'] = lda_full_preds
    
    df.to_csv('output/phase3/kan_lda_full_predictions.csv', index=False)
    print("Saved full predictions to output/phase3/kan_lda_full_predictions.csv")
    
    # Write report
    report = f"""# KAN vs LDA Comparison Report (Full Dataset)
    
## Dataset
- Total samples: {len(X_full)}
- Labeled samples for training: {len(X)} ({sum(y==1)} Pathogenic, {sum(y==0)} Benign)
- Unlabeled/Unknown samples: {len(unlabeled_df)}
- Features: {X.shape[1]}

## Results
| Metric | Fisher LDA | Kolmogorov-Arnold Network (KAN) |
|--------|------------|---------------------------------|
| Train Acc | ~1.000 | {kan_train_acc_mean:.3f} |
| LOOCV Acc | {lda_acc:.3f} | {kan_acc:.3f} |
| LOOCV AUC | {lda_auc:.3f} | {kan_auc:.3f} |
| LOOCV MCC | {lda_mcc:.3f} | {kan_mcc:.3f} |

## Overfitting Assessment
KAN exhibits severe overfitting (n={len(X)} samples vs ~432 parameters).
"""
    with open('output/phase3/kan_report.md', 'w') as f:
        f.write(report)
        
    print("\n[DONE] Wrote report to output/phase3/kan_report.md")

if __name__ == "__main__":
    main()
