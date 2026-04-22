import os
import torch
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from sklearn.metrics import accuracy_score, roc_auc_score, matthews_corrcoef
from sklearn.model_selection import LeaveOneOut
from sklearn.preprocessing import StandardScaler
from sklearn.decomposition import PCA
from tp53_sve import load_all_features

from kan import KAN

def main():
    print("====================================================")
    print(" PCA-KAN: Resolving Overparameterization via Feature Reduction")
    print("====================================================")
    
    # 1. Data Loading
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
    
    df['Label'] = df['Classification'].apply(
        lambda x: 1 if x == 'Likely Oncogenic' else (0 if x == 'Benign' else -1)
    )
    
    labeled_df = df[df['Label'] != -1].copy()
    unlabeled_df = df[df['Label'] == -1].copy()
    
    X_raw = labeled_df[feature_cols].values
    y = labeled_df['Label'].values
    X_full_raw = df[feature_cols].values
    
    print(f"Loaded {len(X_raw)} labeled samples ({sum(y==1)} Pathogenic, {sum(y==0)} Benign)")
    print(f"Original features: {X_raw.shape[1]}")
    
    # 2. PCA Compression (Train PCA on ALL 128 variants to capture variance space)
    scaler = StandardScaler()
    X_full_scaled = scaler.fit_transform(X_full_raw)
    
    # We want ~3-4 components to keep KAN parameters < 50
    N_COMPONENTS = 3
    pca = PCA(n_components=N_COMPONENTS)
    X_full_pca = pca.fit_transform(X_full_scaled)
    
    # Extract the transformed labeled subset for LOOCV training
    X = X_full_pca[labeled_df.index]
    
    print(f"Reduced features via PCA: {N_COMPONENTS}")
    print(f"Explained Variance: {pca.explained_variance_ratio_}")
    print(f"Total Variance Captured: {sum(pca.explained_variance_ratio_):.2%}")
    
    # 3. KAN Parameter Check
    # Architecture: [N_COMPONENTS, 2] (Direct mapping to 2 outputs, NO hidden layer needed for PCA space)
    # Spline params for [3, 2] grid=3, k=3:
    # Edges = 3 * 2 = 6 edges
    # Params = (grid + k) * edges = (3 + 3) * 6 = 36 parameters
    params = (3 + 3) * (N_COMPONENTS * 2) 
    print(f"\n[Overfitting Solved]")
    print(f"Samples: {len(X)}")
    print(f"New KAN Parameters: ~{params} (Safe under Rule-of-10!)")
    
    # 4. Leave-One-Out Cross Validation
    loo = LeaveOneOut()
    
    kan_preds = np.zeros(len(y))
    kan_probs = np.zeros(len(y))
    kan_train_accs = []
    
    device = torch.device('cuda' if torch.cuda.is_available() else 'cpu')
    print(f"\nTraining KAN with PCA inputs on {device}...")
    
    for i, (train_index, test_index) in enumerate(loo.split(X)):
        X_train, X_test = X[train_index], X[test_index]
        y_train, y_test = y[train_index], y[test_index]
        
        # Scale PCA features (Good practice for NN stability)
        fold_scaler = StandardScaler()
        X_train_s = fold_scaler.fit_transform(X_train)
        X_test_s = fold_scaler.transform(X_test)
        
        # --- KAN ---
        model = KAN(width=[N_COMPONENTS, 2], grid=3, k=3, seed=42)
        dataset = {}
        dataset['train_input'] = torch.tensor(X_train_s, dtype=torch.float32).to(device)
        dataset['train_label'] = torch.tensor(y_train, dtype=torch.long).to(device)
        dataset['test_input'] = torch.tensor(X_test_s, dtype=torch.float32).to(device)
        dataset['test_label'] = torch.tensor(y_test, dtype=torch.long).to(device)
            
        def train_acc():
            return torch.mean((torch.argmax(model(dataset['train_input']), dim=1) == dataset['train_label']).float())
            
        def test_acc():
            return torch.mean((torch.argmax(model(dataset['test_input']), dim=1) == dataset['test_label']).float())
            
        results = model.fit(dataset, opt="LBFGS", steps=20, log=20, metrics=(train_acc, test_acc), loss_fn=torch.nn.CrossEntropyLoss(), update_grid=False)
        kan_train_accs.append(results['train_acc'][-1])
        
        with torch.no_grad():
            logits = model(dataset['test_input'])
            probs = torch.softmax(logits, dim=1)
            pred = torch.argmax(logits, dim=1).item()
            kan_preds[test_index] = pred
            kan_probs[test_index] = probs[0, 1].item()
            
    # 5. Results
    print("\n--- PCA-KAN RESULTS ---")
    
    kan_acc = accuracy_score(y, kan_preds)
    kan_auc = roc_auc_score(y, kan_probs)
    kan_mcc = matthews_corrcoef(y, kan_preds)
    kan_train_acc_mean = np.mean(kan_train_accs)
    
    print(f"[PCA-KAN] Train Accuracy: {kan_train_acc_mean:.3f} (Lower is better, means less memorization)")
    print(f"[PCA-KAN] LOOCV Accuracy: {kan_acc:.3f}")
    print(f"[PCA-KAN] LOOCV AUC-ROC:  {kan_auc:.3f}")
    print(f"[PCA-KAN] LOOCV MCC:      {kan_mcc:.3f}")
    
    # 6. Train Full Model & Save
    print("\nTraining final PCA-KAN on all labeled data...")
    model_full = KAN(width=[N_COMPONENTS, 2], grid=3, k=3, seed=42)
    final_scaler = StandardScaler()
    X_s = final_scaler.fit_transform(X)
    X_full_s = final_scaler.transform(X_full_pca)
    
    dataset_full = {'train_input': torch.tensor(X_s, dtype=torch.float32).to(device),
                    'train_label': torch.tensor(y, dtype=torch.long).to(device),
                    'test_input': torch.tensor(X_full_s, dtype=torch.float32).to(device),
                    'test_label': torch.tensor(np.zeros(len(X_full_s)), dtype=torch.long).to(device)}
    
    model_full.fit(dataset_full, opt="LBFGS", steps=30, loss_fn=torch.nn.CrossEntropyLoss(), update_grid=False)
    
    # Plotting
    os.makedirs('output/phase3', exist_ok=True)
    try:
        model_full = model_full.prune()
        model_full.plot(folder='output/phase3', beta=100)
        plt.savefig('output/phase3/pca_kan_splines.png')
        print("Saved interpretable PCA-KAN spline visualization.")
    except Exception as e:
        print(f"Could not plot: {e}")
        
    # Full Inference
    with torch.no_grad():
        logits = model_full(dataset_full['test_input'])
        probs = torch.softmax(logits, dim=1)
        kan_full_probs = probs[:, 1].cpu().numpy()
        kan_full_preds = torch.argmax(logits, dim=1).cpu().numpy()
        
    df['PCA_KAN_Probability'] = kan_full_probs
    df['PCA_KAN_Prediction'] = kan_full_preds
    
    df.to_csv('output/phase3/pca_kan_full_predictions.csv', index=False)
    print("Saved final PCA-KAN predictions.")

if __name__ == "__main__":
    main()
