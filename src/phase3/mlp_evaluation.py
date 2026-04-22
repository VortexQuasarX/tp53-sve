import os
import torch
import torch.nn as nn
import torch.optim as optim
import numpy as np
import pandas as pd
from sklearn.metrics import accuracy_score, roc_auc_score, matthews_corrcoef
from sklearn.model_selection import LeaveOneOut
from sklearn.preprocessing import StandardScaler
from sklearn.decomposition import PCA
from tp53_sve import load_all_features

# Standard MLP Architecture
class SimpleMLP(nn.Module):
    def __init__(self, input_dim):
        super(SimpleMLP, self).__init__()
        self.fc1 = nn.Linear(input_dim, 8)
        self.relu = nn.ReLU()
        self.fc2 = nn.Linear(8, 2)
        
    def forward(self, x):
        x = self.fc1(x)
        x = self.relu(x)
        x = self.fc2(x)
        return x

def main():
    print("====================================================")
    print(" BASELINE MLP: Standard Neural Network Comparison")
    print("====================================================")
    
    device = torch.device('cuda' if torch.cuda.is_available() else 'cpu')
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
    X_raw = labeled_df[feature_cols].values
    y = labeled_df['Label'].values
    X_full_raw = df[feature_cols].values
    
    print(f"Loaded {len(X_raw)} labeled samples.")
    
    # We will use the exact same 3 PCA components used for KAN to make it a perfectly fair comparison
    scaler = StandardScaler()
    X_full_scaled = scaler.fit_transform(X_full_raw)
    
    pca = PCA(n_components=3)
    X_full_pca = pca.fit_transform(X_full_scaled)
    X = X_full_pca[labeled_df.index]
    
    mlp = SimpleMLP(input_dim=3)
    params = sum(p.numel() for p in mlp.parameters() if p.requires_grad)
    print(f"MLP Parameters: {params} (Comparable to PCA-KAN's ~36)")
    
    loo = LeaveOneOut()
    mlp_preds = np.zeros(len(y))
    mlp_probs = np.zeros(len(y))
    train_accs = []
    
    for i, (train_index, test_index) in enumerate(loo.split(X)):
        X_train, X_test = X[train_index], X[test_index]
        y_train, y_test = y[train_index], y[test_index]
        
        # Scale
        fold_scaler = StandardScaler()
        X_train_s = fold_scaler.fit_transform(X_train)
        X_test_s = fold_scaler.transform(X_test)
        
        # Convert to tensors
        X_train_t = torch.tensor(X_train_s, dtype=torch.float32).to(device)
        y_train_t = torch.tensor(y_train, dtype=torch.long).to(device)
        X_test_t = torch.tensor(X_test_s, dtype=torch.float32).to(device)
        
        model = SimpleMLP(input_dim=3).to(device)
        criterion = nn.CrossEntropyLoss()
        optimizer = optim.Adam(model.parameters(), lr=0.01)
        
        # Train
        model.train()
        for epoch in range(100):  # 100 epochs is plenty for this tiny network
            optimizer.zero_grad()
            outputs = model(X_train_t)
            loss = criterion(outputs, y_train_t)
            loss.backward()
            optimizer.step()
            
        # Eval Train Acc
        model.eval()
        with torch.no_grad():
            train_preds = torch.argmax(model(X_train_t), dim=1)
            train_acc = (train_preds == y_train_t).float().mean().item()
            train_accs.append(train_acc)
            
            # Eval Test
            test_logits = model(X_test_t)
            test_probs = torch.softmax(test_logits, dim=1)
            mlp_preds[test_index] = torch.argmax(test_logits, dim=1).item()
            mlp_probs[test_index] = test_probs[0, 1].item()
            
    print("\n--- MLP RESULTS (LOOCV) ---")
    mlp_acc = accuracy_score(y, mlp_preds)
    mlp_auc = roc_auc_score(y, mlp_probs)
    mlp_mcc = matthews_corrcoef(y, mlp_preds)
    
    print(f"[MLP] Train Accuracy: {np.mean(train_accs):.3f}")
    print(f"[MLP] LOOCV Accuracy: {mlp_acc:.3f}")
    print(f"[MLP] LOOCV AUC-ROC:  {mlp_auc:.3f}")
    print(f"[MLP] LOOCV MCC:      {mlp_mcc:.3f}")
    
    # Train final and save
    model_full = SimpleMLP(input_dim=3).to(device)
    optimizer = optim.Adam(model_full.parameters(), lr=0.01)
    
    final_scaler = StandardScaler()
    X_s = final_scaler.fit_transform(X)
    X_full_s = final_scaler.transform(X_full_pca)
    
    X_s_t = torch.tensor(X_s, dtype=torch.float32).to(device)
    y_t = torch.tensor(y, dtype=torch.long).to(device)
    X_full_t = torch.tensor(X_full_s, dtype=torch.float32).to(device)
    
    model_full.train()
    for _ in range(100):
        optimizer.zero_grad()
        loss = criterion(model_full(X_s_t), y_t)
        loss.backward()
        optimizer.step()
        
    model_full.eval()
    with torch.no_grad():
        full_logits = model_full(X_full_t)
        full_probs = torch.softmax(full_logits, dim=1)[:, 1].cpu().numpy()
        full_preds = torch.argmax(full_logits, dim=1).cpu().numpy()
        
    df['MLP_Probability'] = full_probs
    df['MLP_Prediction'] = full_preds
    
    os.makedirs('output/phase3', exist_ok=True)
    df.to_csv('output/phase3/mlp_full_predictions.csv', index=False)
    print("\nSaved final MLP predictions to output/phase3/mlp_full_predictions.csv")

if __name__ == "__main__":
    main()
