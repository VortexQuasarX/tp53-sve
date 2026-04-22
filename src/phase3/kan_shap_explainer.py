import os
import sys
import torch  # MUST IMPORT TORCH FIRST TO PREVENT WINDOWS DLL CLASH
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from sklearn.preprocessing import StandardScaler
import shap

from tp53_sve import load_all_features
from kan import KAN

def main():
    print("====================================================")
    print(" EXPLAINABLE AI (XAI) - SHAP ON KAN OUTPUT")
    print("====================================================")
    
    # 1. Data Loading
    try:
        df = load_all_features()
    except Exception as e:
        print(f"Error loading features: {e}")
        return
        
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
    
    # Filter labeled data
    df['Label'] = df['Classification'].apply(
        lambda x: 1 if x == 'Likely Oncogenic' else (0 if x == 'Benign' else -1)
    )
    
    labeled_df = df[df['Label'] != -1].copy()
    
    X = labeled_df[feature_cols].values
    y = labeled_df['Label'].values
    
    print(f"Loaded {len(X)} labeled samples for SHAP analysis.")
    
    # Scale features
    scaler = StandardScaler()
    X_s = scaler.fit_transform(X)
    
    device = torch.device('cuda' if torch.cuda.is_available() else 'cpu')
    
    # 2. Train KAN Model
    print("Training KAN Model...")
    model = KAN(width=[X.shape[1], 2, 2], grid=3, k=3, seed=42)
    dataset = {
        'train_input': torch.tensor(X_s, dtype=torch.float32).to(device),
        'train_label': torch.tensor(y, dtype=torch.long).to(device),
        'test_input': torch.tensor(X_s, dtype=torch.float32).to(device),
        'test_label': torch.tensor(y, dtype=torch.long).to(device)
    }
    
    model.fit(dataset, opt="LBFGS", steps=20, loss_fn=torch.nn.CrossEntropyLoss(), update_grid=False)
    
    # Prune model for simpler inference
    model = model.prune()
    
    # 3. Create a wrapper function for SHAP (KernelExplainer expects numpy arrays)
    def kan_predict_proba(X_numpy):
        # SHAP provides numpy arrays, we need to return probabilities for the positive class
        X_tensor = torch.tensor(X_numpy, dtype=torch.float32).to(device)
        with torch.no_grad():
            logits = model(X_tensor)
            probs = torch.softmax(logits, dim=1)
            # Return probability of class 1 (Pathogenic)
            return probs[:, 1].cpu().numpy()
            
    # 4. Generate SHAP Values
    print("Computing SHAP values (KernelExplainer)...")
    try:
        # Use kmeans to summarize the background dataset if too large, but 67 is small enough
        background = shap.kmeans(X_s, 10) # 10 summary points
        explainer = shap.KernelExplainer(kan_predict_proba, background)
        
        # Calculate SHAP values for all labeled instances
        shap_values = explainer.shap_values(X_s, nsamples=100)
        
        # 5. Generate Visualizations
        os.makedirs("output/phase3", exist_ok=True)
        
        plt.figure(figsize=(10, 8))
        shap.summary_plot(shap_values, X_s, feature_names=feature_cols, show=False)
        plt.tight_layout()
        plt.savefig("output/phase3/kan_shap_summary.png", dpi=300, bbox_inches='tight')
        print("[SUCCESS] Saved SHAP Summary Plot to output/phase3/kan_shap_summary.png")
        plt.clf()
        
        # Feature Importance Bar Plot
        plt.figure(figsize=(10, 8))
        shap.summary_plot(shap_values, X_s, feature_names=feature_cols, plot_type="bar", show=False)
        plt.tight_layout()
        plt.savefig("output/phase3/kan_shap_bar.png", dpi=300, bbox_inches='tight')
        print("[SUCCESS] Saved SHAP Bar Plot to output/phase3/kan_shap_bar.png")
        
    except Exception as e:
        print(f"SHAP Error: {e}")

if __name__ == "__main__":
    main()
