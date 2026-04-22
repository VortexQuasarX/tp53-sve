import os
import sys
import torch  # MUST BE FIRST
import numpy as np
import pandas as pd
from sklearn.preprocessing import StandardScaler
from kan import KAN

def run_counterfactual_rescue(target_mut='R175H'):
    # Load Data
    try:
        df = pd.read_csv("output/phase3/plddt_weighted_features.csv")
    except FileNotFoundError:
        return pd.DataFrame(), "Weighted features not found!"
        
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
    
    df['Label'] = df['Classification'].apply(lambda x: 1 if x == 'Likely Oncogenic' else (0 if x == 'Benign' else -1))
    labeled_df = df[df['Label'] != -1].copy()
    
    X = labeled_df[feature_cols].values
    y = labeled_df['Label'].values
    scaler = StandardScaler()
    X_s = scaler.fit_transform(X)
    device = torch.device('cuda' if torch.cuda.is_available() else 'cpu')
    
    # Train KAN Model
    model = KAN(width=[X.shape[1], 2, 2], grid=3, k=3, seed=42)
    dataset = {
        'train_input': torch.tensor(X_s, dtype=torch.float32).to(device),
        'train_label': torch.tensor(y, dtype=torch.long).to(device),
        'test_input': torch.tensor(X_s, dtype=torch.float32).to(device),
        'test_label': torch.tensor(y, dtype=torch.long).to(device)
    }
    model.fit(dataset, opt="LBFGS", steps=15, loss_fn=torch.nn.CrossEntropyLoss(), update_grid=False)
    
    target_row = df[df['Mutation'] == target_mut]
    if target_row.empty:
        return pd.DataFrame(), f"Target {target_mut} not found."
        
    x_orig_unscaled = target_row[feature_cols].values[0]
    x_orig = scaler.transform([x_orig_unscaled])[0]
    
    # Gradient Descent Optimization
    delta = torch.zeros(1, X.shape[1], requires_grad=True, dtype=torch.float32, device=device)
    x_tensor = torch.tensor([x_orig], dtype=torch.float32).to(device)
    optimizer = torch.optim.Adam([delta], lr=0.1)
    
    L1_LAMBDA = 0.01
    
    for epoch in range(1500):
        optimizer.zero_grad()
        x_new = x_tensor + delta
        logits = model(x_new)
        loss = -logits[0, 0] + logits[0, 1] + L1_LAMBDA * torch.sum(torch.abs(delta))
        loss.backward()
        optimizer.step()
        with torch.no_grad():
            delta.clamp_(-10.0, 10.0)
            
    with torch.no_grad():
        final_x = x_tensor + delta
        delta_unscaled = scaler.inverse_transform(final_x.cpu().numpy())[0] - x_orig_unscaled
        
    actionable_targets = []
    for i, col in enumerate(feature_cols):
        shift_scaled = delta[0, i].item()
        shift_raw = delta_unscaled[i]
        if abs(shift_scaled) > 0.1:
            actionable_targets.append({
                'Biological Feature': col,
                'Current Damaged Value': x_orig_unscaled[i],
                'Target Rescued Value': x_orig_unscaled[i] + shift_raw,
                'Required Delta': shift_raw,
                'Scaled Impact': shift_scaled
            })
            
    results_df = pd.DataFrame(actionable_targets)
    if not results_df.empty:
        results_df['Abs_Impact'] = results_df['Scaled Impact'].abs()
        results_df = results_df.sort_values(by='Abs_Impact', ascending=False).drop(columns=['Abs_Impact'])
        return results_df, "Success"
    return pd.DataFrame(), "No rescue path found."

def main():
    print("Run kan_counterfactual via the Streamlit app module now.")

if __name__ == "__main__":
    main()
