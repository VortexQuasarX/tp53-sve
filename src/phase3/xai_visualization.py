import os
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
from sklearn.preprocessing import StandardScaler
from sklearn.discriminant_analysis import LinearDiscriminantAnalysis

def main():
    print("====================================================")
    print(" EXPLAINABLE AI (XAI) - FEATURE IMPORTANCE EXTRACTION")
    print("====================================================")
    
    # 1. Load the Data
    try:
        features_df = pd.read_csv("output/phase2/combined_features.csv")
    except Exception as e:
        print(f"Error loading data: {e}")
        return
        
    # Drop variants with missing class labels for training the XAI explainer
    train_df = features_df.dropna(subset=['Class'])
    
    # 2. Prepare Features
    # Exclude metadata columns
    exclude_cols = ['Mutation', 'Sequence', 'WT_Residue', 'Position', 'Mut_Residue', 'Class']
    feature_cols = [col for col in train_df.columns if col not in exclude_cols]
    
    X = train_df[feature_cols].values
    
    # Map Class: Pathogenic -> 1, Benign -> 0
    y = train_df['Class'].map({'Pathogenic': 1, 'Benign': 0}).values
    
    # 3. Standardize Features (Crucial for comparing LDA weights accurately)
    scaler = StandardScaler()
    X_scaled = scaler.fit_transform(X)
    
    # 4. Train the AI Model (Fisher's LDA)
    lda = LinearDiscriminantAnalysis()
    lda.fit(X_scaled, y)
    
    # 5. Extract the "Brain" of the AI (The Coefficients)
    # lda.coef_[0] gives the weight for each feature. 
    # Positive weight = pushes prediction towards Pathogenic (Class 1)
    # Negative weight = pushes prediction towards Benign (Class 0)
    coefficients = lda.coef_[0]
    
    # Create a DataFrame for the explainability
    importance_df = pd.DataFrame({
        'Feature': feature_cols,
        'Coefficient': coefficients,
        'Absolute_Importance': np.abs(coefficients)
    })
    
    # Sort by absolute importance to find the most critical features
    importance_df = importance_df.sort_values(by='Absolute_Importance', ascending=False)
    
    # 6. Console Output (The XAI Text Readout)
    print("\n[+] Top 15 Structural Features Driving the AI's Decision:\n")
    print(f"{'Feature Rank':<5} | {'Structural Metric Name':<35} | {'AI Weight (Directional Impact)'}")
    print("-" * 75)
    
    for i, row in enumerate(importance_df.head(15).itertuples(), 1):
        direction = "Pathogenic" if row.Coefficient > 0 else "Benign"
        print(f"#{i:<4} | {row.Feature:<35} | {row.Coefficient:>6.3f} (Drives towards {direction})")
        
    # 7. Generate Visualizable AI (The Graph)
    plt.figure(figsize=(12, 10))
    sns.set_theme(style="whitegrid")
    
    # Plot the top 15 features
    top_features = importance_df.head(15)
    
    # Color code by direction
    colors = ['#cc1100' if c > 0 else '#0066cc' for c in top_features['Coefficient']]
    
    ax = sns.barplot(
        x='Coefficient', 
        y='Feature', 
        data=top_features, 
        palette=colors
    )
    
    plt.title("Explainable AI (XAI): Internal Feature Weights for TP53 Pathogenicity Diagnosis", fontsize=16, pad=20, weight='bold')
    plt.xlabel("Directional Impact Weight (Standardized)\n<-- Drives toward Benign | Drives toward Pathogenic -->", fontsize=12, labelpad=15)
    plt.ylabel("Biophysical 3D Structural Features", fontsize=12)
    
    # Add a vertical line at 0
    plt.axvline(x=0, color='black', linestyle='-', linewidth=1.5)
    
    plt.tight_layout()
    
    os.makedirs("output/phase3", exist_ok=True)
    plt.savefig("output/phase3/xai_feature_importance.png", dpi=300, bbox_inches='tight')
    print("\n[SUCCESS] Saved Visualizable AI Graph to: output/phase3/xai_feature_importance.png")

if __name__ == "__main__":
    main()
