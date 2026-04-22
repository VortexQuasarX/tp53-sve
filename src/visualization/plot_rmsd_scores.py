import pandas as pd
import matplotlib.pyplot as plt
import os

def main():
    csv_file = "output/rmsd_scores.csv"
    output_file = "output/rmsd_ranking.png"
    
    if not os.path.exists(csv_file):
        print("RMSD Scores file not found. Run calculation first.")
        return

    df = pd.read_csv(csv_file)
    
    # Sort by RMSD to make the chart readable
    df = df.sort_values(by="RMSD (Angstroms)", ascending=True)
    
    # Plot
    plt.figure(figsize=(12, 8))
    
    # Create horizontal bar chart
    bars = plt.barh(df['Mutation'], df['RMSD (Angstroms)'], color='salmon')
    
    # Highlight the "Top Destructive" ones in Red
    # (Just a visual flair, making top 3 dark red)
    for i in range(len(bars)-3, len(bars)):
        bars[i].set_color('darkred')

    plt.xlabel('RMSD (Angstroms) - Higher means more structural devastation', fontsize=12)
    plt.title('Structural Impact of Each Mutation (RMSD)', fontsize=16)
    plt.grid(axis='x', linestyle='--', alpha=0.7)
    
    # Add values to the end of bars
    for index, value in enumerate(df['RMSD (Angstroms)']):
        plt.text(value + 0.01, index, str(value), va='center')

    plt.tight_layout()
    plt.savefig(output_file, dpi=150)
    print(f"Plot saved to {output_file}")

if __name__ == "__main__":
    main()
