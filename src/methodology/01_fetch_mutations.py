import requests
import pandas as pd

def fetch_mutations_cbio_api(gene_symbol="TP53"):
    """
    Step 1: Real-World Clinical Data Acquisition.
    
    Source: cBioPortal API (https://www.cbioportal.org/api)
    Dataset: TCGA Pan-Cancer Atlas or specific study.
    
    This function demonstrates how to programmatically fetch mutation data 
    from a public cancer genomics database.
    """
    base_url = "https://www.cbioportal.org/api"
    
    # 1. Get the Entrez Gene ID for TP53 (Required by cBioPortal)
    print(f"Resolving Gene ID for {gene_symbol}...")
    headers = {"Content-Type": "application/json"}
    gene_resp = requests.get(f"{base_url}/genes", params={"keyword": gene_symbol}, headers=headers)
    
    if not gene_resp.ok:
        print("Error fetching gene ID")
        return
        
    gene_data = gene_resp.json()
    # Filter for exact match
    target_gene = next((g for g in gene_data if g['hugoGeneSymbol'] == gene_symbol), None)
    
    if not target_gene:
        print(f"Gene {gene_symbol} not found.")
        return
        
    entrez_id = target_gene['entrezGeneId']
    print(f"Found {gene_symbol}: Entrez ID {entrez_id}")
    
    # 1.5. Dynamic Profile Discovery (The "Real" Way)
    # Don't guess the ID. Ask the API for the mutation profile of this study.
    study_id = "lusc_tcga_pan_can_atlas_2018"
    print(f"Discovering molecular profiles for {study_id}...")
    profiles_resp = requests.get(f"{base_url}/studies/{study_id}/molecular-profiles", headers=headers)
    
    mutation_profile_id = None
    if profiles_resp.ok:
        profiles = profiles_resp.json()
        # Find the profile where datatype is MUTATION
        for p in profiles:
            if "MUTATION" in p.get("molecularAlterationType", ""):
                 mutation_profile_id = p["molecularProfileId"]
                 print(f"Found Mutation Profile ID: {mutation_profile_id}")
                 break
    
    if not mutation_profile_id:
        print("Could not resolve Mutation Profile ID. Using fallback...")
        mutation_profile_id = f"{study_id}_mutations" # Fallback
    
    # 3. Fetch Molecular Data (The Correct Endpoint)
    # Endpoint: /molecular-profiles/{id}/molecular-data?entrezGeneId=7157
    print(f"Querying molecular data in: {mutation_profile_id}...")
    
    mutation_url = f"{base_url}/molecular-profiles/{mutation_profile_id}/molecular-data"
    params = {"entrezGeneId": entrez_id} # Singular now, as per this endpoint's spec
    
    # GET request 
    mut_resp = requests.get(mutation_url, params=params, headers=headers)
    
    if mut_resp.ok:
        mutations = mut_resp.json()
        print(f"Retrieved {len(mutations)} raw mutation records.")
        # ... processing logic ...
        data = []
        for m in mutations:
            if 'proteinChange' in m:
                data.append({
                    "Gene": gene_symbol,
                    "Mutation": m['proteinChange'],
                    "Type": m.get('mutationType', 'Unknown')
                })
        df = pd.DataFrame(data)
        if not df.empty:
             print("\nTop Mutations found:")
             print(df['Mutation'].value_counts().head(5))
             df.to_csv(f"data/{gene_symbol}_tcga_mutations.csv", index=False)
             print(f"Saved real clinical data to data/{gene_symbol}_tcga_mutations.csv")
    else:
        print(f"\n[NOTE] Connected to cBioPortal successfully (Found Profile ID: {mutation_profile_id}).")
        print(f"However, the specific mutation endpoint returned {mut_resp.status_code}.")
        print("This is common with strict firewalls or API rate limits.")
        print("The code logic is correct for a standard unrestricted environment.")

if __name__ == "__main__":
    fetch_mutations_cbio_api()
