from Bio import Entrez, SeqIO

def fetch_sequence_biopython(accession_id="NP_000537.3"):
    """
    Step 2: Real-World Sequence Acquisition.
    
    Library: Biopython (Bio.Entrez)
    Source: NCBI GenBank
    Target: NM_000546.6 (TP53 Transcript) -> NP_000537.3 (Protein Isoform a)
    
    This demonstrates the standard bioinformatics method to fetch 
    official sequences.
    """
    # 1. NCBI requires an email for API usage
    Entrez.email = "your.email@example.com"
    
    print(f"Querying NCBI Entrez for {accession_id}...")
    
    try:
        # 2. Fetch the record using efetch
        handle = Entrez.efetch(db="protein", id=accession_id, rettype="gb", retmode="text")
        record = SeqIO.read(handle, "genbank")
        handle.close()
        
        print(f"Successfully fetched: {record.description}")
        print(f"Sequence Length: {len(record.seq)} amino acids")
        
        # 3. Save as FASTA (Standard Format)
        output_file = "data/tp53_biopython_wt.fasta"
        SeqIO.write(record, output_file, "fasta")
        print(f"Saved FASTA to {output_file}")
        
    except Exception as e:
        print(f"Error connecting to NCBI: {e}")

if __name__ == "__main__":
    fetch_sequence_biopython()
