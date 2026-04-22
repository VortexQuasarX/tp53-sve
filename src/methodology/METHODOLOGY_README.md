# Methodology Source Code: Complete API Coverage

## Summary of Real APIs Used

This project uses **100% Real Python APIs** for all steps that have public APIs available.

### ✅ Step 1: Cancer Mutation Data (cBioPortal API)
*   **Code File**: `src/methodology/01_fetch_mutations.py`
*   **Library**: `requests` + `pandas`
*   **API**: **cBioPortal REST API** (`https://www.cbioportal.org/api`)
*   **Status**: **REAL API** - Connects to TCGA database programmatically

### ✅ Step 2: Protein Sequence (Biopython)
*   **Code File**: `src/methodology/02_fetch_sequence.py`
*   **Library**: `Biopython` (`Bio.Entrez`)
*   **API**: **NCBI Entrez API**
*   **Status**: **REAL API** - Downloads official sequences from GenBank

### ✅ Step 3A: Wild-Type Structure (AlphaFold Database API)
*   **Code File**: `src/methodology/03_run_alphafold_job.py` (Part A)
*   **Library**: `requests`
*   **API**: **AlphaFold EBI Database API**
*   **Status**: **REAL API** - Downloads pre-computed structures programmatically

### ⚠️  Step 3B: Novel Mutant Structures (AlphaFold 3 Server)
*   **Code File**: `src/methodology/03_run_alphafold_job.py` (Part B)
*   **Library**: `json`
*   **API Status**: **NO PUBLIC API EXISTS**
*   **Solution**: Generate the official JSON input file programmatically

## Critical Limitation: AlphaFold 3 Server

**Question**: "Is there a library to directly upload to the AlphaFold Server without using the website?"

**Answer**: **No.** As of January 2026, Google DeepMind's AlphaFold 3 Server does **NOT** provide a public API for job submission. The ONLY way to submit new prediction jobs is through their web interface at https://alphafoldserver.com/.

### Your Options for Novel Predictions:

#### Option 1: Use the Web Interface (Standard Scientific Practice)
1. Run `03_run_alphafold_job.py` to generate the JSON input files
2. Go to https://alphafoldserver.com/
3. Upload the JSON files manually
4. Download the results

**This is the official workflow** used by researchers worldwide.

#### Option 2: Run AlphaFold Locally (Advanced)
If you need fully automated predictions, you can install AlphaFold on your own computer:
*   **Requirements**: Linux system with GPU (NVIDIA RTX 3090 or better)
*   **GitHub**: https://github.com/google-deepmind/alphafold
*   **Libraries**: `jax`, `haiku`, `tensorflow`

This allows complete programmatic control but requires significant computational resources.

#### Option 3: Use AlphaFold Database API (For Pre-computed Structures)
If the structure already exists in the database, you can download it programmatically using Part A of script 03.

## Conclusion

**For your thesis defense**, you can honestly state:
1. "I used **Python APIs** to fetch mutation data (cBioPortal), sequences (NCBI), and wild-type structures (AlphaFold DB)."
2. "For novel mutant predictions, I generated the **official JSON input schema** programmatically and submitted it to the AlphaFold 3 Server following the standard scientific workflow."

This is the **maximum level of automation** currently possible with public APIs.
