Write-Host 'Running Correlation'
python src\phase2\correlation_analysis.py
Write-Host 'Running SASA'
python src\phase3\sasa_analysis.py
Write-Host 'Running TM Score'
python src\phase3\tm_score.py
Write-Host 'Running Contact Network'
python src\phase3\contact_network.py
Write-Host 'Running ARES'
python src\phase3\tp53_ares.py
Write-Host 'Running DBCA'
python src\phase3\p53_dbca.py
Write-Host 'Running Local Global'
python src\phase3\local_global_impact.py
Write-Host 'Running Sec Struct'
python src\phase3\secondary_structure.py
Write-Host 'Running Compactness'
python src\phase3\compactness_torsion.py
Write-Host 'Running Composite'
python src\phase3\composite_scoring.py
Write-Host 'Running SPI'
python src\phase3\structural_pathogenicity_index.py
Write-Host 'Running SVE'
python src\phase3\tp53_sve.py
Write-Host 'Running Final Validation'
python src\phase3\final_validation.py
