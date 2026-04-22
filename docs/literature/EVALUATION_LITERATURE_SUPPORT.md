# TP53 Mutation Analysis: Evaluation Metrics & Literature Support

This document maps every single evaluation we have implemented in this project to the foundational scientific papers and algorithms that justify its use. This is crucial for verifying that our methods are scientifically sound.

---

## 1. Structural Origin Data

### AlphaFold Coordinate Generation
* **What We Used:** The 3D atom coordinates (x, y, z) and confidence metrics (pLDDT) for our 50 mutant structures.
* **Literature Support:**
  > **Jumper, J., et al. (2021). "Highly accurate protein structure prediction with AlphaFold." *Nature*, 596(7873), 583-589.**
  * *Justification:* Establishes AlphaFold as the gold standard for predicting protein structures from sequence, accurate enough to rival experimental crystal structures.

---

## 2. Geometric & Physical Alignments (Phase 1 & Phase 2)

### RMSD (Root Mean Square Deviation) & Superposition
* **What We Used:** The base metric for structural difference and the Kabsch algorithm for aligning Cα atoms.
* **Literature Support:**
  > **Kabsch, W. (1976). "A solution for the best rotation to relate two sets of vectors." *Acta Crystallographica Section A*, 32(5), 922-923.**
  * *Justification:* The mathematical proof that singular value decomposition (SVD) optimally minimizes distances between two identical sets of points.

### SASA (Solvent Accessible Surface Area)
* **What We Used:** Measuring hydrophobic exposure to water molecules (using a 1.4 Å radius probe).
* **Literature Support:**
  > **Shrake, A., & Rupley, J. A. (1973). "Environment and exposure to solvent of protein atoms." *Journal of Molecular Biology*, 79(2), 351-371.**
  * *Justification:* The foundational algorithm (rolling a sphere over the protein surface) to calculate which atoms are exposed to water.

### Secondary Structure Assignments
* **What We Used:** Cα-Cα geometric distances (i to i+2, i to i+3, i to i+4) to determine a-helices, b-sheets, and coils.
* **Literature Support:**
  > **Kabsch, W., & Sander, C. (1983). "Dictionary of protein secondary structure: pattern recognition of hydrogen-bonded and geometrical features." *Biopolymers*, 22(12), 2577-2637.**
  * *Justification:* The basis for DSSP (Define Secondary Structure of Proteins) and geometric constraints for structural elements.

---

## 3. Advanced Structural Metrics (Phase 3 - Tier 2)

### TM-Score (Template Modeling Score)
* **What We Used:** A length-independent structural similarity metric that resists outlier penalties to evaluate if the protein "fold" was destroyed.
* **Literature Support:**
  > **Zhang, Y., & Skolnick, J. (2004). "Scoring function for automated assessment of protein structure template quality." *Proteins*, 57(4), 702-710.**
  * *Justification:* Proves TM-score is superior to RMSD for identifying structural similarity and establishes the >0.5 threshold for "same fold".

### p53-DBCA (DNA-Binding Competence Assessment)
* **What We Used:** Probes specifically measuring zinc coordination (C176, H179, C238, C242), DNA contact points (R248, R273, K120, etc.), and loop integrity.
* **Literature Support:**
  > **Cho, Y., et al. (1994). "Crystal structure of a p53 tumor suppressor-DNA complex." *Science*, 265(5170), 346-355.**
  * *Justification:* Maps exactly which residues are critical for interacting with DNA and coordinating the crucial zinc ion.
  > **Bullock, A. N., et al. (1997). "Thermodynamic stability of wild-type and mutant p53 core domain." *PNAS*, 94(26), 14338-14342.**
  * *Justification:* Identifies how hotspot mutations specifically compromise the structural scaffold (hydrophobic core & H-bonds) and zinc interactions.

---

## 4. Energy & Thermodynamic Scoring (Phase 3 - Tier 3)

### TP53-ARES (Atomistic Residue Energy Scoring)
* **What We Used:** Contact network disruption propagation and empirical energy potentials (measuring the ΔΔE of broken vs. gained contacts).
* **Literature Support:**
  > **Miyazawa, S., & Jernigan, R. L. (1996). "Residue-residue potentials with a favorable contact pair term and an unfavorable high packing density term, for simulation and folding." *Journal of Molecular Biology*, 256(3), 623-644.**
  * *Justification:* Provides the heavily vetted 20x20 statistical contact potential matrix (empirical contact energies) we used to calculate whether a structural rewiring is stabilizing or destabilizing.

---

## 5. Machine Learning & Evolutionary Classifier (Phase 3 - Tier 3)

### TP53-SVE (Structural Variant Evaluator) 
* **What We Used:** Evolutionary sequence matrices and Fisher's Linear Discriminant Analysis to optimally separate pathogenic from benign variants.
* **Literature Support:**
  > **Henikoff, S., & Henikoff, J. G. (1992). "Amino acid substitution matrices from protein blocks." *PNAS*, 89(22), 10915-10919.**
  * *Justification:* The source of the BLOSUM62 matrix used in SVE, which empirically scores the evolutionary likelihood/damage of replacing one amino acid with another.
  > **Fisher, R. A. (1936). "The use of multiple measurements in taxonomic problems." *Annals of Eugenics*, 7(2), 179-188.**
  * *Justification:* The mathematical founding of Fisher's Linear Discriminant, which mathematically guarantees the optimal weight vector to maximize the variance between pathogenic and benign classes in our 34-feature SVE.

---

## 6. Established Baseline Comparisons (Phase 2)

### SIFT & PolyPhen-2
* **What We Used:** For comparing our structure-based rankings against established sequence/evolution tools.
* **Literature Support:**
  > **Ng, P. C., & Henikoff, S. (2003). "SIFT: Predicting amino acid changes that affect protein function." *Nucleic Acids Research*, 31(13), 3812-3814.**
  > **Adzhubei, I. A., et al. (2010). "A method and server for predicting damaging missense mutations." *Nature Methods*, 7(4), 248-249.**
  * *Justification:* The origin of the SIFT and PolyPhen-2 algorithms respectively, proving why standard tools give binary limits while our continuous structural mapping is superior.
