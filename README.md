# The "MmpL Collapse": Genomic Entropy and Evolutionary Game Theory Reveal the Mechanism of Mycobacterial Drug Resistance

![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg) ![Status: Validated](https://img.shields.io/badge/Status-Validated-green.svg) ![Language: Python](https://img.shields.io/badge/Language-Python%203.10-blue.svg)

## 🧬 Abstract

The irreversible transition from the drug-sensitive "Smooth" (S) to the hyper-virulent, drug-resistant "Rough" (R) morphotype in *Mycobacterium abscessus* is the primary driver of treatment failure in cystic fibrosis patients. While historically viewed as a stochastic mutational event, the biophysical mechanism governing this switch has remained undefined.

In this study, we present a unified computational framework that decodes the S-to-R transition as a deterministic **thermodynamic phase shift**. By integrating **Genomic Entropy Scanning** with **Topological Data Analysis (TDA)** across 260 clinical isolates, we identify the causal mechanism: **"Quantum Slippage,"** a hyper-unstable Poly-G tract within the *mmpL* lipid transporter gene. This locus exhibits bimodal structural collapse ($Entropy \approx 1.15$), oscillating between a functional 18bp wild-type state and a collapsed 4bp frameshift variant.

Topological modeling reveals that this genomic collapse triggers a macroscopic structural reorganization, where the R-variant adopts a rigid **"Fortress" topology** characterized by high-persistence $H_1$ homology loops. Furthermore, **Evolutionary Game Theory** modeling confirms that this switch is thermodynamically inevitable: the high metabolic cost of synthesizing Smooth surface lipids ($c = 1.7$ ATP/unit) creates an evolutionary trap under stress, forcing the population to defect to the energy-efficient Rough state.

---

## 📂 Repository Structure

```text
M_abscessus_Rough_Switch/
│
├── 📜 README.md
│
├── 📁 01_Genomic_Evidence          <-- PHASE 1: The Mutation
│   ├── Table_1_Genomic_Validation.csv      # BLAST-validated MmpL/TetR lesions
│   ├── Fig_1_Genomic_PCA_Clusters.png      # PCA projection of entropy (S vs R)
│   ├── Fig_2_Comparative_Persistence.png   # Validation: M. abscessus Loop vs M. tb Tree
│   ├── MW2G8DCH016-Alignment.json          # Raw alignment data source
│   ├── 01_Genotype_Mechanisms.ipynb        # Source Notebook
│   └── 01_genotype_mechanisms.py           # Python Script
│
├── 📁 02_Biophysics_Topology       <-- PHASE 2: The Structure
│   ├── Fig_3_Topological_Barcode.png       # Persistent Homology (The "Fortress")
│   ├── Fig_4_Evolutionary_Collapse.png     # Dynamics of the structural failure
│   ├── Fig_5_Evolutionary_Phase_Space.png  # Heatmap of metabolic cost ($c=1.7$)
│   ├── Fig_6_Bifurcation_TippingPoint.png  # The S-Curve switch threshold
│   ├── 02_Phenotype_Topology.ipynb         # Source Notebook
│   └── 02_phenotype_topology.py            # Python Script
│
└── 📁 03_Clinical_Synthesis        <-- PHASE 3: The Application
    ├── Fig_7_Clinical_Risk_Dashboard.png   # Probability model for patient risk
    ├── 03_Evolutionary_Synthesis.ipynb     # Source Notebook
    └── 03_evolutionary_synthesis.py        # Python Script
