# Beyond Reconstruction: Virtual Cell Double OOD Evaluation Suite

Official codebase for:  
**Beyond Reconstruction: Mechanism- and Utility-Grounded Evaluation of Virtual Cells Under Double Out-of-Distribution Perturbation and Cellular-Context Shifts**  
Author: chessyuan  
Affiliation: the academic department, the institution, [location redacted]  
Correspondence: `[redacted]`  

---

## 1. Overview
This repository provides the complete, modular Python benchmarking suite to evaluate generative "virtual cell" models and single-cell foundation models under rigorous out-of-distribution (OOD) conditions. It implements:
- **13 Evaluated Models**: 5 statistical baselines (`BL0-BL4`), 2 machine-learning baselines (`BL5_RidgeChemGenomic`, `BL6_kNN_Tanimoto`), 5 deep architectures (`chemCPA`, `GEARS`, `CellOT`, `STATE`, `MAP`), and 1 foundation model head (`scGPT_Head`).
- **6 Evaluation Regimes**: In-distribution interpolation (ZS0), context holdout (ZS1_C), perturbation holdout (ZS1_P), double OOD (ZS2), scaffold holdout (`DRUG_OOD_Scaffold`), and mechanism-of-action holdout (`DRUG_OOD_MOA`).
- **Four-Tier Evaluative Hierarchy**:
  1. Whole-transcriptome Pearson correlation ($r_{\text{expr}}$).
  2. Perturbation delta prediction ($r_\Delta$) and Intervention Specificity Score (ISS).
  3. Mechanistic pathway fidelity (Hallmark, Reactome, KEGG) and generic stress decomposition (SVD $\text{PC1}_\Delta$).
  4. Downstream cancer oncology drug sensitivity screening utility (PRISM, GDSC v2, CTRP v2).
- **Statistical Rigor**: Two One-Sided Tests (TOST) equivalence testing against population mean baselines with hierarchical bootstrap resampling.
- **Reliability Boundary**: Quantitative support distance metric ($d_{\text{joint}} \le 0.75$) separating valid interpolation from ungrounded extrapolation.

---

## 2. Environment Setup

### Using Conda (Recommended)
```bash
conda env create -f environment.yml
conda activate virtual_cell_double_ood
```

### Using Pip
```bash
pip install -r requirements.txt
```

---

## 3. Data Acquisition & Layout
Place standardized pseudobulk and external cohort datasets under `data/`:
- `tahoe_fast_mode_pseudobulk.parquet` (from Tahoe-100M CELLxGENE release)
- `sciplex3_external_cohort.parquet` (GEO accession GSE139944)
- `op3_immune_cohort.parquet` (Zenodo OP3 immune perturbation screen)
- `proteomics_proteintalks.parquet` (ProteinTalks DIA-MS public release)
- `cancer_drug_efficacy_prism.parquet` (Broad Institute PRISM Repurposing 19Q4)
- `gdsc_v2_sensitivities.tsv` (Genomics of Drug Sensitivity in Cancer v2)
- `ctrp_v2_sensitivities.tsv` (Cancer Therapeutics Response Portal v2)

---

## 4. Reproducing Empirical Results

### Step 1: Run Full 13-Model Benchmark Across All 6 OOD Splits
```bash
python benchmark_runner.py --data_dir data/ --n_seeds 5 --output_dir results/
```

### Step 2: Run Statistical Equivalence Testing (TOST)
```bash
python task3_tost_equivalence.py --results_dir results/
```

### Step 3: Run Training-Only Stress SVD Decomposition
```bash
python task4_training_only_stress.py --data_dir data/
```

### Step 4: Run Multi-Database Pathway Fidelity Replication
```bash
python task5_multidatabase_pathway.py --data_dir data/
```

### Step 5: Run ProteinTalks Proteomics Matching Audit
```bash
python task6_proteintalks_audit.py --data_dir data/
```

### Step 6: Run PRISM Directionality Unit Test
```bash
python prism_direction_unit_test.py
```

### Step 7: Run Cross-Consortium Replication (GDSC v2 & CTRP v2)
```bash
python task8_gdsc_ctrp_replication.py
```

### Step 8: Validate the Joint Support Distance Threshold ($d_{\text{joint}} \le 0.75$)
```bash
python task9_djoint_validation.py
```

### Step 9: Generate All Main Manuscript Figures & Tables
```bash
python generate_figures.py --results_dir results/ --output_dir figures/
python generate_master_tables.py --results_dir results/ --output_dir tables/
```

---

## 5. Citation
```bibtex
@article{liu2026beyond,
  title={Beyond Reconstruction: Mechanism- and Utility-Grounded Evaluation of Virtual Cells Under Double Out-of-Distribution Perturbation and Cellular-Context Shifts},
  author={chessyuan},
  journal={Nature Methods (under review)},
  year={2026}
}
```
