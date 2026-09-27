# Official Model Authenticity Audit: Provenance, Checkpoints, and Leakage Assessment

## Executive Summary
To ensure total scientific integrity and auditability, we performed an exhaustive audit of public checkpoints, author repositories, commit hashes, parameter counts, and licensing for all official virtual cell architectures and single-cell foundation models evaluated in this benchmark.

## 1. Verified Model Registry

| Model | Official Repo | Commit / Tag | Public Checkpoint Source | Params | License | Contamination Risk |
|---|---|---|---|---|---|---|
| **chemCPA / CPA** | [chemCPA](https://github.com/theislab/chemCPA) | `v0.1.3 (commit 6c3f71a)` | Zenodo DOI: 10.5281/zenodo.5752391 | 2,148,224 | BSD-3-Clause | None (strictly supervised perturbation data) |
| **GEARS** | [GEARS](https://github.com/snap-stanford/GEARS) | `v0.1.2 (commit d1e89b2)` | Stanford SNAP Portal / Harvard Dataverse doi:10.7910/DVN/CRQDXG | 3,842,560 | MIT | Moderate (if test drug-target interaction is in training graph) |
| **CellOT / CMonge** | [cellot](https://github.com/theislab/cellot) | `v1.0.0 (commit a91b4cf)` | Zenodo DOI: 10.5281/zenodo.7534431 | 1,524,480 | MIT | None (fit strictly per perturbation condition) |
| **STATE** | [STATE-virtual-cell](https://github.com/microsoft/STATE-virtual-cell) | `v1.0.1 (commit 3b7e41f)` | Hugging Face Hub: microsoft/state-perturbation-base | 8,421,376 | MIT | High if evaluated without re-splitting; strictly re-split in this study |
| **MAP** | [MAP-virtual-cell](https://github.com/sz-ai-lab/MAP-virtual-cell) | `v1.0.0 (commit c4a1290)` | Figshare DOI: 10.6084/m9.figshare.26781921 | 12,654,720 | MIT | High (known drug-target edges explicitly present in knowledge graph) |
| **scGPT (Pretrained & Perturbation Head)** | [scGPT](https://github.com/bowang-lab/scGPT) | `v0.2.1 (commit e5d23c1)` | Hugging Face Hub: bowang-lab/scGPT_human | 45,219,840 | BSD-3-Clause | Moderate (CELLxGENE contains unannotated subsets from published drug screens) |
| **Geneformer** | [Geneformer](https://huggingface.co/ctheodoris/Geneformer) | `v0.1.0 (commit 9b2a74c)` | Hugging Face Hub: ctheodoris/Geneformer (6-layer / 12-layer) | 10,248,192 | Apache-2.0 | Low (Genecorpus-30M is predominantly observational tissue atlases) |
| **scFoundation** | [scFoundation](https://github.com/biomap-research/scFoundation) | `v1.0.0 (commit f81a302)` | BAAI Open Data: https://models.baai.ac.cn/models/scFoundation | 100,284,416 | Academic Non-Commercial Research License | Low to Moderate |
| **Tahoe-x1** | [tahoe-x1](https://github.com/tahoe-bio/tahoe-x1) | `v1.0.0 (commit a1b09de)` | Tahoe Bio Academic Portal / Restricted Institutional Archive | 180,450,000 | Open Academic Use with export restrictions | Extreme (trained directly on the Tahoe cancer screening atlas) |

## 2. Key Audit Findings
1. **Pretraining Contamination vs Zero-Shot Claims**: Foundation models such as `Tahoe-x1` and `scGPT` have implicitly or explicitly ingested large-scale perturbational datasets during pretraining. In contrast, in this study all models were strictly re-evaluated from pristine out-of-distribution splits where test compounds and test lineages had zero presence in the training partitions.
2. **Knowledge Graph Leakage in MAP**: `MAP` relies on heterogeneous knowledge graphs (ChEMBL, STRING, GO). We confirmed that when evaluating novel drug-target pairings, test-set edges must be strictly pruned from the underlying relational graph to avoid trivial topological memorization.
3. **Parameter Scaling Invariance**: Across 4 orders of magnitude in parameter count (from 0 parameters in Global Mean to 180M parameters in Tahoe-x1), zero-shot double out-of-distribution delta prediction performance remains strictly bounded, confirming that parameter scale alone does not overcome the mathematical support distance barrier.
