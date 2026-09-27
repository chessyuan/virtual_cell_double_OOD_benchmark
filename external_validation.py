import os, json
import numpy as np
import pandas as pd
from sklearn.linear_model import Ridge

from chemical_genomics import smiles_to_ecfp4
from baselines import Baseline1_GlobalMean
from metrics import compute_level2_delta_metrics_vectorized, compute_level3_mechanisms_vectorized

base_dir = 'project_05_virtual_cell_double_OOD'
genes_df = pd.read_csv(os.path.join(base_dir, '06_metadata/genes_universe.tsv'), sep='	')
gene_cols = genes_df['gene_symbol'].tolist()
n_genes = len(gene_cols)

drug_df = pd.read_csv(os.path.join(base_dir, '06_metadata/drug_metadata.tsv'), sep='	')
drug_fps = {row['drug_name']: smiles_to_ecfp4(row['canonical_smiles']) for _, row in drug_df.iterrows()}

with open(os.path.join(base_dir, '06_metadata/hallmark_pathways.json')) as f:
    hallmark_pathways = json.load(f)

atlas_df = pd.read_parquet(os.path.join(base_dir, '05_processed_data/tahoe_fast_mode_pseudobulk.parquet'))
ctrl_df = atlas_df[atlas_df['is_control'] == True]
basal_dict = {cell: group[gene_cols].values[0] for cell, group in ctrl_df.groupby('cell_line')}

# Train on Tahoe solid cancer cell lines
train_cancer = atlas_df[(atlas_df['is_control'] == False) & (~atlas_df['cell_line'].isin(['CD4_T_cells', 'Monocytes']))]
train_deltas = np.array([r[gene_cols].values.astype(np.float32) - basal_dict[r['cell_line']] for _, r in train_cancer.iterrows()])
train_basals = np.array([basal_dict[r['cell_line']] for _, r in train_cancer.iterrows()])
train_fps = np.array([drug_fps.get(r['perturbation'], np.zeros(1024)) for _, r in train_cancer.iterrows()])
train_doses = np.array(train_cancer['dose_uM'].values, dtype=np.float32)
train_X = np.concatenate([train_basals, train_fps, train_doses[:, None]], axis=1)

u, s, vt = np.linalg.svd(train_deltas - np.mean(train_deltas, axis=0), full_matrices=False)
generic_pc1 = vt[0]

bl1 = Baseline1_GlobalMean(); bl1.fit(train_deltas)
ridge = Ridge(alpha=25.0); ridge.fit(train_X, train_deltas)

ext_results = []

# 1. External Cohort 1: SciPlex3
sciplex_df = pd.read_parquet(os.path.join(base_dir, '05_processed_data/sciplex3_external_cohort.parquet'))
sciplex_treated = sciplex_df[sciplex_df['is_control'] == False]
s_deltas = np.array([r[gene_cols].values.astype(np.float32) - basal_dict[r['cell_line']] for _, r in sciplex_treated.iterrows()])
s_basals = np.array([basal_dict[r['cell_line']] for _, r in sciplex_treated.iterrows()])
s_fps = np.array([drug_fps.get(r['perturbation'], np.zeros(1024)) for _, r in sciplex_treated.iterrows()])
s_doses = np.array(sciplex_treated['dose_uM'].values, dtype=np.float32)
s_X = np.concatenate([s_basals, s_fps, s_doses[:, None]], axis=1)

pred_bl1_s = np.tile(bl1.mean_delta, (len(sciplex_treated), 1))
pred_ridge_s = ridge.predict(s_X)

l2_bl1 = compute_level2_delta_metrics_vectorized(pred_bl1_s, s_deltas)
l3_bl1 = compute_level3_mechanisms_vectorized(pred_bl1_s, s_deltas, gene_cols, hallmark_pathways, generic_pc1)
l2_rg = compute_level2_delta_metrics_vectorized(pred_ridge_s, s_deltas)
l3_rg = compute_level3_mechanisms_vectorized(pred_ridge_s, s_deltas, gene_cols, hallmark_pathways, generic_pc1)

ext_results.append({'cohort': 'SciPlex3 (Cancer)', 'model': 'GlobalMean', 'pcc_delta': l2_bl1['pcc_delta'], 'pcc_residual_delta': l3_bl1['pcc_residual_delta'], 'pathway_fidelity': l3_bl1['pathway_fidelity']})
ext_results.append({'cohort': 'SciPlex3 (Cancer)', 'model': 'Ridge', 'pcc_delta': l2_rg['pcc_delta'], 'pcc_residual_delta': l3_rg['pcc_residual_delta'], 'pathway_fidelity': l3_rg['pathway_fidelity']})

# 2. External Cohort 2: OP3 (Immune Cells - Extreme Context Shift)
op3_df = pd.read_parquet(os.path.join(base_dir, '05_processed_data/op3_immune_cohort.parquet'))
op3_treated = op3_df[op3_df['is_control'] == False]
op_deltas = np.array([r[gene_cols].values.astype(np.float32) - basal_dict[r['cell_line']] for _, r in op3_treated.iterrows()])
op_basals = np.array([basal_dict[r['cell_line']] for _, r in op3_treated.iterrows()])
op_fps = np.array([drug_fps.get(r['perturbation'], np.zeros(1024)) for _, r in op3_treated.iterrows()])
op_doses = np.array(op3_treated['dose_uM'].values, dtype=np.float32)
op_X = np.concatenate([op_basals, op_fps, op_doses[:, None]], axis=1)

pred_bl1_op = np.tile(bl1.mean_delta, (len(op3_treated), 1))
pred_ridge_op = ridge.predict(op_X)

l2_bl1_op = compute_level2_delta_metrics_vectorized(pred_bl1_op, op_deltas)
l3_bl1_op = compute_level3_mechanisms_vectorized(pred_bl1_op, op_deltas, gene_cols, hallmark_pathways, generic_pc1)
l2_rg_op = compute_level2_delta_metrics_vectorized(pred_ridge_op, op_deltas)
l3_rg_op = compute_level3_mechanisms_vectorized(pred_ridge_op, op_deltas, gene_cols, hallmark_pathways, generic_pc1)

ext_results.append({'cohort': 'OP3 (Primary Immune)', 'model': 'GlobalMean', 'pcc_delta': l2_bl1_op['pcc_delta'], 'pcc_residual_delta': l3_bl1_op['pcc_residual_delta'], 'pathway_fidelity': l3_bl1_op['pathway_fidelity']})
ext_results.append({'cohort': 'OP3 (Primary Immune)', 'model': 'Ridge', 'pcc_delta': l2_rg_op['pcc_delta'], 'pcc_residual_delta': l3_rg_op['pcc_residual_delta'], 'pathway_fidelity': l3_rg_op['pathway_fidelity']})

# 3. External Cohort 3: Replogle Genetic Knockdowns
gen_df = pd.read_parquet(os.path.join(base_dir, '05_processed_data/genetic_perturbation_cohort.parquet'))
gen_treated = gen_df[gen_df['is_control'] == False]
k562_ctrl = basal_dict['K562']
g_deltas = np.array([r[gene_cols].values.astype(np.float32) - k562_ctrl for _, r in gen_treated.iterrows()])

pred_bl1_g = np.tile(bl1.mean_delta, (len(gen_treated), 1))
l2_bl1_g = compute_level2_delta_metrics_vectorized(pred_bl1_g, g_deltas)
l3_bl1_g = compute_level3_mechanisms_vectorized(pred_bl1_g, g_deltas, gene_cols, hallmark_pathways, generic_pc1)

ext_results.append({'cohort': 'Replogle (Genetic CRISPRi)', 'model': 'GlobalMean', 'pcc_delta': l2_bl1_g['pcc_delta'], 'pcc_residual_delta': l3_bl1_g['pcc_residual_delta'], 'pathway_fidelity': l3_bl1_g['pathway_fidelity']})

df_ext = pd.DataFrame(ext_results)
ext_out = os.path.join(base_dir, '20_cross_dataset/external_cohort_validation.tsv')
df_ext.to_csv(ext_out, sep='	', index=False)
print(df_ext.round(4))
