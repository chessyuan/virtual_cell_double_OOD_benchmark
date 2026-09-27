import os, json
import numpy as np
import pandas as pd
from scipy import stats
from sklearn.metrics import roc_auc_score, precision_recall_curve, auc
from sklearn.linear_model import Ridge

from chemical_genomics import smiles_to_ecfp4, compute_tanimoto_similarity
from models import STATE_Module, ChemCPA_Module
from metrics import compute_level2_delta_metrics_vectorized, compute_level3_mechanisms_vectorized

base_dir = 'project_05_virtual_cell_double_OOD'
atlas_df = pd.read_parquet(os.path.join(base_dir, '05_processed_data/tahoe_fast_mode_pseudobulk.parquet'))
genes_df = pd.read_csv(os.path.join(base_dir, '06_metadata/genes_universe.tsv'), sep='	')
drug_df = pd.read_csv(os.path.join(base_dir, '06_metadata/drug_metadata.tsv'), sep='	')
cell_df = pd.read_csv(os.path.join(base_dir, '06_metadata/cell_line_metadata.tsv'), sep='	')
prism_df = pd.read_parquet(os.path.join(base_dir, '05_processed_data/cancer_drug_efficacy_prism.parquet'))
prot_df = pd.read_parquet(os.path.join(base_dir, '05_processed_data/proteomics_proteintalks.parquet'))
with open(os.path.join(base_dir, '06_metadata/hallmark_pathways.json')) as f:
    hallmark_pathways = json.load(f)

gene_cols = genes_df['gene_symbol'].tolist()
n_genes = len(gene_cols)
gene2idx = {g: i for i, g in enumerate(gene_cols)}

drug_fps = {row['drug_name']: smiles_to_ecfp4(row['canonical_smiles']) for _, row in drug_df.iterrows()}

ctrl_df = atlas_df[atlas_df['is_control'] == True]
basal_dict = {cell: group[gene_cols].values[0] for cell, group in ctrl_df.groupby('cell_line')}

# -------------------------------------------------------------
# 1. UNCERTAINTY & SELECTIVE ABSTENTION (COVERAGE-RISK CURVES)
# -------------------------------------------------------------
print('--> Step 16: Uncertainty Quantification & Selective Abstention', flush=True)
treated_df = atlas_df[atlas_df['is_control'] == False]
all_deltas = np.array([row[gene_cols].values.astype(np.float32) - basal_dict[row['cell_line']] for _, row in treated_df.iterrows()])
u, s, vt = np.linalg.svd(all_deltas - np.mean(all_deltas, axis=0), full_matrices=False)
generic_pc1 = vt[0]

# Train STATE on seen drugs and seen lines, test on Double OOD (seed 42)
with open(os.path.join(base_dir, '07_splits/split_protocol_seed_42.json')) as f:
    split_info = json.load(f)
test_drugs = split_info['ZS2_Double_OOD']['test_drugs']
test_cells = split_info['ZS2_Double_OOD']['test_cells']

train_df = treated_df[(~treated_df['perturbation'].isin(test_drugs)) & (~treated_df['cell_line'].isin(test_cells))]
test_df = treated_df[(treated_df['perturbation'].isin(test_drugs)) & (treated_df['cell_line'].isin(test_cells))]

train_deltas = np.array([r[gene_cols].values.astype(np.float32) - basal_dict[r['cell_line']] for _, r in train_df.iterrows()])
test_deltas = np.array([r[gene_cols].values.astype(np.float32) - basal_dict[r['cell_line']] for _, r in test_df.iterrows()])

train_fps = np.array([drug_fps[r['perturbation']] for _, r in train_df.iterrows()])
test_fps = np.array([drug_fps[r['perturbation']] for _, r in test_df.iterrows()])

train_basals = np.array([basal_dict[r['cell_line']] for _, r in train_df.iterrows()])
test_basals = np.array([basal_dict[r['cell_line']] for _, r in test_df.iterrows()])

train_doses = np.array(train_df['dose_uM'].values, dtype=np.float32)
test_doses = np.array(test_df['dose_uM'].values, dtype=np.float32)

import torch
import torch.optim as optim
torch.manual_seed(42)
state_model = STATE_Module()
opt = optim.AdamW(state_model.parameters(), lr=1e-3)
crit = torch.nn.MSELoss()
t_b = torch.tensor(train_basals, dtype=torch.float32)
t_c = torch.tensor(train_fps, dtype=torch.float32)
t_d = torch.tensor(train_doses, dtype=torch.float32)
t_y = torch.tensor(train_deltas, dtype=torch.float32)

for ep in range(30):
    opt.zero_grad()
    pred, mu, logvar = state_model(t_b, t_c, t_d)
    kl = -0.5 * torch.sum(1 + logvar - mu.pow(2) - logvar.exp()) / t_b.shape[0]
    loss = crit(pred, t_y) + 0.01 * kl
    loss.backward()
    opt.step()

# MC dropout sampling (20 passes)
test_b = torch.tensor(test_basals, dtype=torch.float32)
test_c = torch.tensor(test_fps, dtype=torch.float32)
test_d = torch.tensor(test_doses, dtype=torch.float32)

state_model.eval()
mc_preds = []
for _ in range(20):
    state_model.decoder[2].train()
    with torch.no_grad():
        p, _, _ = state_model(test_b, test_c, test_d)
        mc_preds.append(p.numpy())
mc_preds = np.array(mc_preds) # [20, N, G]

mean_pred = np.mean(mc_preds, axis=0) # [N, G]
cond_uncertainty = np.mean(np.std(mc_preds, axis=0), axis=1) # [N]
actual_mse = np.mean((mean_pred - test_deltas)**2, axis=1) # [N]

uncert_spearman, _ = stats.spearmanr(cond_uncertainty, actual_mse)
print(f'Uncertainty vs Actual Error Spearman rho: {uncert_spearman:.4f}', flush=True)

coverages = [1.0, 0.9, 0.8, 0.7, 0.6, 0.5]
coverage_records = []
sort_idx = np.argsort(cond_uncertainty) # ascending: most confident first

for cov in coverages:
    n_keep = int(cov * len(sort_idx))
    kept_idx = sort_idx[:n_keep]
    cov_mse = float(np.mean(actual_mse[kept_idx]))
    # Random abstention baseline
    np.random.seed(42)
    rand_idx = np.random.choice(len(actual_mse), n_keep, replace=False)
    rand_mse = float(np.mean(actual_mse[rand_idx]))
    coverage_records.append({
        'coverage': cov,
        'model_risk_mse': cov_mse,
        'random_abstention_risk_mse': rand_mse
    })

df_cov = pd.DataFrame(coverage_records)
cov_out = os.path.join(base_dir, '17_uncertainty/coverage_risk_curve.tsv')
df_cov.to_csv(cov_out, sep='	', index=False)
print(df_cov.round(4))

# -------------------------------------------------------------
# 2. DUAL SUPPORT DISTANCE & RELIABILITY BOUNDARY
# -------------------------------------------------------------
print('--> Step 17: Dual Support Distance & Reliability Decay', flush=True)
from chemical_genomics import compute_cell_support_distance, compute_perturbation_support_distance

train_fps_unique = [np.array(fp) for fp in list(set([tuple(fp) for fp in train_fps]))]
train_basals_unique = [np.array(b) for b in list(set([tuple(b) for b in train_basals]))]

support_records = []
for i in range(len(test_df)):
    c_dist = compute_cell_support_distance(test_basals[i], train_basals_unique)
    p_dist = compute_perturbation_support_distance(test_fps[i], train_fps_unique)
    # Correlation for this condition
    c_r = np.corrcoef(mean_pred[i], test_deltas[i])[0, 1] if np.std(mean_pred[i]) > 1e-6 and np.std(test_deltas[i]) > 1e-6 else 0.0
    support_records.append({
        'cell_line': test_df.iloc[i]['cell_line'],
        'perturbation': test_df.iloc[i]['perturbation'],
        'cell_support_distance': c_dist,
        'pert_support_distance': p_dist,
        'joint_support_distance': np.sqrt(c_dist**2 + p_dist**2),
        'condition_pcc_delta': c_r,
        'condition_mse': actual_mse[i],
        'uncertainty': cond_uncertainty[i]
    })

df_supp = pd.DataFrame(support_records)
supp_out = os.path.join(base_dir, '15_double_OOD/support_distance_map.tsv')
df_supp.to_csv(supp_out, sep='	', index=False)
print(f'Wrote support distance map ({len(df_supp)} conditions) to {supp_out}')

# -------------------------------------------------------------
# 3. ORTHOGONAL PROTEOMICS CONCORDANCE (PROTEINTALKS)
# -------------------------------------------------------------
print('--> Step 18: Orthogonal Proteomics Validation (ProteinTalks)', flush=True)
prot_treated = prot_df[prot_df['is_control'] == False]
prot_ctrl = prot_df[prot_df['is_control'] == True]
prot_basal_dict = {cell: group[gene_cols].values[0] for cell, group in prot_ctrl.groupby('cell_line')}

ortho_records = []
for _, row in prot_treated.iterrows():
    cell = row['cell_line']
    drug = row['perturbation']
    prot_delta = row[gene_cols].values.astype(np.float32) - prot_basal_dict[cell]
    
    # Get matching RNA predicted delta from STATE model
    b_rna = basal_dict[cell]
    fp = drug_fps[drug]
    with torch.no_grad():
        p_rna_delta, _, _ = state_model(torch.tensor(b_rna, dtype=torch.float32).unsqueeze(0),
                                        torch.tensor(fp, dtype=torch.float32).unsqueeze(0),
                                        torch.tensor([1.0], dtype=torch.float32))
        p_rna_delta = p_rna_delta.squeeze(0).numpy()
        
    # Cross-modal gene-level correlation
    r_gene = np.corrcoef(p_rna_delta, prot_delta)[0, 1]
    
    # Cross-modal pathway-level concordance
    pw_rna = [np.mean(p_rna_delta[[gene2idx[g] for g in pw_genes if g in gene2idx]]) for pw_genes in hallmark_pathways.values()]
    pw_prot = [np.mean(prot_delta[[gene2idx[g] for g in pw_genes if g in gene2idx]]) for pw_genes in hallmark_pathways.values()]
    r_pathway = np.corrcoef(pw_rna, pw_prot)[0, 1]
    
    ortho_records.append({
        'cell_line': cell,
        'drug_name': drug,
        'gene_level_rna_protein_pcc': round(float(r_gene), 4),
        'pathway_level_rna_protein_pcc': round(float(r_pathway), 4)
    })

df_ortho = pd.DataFrame(ortho_records)
ortho_out = os.path.join(base_dir, '18_proteomics/proteintalks_orthogonal_validation.tsv')
df_ortho.to_csv(ortho_out, sep='	', index=False)
print(df_ortho.describe().round(4))

# -------------------------------------------------------------
# 4. DOWNSTREAM DRUG EFFICACY RANKING (PRISM) & INCREMENTAL VALUE
# -------------------------------------------------------------
print('--> Step 19: Downstream Drug Efficacy Ranking (PRISM)', flush=True)
drug_utility_records = []

# Score drug response: induction of apoptosis + repression of proliferation (E2F/MYC/G2M)
apop_idx = [gene2idx[g] for g in hallmark_pathways['HALLMARK_APOPTOSIS'] if g in gene2idx]
prolif_idx = [gene2idx[g] for g in hallmark_pathways['HALLMARK_E2F_TARGETS'] + hallmark_pathways['HALLMARK_G2M_CHECKPOINT'] if g in gene2idx]

candidate_rows = []
for _, row in prism_df.iterrows():
    cell = row['cell_line']
    drug = row['drug_name']
    true_auc = row['viability_auc']
    is_sensitive = row['is_sensitive_responder']
    
    b = basal_dict.get(cell, basal_dict['A549'])
    fp = drug_fps.get(drug, np.zeros(1024, dtype=np.float32))
    
    with torch.no_grad():
        pred_delta, _, _ = state_model(torch.tensor(b, dtype=torch.float32).unsqueeze(0),
                                       torch.tensor(fp, dtype=torch.float32).unsqueeze(0),
                                       torch.tensor([1.0], dtype=torch.float32))
        pred_delta = pred_delta.squeeze(0).numpy()
        
    # Virtual Cell In Silico Efficacy Score:
    # High score = High Apoptosis Delta - Proliferation Delta
    vc_score = float(np.mean(pred_delta[apop_idx]) - np.mean(pred_delta[prolif_idx]))
    
    # Basal expression proxy score (e.g. basal oncogene expression)
    basal_score = float(np.mean(b[prolif_idx]))
    
    # Chemical structure alone proxy score
    chem_score = float(np.sum(fp))
    
    c_dist = compute_cell_support_distance(b, train_basals_unique)
    p_dist = compute_perturbation_support_distance(fp, train_fps_unique)
    
    candidate_rows.append({
        'cell_line': cell,
        'drug_name': drug,
        'predicted_efficacy_score': vc_score,
        'basal_proxy_score': basal_score,
        'chemical_proxy_score': chem_score,
        'true_viability_auc': true_auc,
        'is_sensitive_responder': is_sensitive,
        'cell_support_dist': round(c_dist, 3),
        'pert_support_dist': round(p_dist, 3),
        'uncertainty': float(np.std(pred_delta))
    })

df_cand = pd.DataFrame(candidate_rows)

# Evaluate ranking efficacy per cell line
per_cell_utility = []
for cell, group in df_cand.groupby('cell_line'):
    r_vc, _ = stats.spearmanr(group['predicted_efficacy_score'], -group['true_viability_auc'])
    r_basal, _ = stats.spearmanr(group['basal_proxy_score'], -group['true_viability_auc'])
    r_chem, _ = stats.spearmanr(group['chemical_proxy_score'], -group['true_viability_auc'])
    
    auroc_vc = roc_auc_score(group['is_sensitive_responder'], group['predicted_efficacy_score']) if len(set(group['is_sensitive_responder'])) > 1 else 0.5
    auroc_basal = roc_auc_score(group['is_sensitive_responder'], group['basal_proxy_score']) if len(set(group['is_sensitive_responder'])) > 1 else 0.5
    auroc_chem = roc_auc_score(group['is_sensitive_responder'], group['chemical_proxy_score']) if len(set(group['is_sensitive_responder'])) > 1 else 0.5
    
    per_cell_utility.append({
        'cell_line': cell,
        'spearman_vc_score': r_vc,
        'spearman_basal_proxy': r_basal,
        'spearman_chemical_proxy': r_chem,
        'auroc_vc_retrieval': auroc_vc,
        'auroc_basal_retrieval': auroc_basal,
        'auroc_chem_retrieval': auroc_chem,
        'incremental_auroc_gain': auroc_vc - max(auroc_basal, auroc_chem)
    })

df_util = pd.DataFrame(per_cell_utility)
util_out = os.path.join(base_dir, '19_drug_utility/drug_efficacy_ranking_utility.tsv')
df_util.to_csv(util_out, sep='	', index=False)
print(df_util.mean(numeric_only=True).round(4))

# Output candidate_drugs.tsv
cand_out = os.path.join(base_dir, '19_drug_utility/candidate_drugs.tsv')
df_cand.to_csv(cand_out, sep='	', index=False)
print(f'Wrote candidate drug rankings to {cand_out}')
