import os, json, time, sys
import numpy as np
import pandas as pd
import torch
import torch.optim as optim

from chemical_genomics import smiles_to_ecfp4, compute_perturbation_support_distance, compute_cell_support_distance
from baselines import (
    Baseline0_ControlIdentity, Baseline1_GlobalMean, Baseline2_ContextMean,
    Baseline3_DrugMean, Baseline4_LinearAdditive, Baseline5_RidgeRegression,
    Baseline6_kNN, Baseline7_MatrixFactorization
)
from models import (
    ChemCPA_Module, GEARS_Module, CellOT_Module, STATE_Module, MAP_Module, scGPT_FM_Module
)
from metrics import (
    compute_level1_reconstruction_vectorized,
    compute_level2_delta_metrics_vectorized,
    compute_level3_mechanisms_vectorized
)

base_dir = 'project_05_virtual_cell_double_OOD'
atlas_df = pd.read_parquet(os.path.join(base_dir, '05_processed_data/tahoe_fast_mode_pseudobulk.parquet'))
drug_df = pd.read_csv(os.path.join(base_dir, '06_metadata/drug_metadata.tsv'), sep='	')
cell_df = pd.read_csv(os.path.join(base_dir, '06_metadata/cell_line_metadata.tsv'), sep='	')
genes_df = pd.read_csv(os.path.join(base_dir, '06_metadata/genes_universe.tsv'), sep='	')
with open(os.path.join(base_dir, '06_metadata/hallmark_pathways.json')) as f:
    hallmark_pathways = json.load(f)

gene_cols = genes_df['gene_symbol'].tolist()
n_genes = len(gene_cols)

drug_fps = {}
for _, row in drug_df.iterrows():
    drug_fps[row['drug_name']] = smiles_to_ecfp4(row['canonical_smiles'])

ctrl_df = atlas_df[atlas_df['is_control'] == True]
basal_dict = {}
for cell, group in ctrl_df.groupby('cell_line'):
    basal_dict[cell] = group[gene_cols].values[0]

treated_df = atlas_df[atlas_df['is_control'] == False]
all_deltas = []
for idx, row in treated_df.iterrows():
    c = row['cell_line']
    expr = row[gene_cols].values.astype(np.float32)
    delta = expr - basal_dict[c]
    all_deltas.append(delta)
all_deltas = np.array(all_deltas)

u, s, vt = np.linalg.svd(all_deltas - np.mean(all_deltas, axis=0), full_matrices=False)
generic_pc1 = vt[0]
print(f'Computed Generic Stress PC1 (explains {s[0]**2 / np.sum(s**2) * 100:.1f}% of variance)', flush=True)

seeds = [42, 43, 44, 45, 46]
ood_tasks = ['ZS0', 'ZS1_P', 'ZS1_C', 'ZS2_Double_OOD', 'DRUG_OOD_Scaffold', 'DRUG_OOD_MOA']

all_benchmark_results = []
compute_benchmark_records = []

for seed in seeds:
    split_path = os.path.join(base_dir, f'07_splits/split_protocol_seed_{seed}.json')
    with open(split_path) as f:
        split_info = json.load(f)
        
    for ood in ood_tasks:
        print(f'--> [Seed {seed}] Task: {ood}', flush=True)
        
        if ood == 'ZS0':
            np.random.seed(seed)
            shuffled_idx = np.random.permutation(len(treated_df))
            split_point = int(0.8 * len(shuffled_idx))
            train_treated = treated_df.iloc[shuffled_idx[:split_point]]
            test_treated = treated_df.iloc[shuffled_idx[split_point:]]
        elif ood == 'ZS1_P':
            test_drugs = split_info['ZS1_P']['test_drugs']
            train_treated = treated_df[~treated_df['perturbation'].isin(test_drugs)]
            test_treated = treated_df[treated_df['perturbation'].isin(test_drugs)]
        elif ood == 'ZS1_C':
            test_cells = split_info['ZS1_C']['test_cells']
            train_treated = treated_df[~treated_df['cell_line'].isin(test_cells)]
            test_treated = treated_df[treated_df['cell_line'].isin(test_cells)]
        elif ood == 'ZS2_Double_OOD':
            test_drugs = split_info['ZS2_Double_OOD']['test_drugs']
            test_cells = split_info['ZS2_Double_OOD']['test_cells']
            train_treated = treated_df[(~treated_df['perturbation'].isin(test_drugs)) & (~treated_df['cell_line'].isin(test_cells))]
            test_treated = treated_df[(treated_df['perturbation'].isin(test_drugs)) & (treated_df['cell_line'].isin(test_cells))]
        elif ood == 'DRUG_OOD_Scaffold':
            test_drugs = split_info['DRUG_OOD_Scaffold']['test_drugs']
            train_treated = treated_df[~treated_df['perturbation'].isin(test_drugs)]
            test_treated = treated_df[treated_df['perturbation'].isin(test_drugs)]
        elif ood == 'DRUG_OOD_MOA':
            test_drugs = split_info['DRUG_OOD_MOA']['test_drugs']
            train_treated = treated_df[~treated_df['perturbation'].isin(test_drugs)]
            test_treated = treated_df[treated_df['perturbation'].isin(test_drugs)]

        if len(test_treated) == 0:
            continue

        train_deltas, train_X_feat, train_basals, train_fps_list, train_doses = [], [], [], [], []
        for _, r in train_treated.iterrows():
            c, d, dose = r['cell_line'], r['perturbation'], r['dose_uM']
            expr = r[gene_cols].values.astype(np.float32)
            b = basal_dict[c]
            train_deltas.append(expr - b)
            fp = drug_fps.get(d, np.zeros(1024, dtype=np.float32))
            train_fps_list.append(fp)
            train_basals.append(b)
            train_doses.append(dose)
            train_X_feat.append(np.concatenate([b, fp, [dose]]))
            
        train_deltas = np.array(train_deltas, dtype=np.float32)
        train_X_feat = np.array(train_X_feat, dtype=np.float32)
        train_basals = np.array(train_basals, dtype=np.float32)
        train_fps_list = np.array(train_fps_list, dtype=np.float32)
        train_doses = np.array(train_doses, dtype=np.float32)

        test_deltas, test_basals, test_fps_list, test_doses = [], [], [], []
        for _, r in test_treated.iterrows():
            c = r['cell_line']
            expr = r[gene_cols].values.astype(np.float32)
            b = basal_dict[c]
            test_deltas.append(expr - b)
            test_basals.append(b)
            test_fps_list.append(drug_fps.get(r['perturbation'], np.zeros(1024, dtype=np.float32)))
            test_doses.append(r['dose_uM'])
            
        test_deltas = np.array(test_deltas, dtype=np.float32)
        test_basals = np.array(test_basals, dtype=np.float32)
        test_fps_list = np.array(test_fps_list, dtype=np.float32)
        test_doses = np.array(test_doses, dtype=np.float32)
        test_expr_true = test_basals + test_deltas

        # 1. Fit baselines
        bl0 = Baseline0_ControlIdentity()
        bl1 = Baseline1_GlobalMean(); bl1.fit(train_deltas)
        bl2 = Baseline2_ContextMean(); bl2.fit(train_treated, gene_cols)
        bl3 = Baseline3_DrugMean(); bl3.fit(train_treated, gene_cols)
        bl4 = Baseline4_LinearAdditive(); bl4.fit(train_treated, gene_cols)
        bl5 = Baseline5_RidgeRegression(alpha=25.0); bl5.fit(train_X_feat, train_deltas)
        bl6 = Baseline6_kNN(k=5); bl6.fit(train_fps_list, train_deltas)

        # 2. Fit deep models (15 epochs is fast and reaches steady state on pseudobulk)
        torch.manual_seed(seed)
        t_tensors = {
            'basal': torch.tensor(train_basals, dtype=torch.float32),
            'chem': torch.tensor(train_fps_list, dtype=torch.float32),
            'dose': torch.tensor(train_doses, dtype=torch.float32),
            'delta': torch.tensor(train_deltas, dtype=torch.float32)
        }
        test_tensors = {
            'basal': torch.tensor(test_basals, dtype=torch.float32),
            'chem': torch.tensor(test_fps_list, dtype=torch.float32),
            'dose': torch.tensor(test_doses, dtype=torch.float32)
        }

        deep_models = {
            'M1_chemCPA': ChemCPA_Module(),
            'M2_GEARS': GEARS_Module(),
            'M3_CellOT': CellOT_Module(),
            'M4_STATE': STATE_Module(),
            'M5_MAP': MAP_Module(),
            'FM1_scGPT_Head': scGPT_FM_Module()
        }

        deep_train_times = {}
        for m_name, m_mod in deep_models.items():
            t_start = time.time()
            m_mod.train()
            optimizer = optim.AdamW(m_mod.parameters(), lr=1e-3, weight_decay=1e-4)
            criterion = torch.nn.MSELoss()
            
            for epoch in range(15):
                optimizer.zero_grad()
                if m_name == 'M4_STATE':
                    pred, mu, logvar = m_mod(t_tensors['basal'], t_tensors['chem'], t_tensors['dose'])
                    kl = -0.5 * torch.sum(1 + logvar - mu.pow(2) - logvar.exp()) / t_tensors['basal'].shape[0]
                    loss = criterion(pred, t_tensors['delta']) + 0.01 * kl
                else:
                    pred = m_mod(t_tensors['basal'], t_tensors['chem'], t_tensors['dose'])
                    loss = criterion(pred, t_tensors['delta'])
                loss.backward()
                optimizer.step()
                
            deep_train_times[m_name] = time.time() - t_start

        # Baseline predictions
        pred_bl0, pred_bl1, pred_bl2, pred_bl3, pred_bl4, pred_bl5, pred_bl6 = [], [], [], [], [], [], []
        for i, (_, r) in enumerate(test_treated.iterrows()):
            c, d, dose = r['cell_line'], r['perturbation'], r['dose_uM']
            b = test_basals[i]
            feat = np.concatenate([b, test_fps_list[i], [dose]])
            
            pred_bl0.append(bl0.predict_delta(c, d, dose, b))
            pred_bl1.append(bl1.predict_delta(c, d, dose, b))
            pred_bl2.append(bl2.predict_delta(c, d, dose, b))
            pred_bl3.append(bl3.predict_delta(c, d, dose, b))
            pred_bl4.append(bl4.predict_delta(c, d, dose, b))
            pred_bl5.append(bl5.predict_delta(feat))
            pred_bl6.append(bl6.predict_delta(test_fps_list[i]))

        model_predictions = {
            'BL0_ControlIdentity': np.array(pred_bl0, dtype=np.float32),
            'BL1_GlobalMean': np.array(pred_bl1, dtype=np.float32),
            'BL2_ContextMean': np.array(pred_bl2, dtype=np.float32),
            'BL3_DrugMean': np.array(pred_bl3, dtype=np.float32),
            'BL4_LinearAdditive': np.array(pred_bl4, dtype=np.float32),
            'BL5_RidgeChemGenomic': np.array(pred_bl5, dtype=np.float32),
            'BL6_kNN_Tanimoto': np.array(pred_bl6, dtype=np.float32)
        }

        # Deep model predictions
        for m_name, m_mod in deep_models.items():
            m_mod.eval()
            with torch.no_grad():
                if m_name == 'M4_STATE':
                    p, _, _ = m_mod(test_tensors['basal'], test_tensors['chem'], test_tensors['dose'])
                    model_predictions[m_name] = p.cpu().numpy()
                else:
                    p = m_mod(test_tensors['basal'], test_tensors['chem'], test_tensors['dose'])
                    model_predictions[m_name] = p.cpu().numpy()

        # Calculate mean support distances for this test split
        unique_train_fps = list(set([tuple(fp) for fp in train_fps_list]))
        unique_train_fps = [np.array(fp) for fp in unique_train_fps]
        unique_train_basals = list(set([tuple(b) for b in train_basals]))
        unique_train_basals = [np.array(b) for b in unique_train_basals]

        c_dists = [compute_cell_support_distance(test_basals[i], unique_train_basals) for i in range(len(test_treated))]
        p_dists = [compute_perturbation_support_distance(test_fps_list[i], unique_train_fps) for i in range(len(test_treated))]
        mean_c_dist = float(np.mean(c_dists))
        mean_p_dist = float(np.mean(p_dists))

        # Vectorized evaluation for all models
        for m_name, p_deltas in model_predictions.items():
            test_expr_pred = test_basals + p_deltas
            l1 = compute_level1_reconstruction_vectorized(test_expr_true, test_expr_pred)
            l2 = compute_level2_delta_metrics_vectorized(p_deltas, test_deltas)
            l3 = compute_level3_mechanisms_vectorized(p_deltas, test_deltas, gene_cols, hallmark_pathways, generic_pc1)
            
            all_benchmark_results.append({
                'seed': seed,
                'ood_task': ood,
                'model': m_name,
                'n_test_conditions': len(test_treated),
                'pcc_expr': l1['pcc_expr'],
                'mse_expr': l1['mse_expr'],
                'pcc_delta': l2['pcc_delta'],
                'iss': l2['iss'],
                'top1_retrieval': l2['top1_retrieval'],
                'top5_retrieval': l2['top5_retrieval'],
                'deg_jaccard_20': l2['deg_jaccard_20'],
                'direction_concordance': l2['direction_concordance'],
                'pathway_fidelity': l3['pathway_fidelity'],
                'pcc_residual_delta': l3['pcc_residual_delta'],
                'generic_stress_fraction': l3['generic_stress_fraction'],
                'mean_cell_support_dist': mean_c_dist,
                'mean_pert_support_dist': mean_p_dist
            })

df_bench = pd.DataFrame(all_benchmark_results)
bench_out = os.path.join(base_dir, '15_double_OOD/double_ood_benchmark_results.tsv')
df_bench.to_csv(bench_out, sep='	', index=False)
print(f'Done! Successfully saved full benchmark ({len(df_bench)} rows) to {bench_out}', flush=True)
