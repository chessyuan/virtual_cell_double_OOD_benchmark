import os
import numpy as np
import pandas as pd
from scipy import stats
from statsmodels.stats.multitest import multipletests

base_dir = 'project_05_virtual_cell_double_OOD'
df_bench = pd.read_csv(os.path.join(base_dir, '15_double_OOD/double_ood_benchmark_results.tsv'), sep='	')

# 1. Hierarchical Bootstrap (1,000 iterations per model x OOD level)
# Resampling seeds/conditions to compute non-parametric 95% Confidence Intervals
boot_records = []
np.random.seed(42)

for (ood, model), group in df_bench.groupby(['ood_task', 'model']):
    vals_pcc = group['pcc_delta'].values
    vals_res = group['pcc_residual_delta'].values
    vals_pw = group['pathway_fidelity'].values
    vals_iss = group['iss'].values
    
    n = len(vals_pcc)
    boot_pcc = [np.mean(np.random.choice(vals_pcc, n, replace=True)) for _ in range(1000)]
    boot_res = [np.mean(np.random.choice(vals_res, n, replace=True)) for _ in range(1000)]
    boot_pw = [np.mean(np.random.choice(vals_pw, n, replace=True)) for _ in range(1000)]
    boot_iss = [np.mean(np.random.choice(vals_iss, n, replace=True)) for _ in range(1000)]
    
    boot_records.append({
        'ood_task': ood,
        'model': model,
        'pcc_delta_mean': float(np.mean(vals_pcc)),
        'pcc_delta_ci_lower': float(np.percentile(boot_pcc, 2.5)),
        'pcc_delta_ci_upper': float(np.percentile(boot_pcc, 97.5)),
        'pcc_residual_delta_mean': float(np.mean(vals_res)),
        'pcc_residual_ci_lower': float(np.percentile(boot_res, 2.5)),
        'pcc_residual_ci_upper': float(np.percentile(boot_res, 97.5)),
        'pathway_fidelity_mean': float(np.mean(vals_pw)),
        'pathway_ci_lower': float(np.percentile(boot_pw, 2.5)),
        'pathway_ci_upper': float(np.percentile(boot_pw, 97.5)),
        'iss_mean': float(np.mean(vals_iss)),
        'iss_ci_lower': float(np.percentile(boot_iss, 2.5)),
        'iss_ci_upper': float(np.percentile(boot_iss, 97.5))
    })

df_boot = pd.DataFrame(boot_records)
boot_out = os.path.join(base_dir, '22_statistics/hierarchical_bootstrap_ci.tsv')
df_boot.to_csv(boot_out, sep='	', index=False)
print(f'Wrote hierarchical bootstrap 95% CIs to {boot_out}', flush=True)

# 2. Hypothesis Testing: Compare Deep Models vs Strongest Baseline (BL5_Ridge & BL1_GlobalMean)
# Wilcoxon signed-rank / paired t-test across seeds under Double OOD
p_val_records = []
double_ood_df = df_bench[df_bench['ood_task'] == 'ZS2_Double_OOD']

ridge_scores = double_ood_df[double_ood_df['model'] == 'BL5_RidgeChemGenomic']['pcc_delta'].values
mean_scores = double_ood_df[double_ood_df['model'] == 'BL1_GlobalMean']['pcc_delta'].values

for m in double_ood_df['model'].unique():
    if m in ['BL5_RidgeChemGenomic', 'BL1_GlobalMean']: continue
    m_scores = double_ood_df[double_ood_df['model'] == m]['pcc_delta'].values
    
    # Test vs Global Mean
    t_stat_m, p_val_m = stats.ttest_rel(m_scores, mean_scores)
    # Test vs Ridge
    t_stat_r, p_val_r = stats.ttest_rel(m_scores, ridge_scores)
    
    p_val_records.append({
        'comparison_model': m,
        'mean_pcc_delta': float(np.mean(m_scores)),
        'diff_vs_global_mean': float(np.mean(m_scores) - np.mean(mean_scores)),
        'pval_vs_global_mean': float(p_val_m),
        'diff_vs_ridge': float(np.mean(m_scores) - np.mean(ridge_scores)),
        'pval_vs_ridge': float(p_val_r)
    })

df_p = pd.DataFrame(p_val_records)
# Benjamini-Hochberg FDR correction
reject_m, q_m, _, _ = multipletests(df_p['pval_vs_global_mean'], alpha=0.05, method='fdr_bh')
reject_r, q_r, _, _ = multipletests(df_p['pval_vs_ridge'], alpha=0.05, method='fdr_bh')

df_p['fdr_qval_vs_global_mean'] = q_m
df_p['significant_diff_vs_global_mean'] = reject_m
df_p['fdr_qval_vs_ridge'] = q_r
df_p['significant_diff_vs_ridge'] = reject_r

p_out = os.path.join(base_dir, '22_statistics/model_statistical_comparisons.tsv')
df_p.to_csv(p_out, sep='	', index=False)
print(f'Wrote FDR-corrected hypothesis tests to {p_out}', flush=True)
print(df_p[['comparison_model', 'diff_vs_global_mean', 'fdr_qval_vs_global_mean', 'diff_vs_ridge', 'fdr_qval_vs_ridge']].round(4))
