import numpy as np
from scipy import stats
from sklearn.metrics import mean_squared_error, r2_score

def compute_level1_reconstruction_vectorized(y_true, y_pred):
    N, G = y_true.shape
    pccs = []
    for i in range(N):
        if np.std(y_true[i]) > 1e-6 and np.std(y_pred[i]) > 1e-6:
            pccs.append(np.corrcoef(y_true[i], y_pred[i])[0, 1])
        else:
            pccs.append(0.0)
    mse = float(np.mean((y_true - y_pred)**2))
    return {
        'pcc_expr': float(np.mean(pccs)),
        'mse_expr': mse
    }

def compute_level2_delta_metrics_vectorized(pred_deltas, true_deltas):
    N, G = true_deltas.shape
    std_true = np.std(true_deltas, axis=1)
    std_pred = np.std(pred_deltas, axis=1)
    valid_mask = (std_true > 1e-6) & (std_pred > 1e-6)
    
    mean_true = np.mean(true_deltas, axis=1, keepdims=True)
    mean_pred = np.mean(pred_deltas, axis=1, keepdims=True)
    
    pcc_deltas = []
    for i in range(N):
        if valid_mask[i]:
            pcc_deltas.append(np.corrcoef(true_deltas[i], pred_deltas[i])[0, 1])
        else:
            pcc_deltas.append(0.0)
            
    std_norm_true = (true_deltas - mean_true) / np.maximum(std_true[:, None] * np.sqrt(G), 1e-7)
    std_norm_pred = (pred_deltas - mean_pred) / np.maximum(std_pred[:, None] * np.sqrt(G), 1e-7)
    
    sim_matrix = np.dot(std_norm_pred, std_norm_true.T)
    diag_sims = np.diag(sim_matrix)
    off_diag_means = (np.sum(sim_matrix, axis=1) - diag_sims) / max(1, N - 1)
    iss = diag_sims - off_diag_means
    
    ranks = np.sum(sim_matrix > diag_sims[:, None], axis=1) + 1
    top1 = (ranks == 1).astype(float)
    top5 = (ranks <= 5).astype(float)
    
    abs_true = np.abs(true_deltas)
    abs_pred = np.abs(pred_deltas)
    jaccard_20_list = []
    sign_concordance_list = []
    
    for i in range(N):
        t20 = set(np.argsort(abs_true[i])[-20:])
        p20 = set(np.argsort(abs_pred[i])[-20:])
        jaccard_20_list.append(len(t20.intersection(p20)) / max(1, len(t20.union(p20))))
        
        t50 = np.argsort(abs_true[i])[-50:]
        sign_concordance_list.append(np.mean(np.sign(true_deltas[i, t50]) == np.sign(pred_deltas[i, t50])))
        
    return {
        'pcc_delta': float(np.mean(pcc_deltas)),
        'iss': float(np.mean(iss)),
        'top1_retrieval': float(np.mean(top1)),
        'top5_retrieval': float(np.mean(top5)),
        'deg_jaccard_20': float(np.mean(jaccard_20_list)),
        'direction_concordance': float(np.mean(sign_concordance_list))
    }

def compute_level3_mechanisms_vectorized(pred_deltas, true_deltas, gene_cols, hallmark_pathways, generic_pc1):
    N, G = true_deltas.shape
    gene2idx = {g: i for i, g in enumerate(gene_cols)}
    
    pw_true_mat = []
    pw_pred_mat = []
    for pw_name, pw_genes in hallmark_pathways.items():
        indices = [gene2idx[g] for g in pw_genes if g in gene2idx]
        if len(indices) >= 3:
            pw_true_mat.append(np.mean(true_deltas[:, indices], axis=1))
            pw_pred_mat.append(np.mean(pred_deltas[:, indices], axis=1))
            
    pw_true_mat = np.array(pw_true_mat).T
    pw_pred_mat = np.array(pw_pred_mat).T
    
    pw_corrs = []
    for i in range(N):
        if np.std(pw_true_mat[i]) > 1e-6 and np.std(pw_pred_mat[i]) > 1e-6:
            pw_corrs.append(np.corrcoef(pw_true_mat[i], pw_pred_mat[i])[0, 1])
        else:
            pw_corrs.append(0.0)
            
    pc1_norm = np.dot(generic_pc1, generic_pc1)
    if pc1_norm > 1e-8:
        betas_true = np.dot(true_deltas, generic_pc1) / pc1_norm
        betas_pred = np.dot(pred_deltas, generic_pc1) / pc1_norm
        res_true = true_deltas - betas_true[:, None] * generic_pc1
        res_pred = pred_deltas - betas_pred[:, None] * generic_pc1
        stress_frac = np.mean(betas_true**2 * pc1_norm / (np.sum(true_deltas**2, axis=1) + 1e-8))
    else:
        res_true = true_deltas
        res_pred = pred_deltas
        stress_frac = 0.0
        
    res_corrs = []
    for i in range(N):
        if np.std(res_true[i]) > 1e-6 and np.std(res_pred[i]) > 1e-6:
            res_corrs.append(np.corrcoef(res_true[i], res_pred[i])[0, 1])
        else:
            res_corrs.append(0.0)
            
    return {
        'pathway_fidelity': float(np.mean(pw_corrs)),
        'pcc_residual_delta': float(np.mean(res_corrs)),
        'generic_stress_fraction': float(stress_frac)
    }
