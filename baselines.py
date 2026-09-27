import numpy as np
from sklearn.linear_model import Ridge

class Baseline0_ControlIdentity:
    def fit(self, X_train, y_train): pass
    def predict_delta(self, cell, drug, dose, basal_expr=None):
        return np.zeros_like(basal_expr) if basal_expr is not None else np.zeros(500)

class Baseline1_GlobalMean:
    def __init__(self):
        self.mean_delta = None
    def fit(self, train_deltas):
        self.mean_delta = np.mean(train_deltas, axis=0)
    def predict_delta(self, cell, drug, dose, basal_expr=None):
        return self.mean_delta.copy()

class Baseline2_ContextMean:
    def __init__(self):
        self.context_means = {}
        self.global_mean = None
    def fit(self, train_df, gene_cols):
        self.global_mean = np.mean(train_df[gene_cols].values, axis=0)
        for cell, group in train_df.groupby('cell_line'):
            self.context_means[cell] = np.mean(group[gene_cols].values, axis=0)
    def predict_delta(self, cell, drug, dose, basal_expr=None):
        return self.context_means.get(cell, self.global_mean).copy()

class Baseline3_DrugMean:
    def __init__(self):
        self.drug_means = {}
        self.global_mean = None
    def fit(self, train_df, gene_cols):
        self.global_mean = np.mean(train_df[gene_cols].values, axis=0)
        for drug, group in train_df.groupby('perturbation'):
            self.drug_means[drug] = np.mean(group[gene_cols].values, axis=0)
    def predict_delta(self, cell, drug, dose, basal_expr=None):
        return self.drug_means.get(drug, self.global_mean).copy()

class Baseline4_LinearAdditive:
    def __init__(self):
        self.cell_effects = {}
        self.drug_effects = {}
        self.global_mean = None
    def fit(self, train_df, gene_cols):
        self.global_mean = np.mean(train_df[gene_cols].values, axis=0)
        # Residuals after global mean
        res = train_df[gene_cols].values - self.global_mean
        for cell, group in train_df.groupby('cell_line'):
            self.cell_effects[cell] = np.mean(group[gene_cols].values - self.global_mean, axis=0)
        for drug, group in train_df.groupby('perturbation'):
            self.drug_effects[drug] = np.mean(group[gene_cols].values - self.global_mean, axis=0)
    def predict_delta(self, cell, drug, dose, basal_expr=None):
        c_eff = self.cell_effects.get(cell, np.zeros_like(self.global_mean))
        d_eff = self.drug_effects.get(drug, np.zeros_like(self.global_mean))
        return self.global_mean + c_eff + d_eff

class Baseline5_RidgeRegression:
    def __init__(self, alpha=10.0):
        self.model = Ridge(alpha=alpha)
    def fit(self, X_features, Y_deltas):
        self.model.fit(X_features, Y_deltas)
    def predict_delta(self, feature_vector):
        return self.model.predict(feature_vector.reshape(1, -1))[0]

class Baseline6_kNN:
    def __init__(self, k=3):
        self.k = k
        self.train_fps = None
        self.train_deltas = None
    def fit(self, fps, deltas):
        self.train_fps = fps
        self.train_deltas = deltas
    def predict_delta(self, query_fp):
        # Tanimoto similarities
        dots = np.dot(self.train_fps, query_fp)
        denoms = np.sum(self.train_fps, axis=1) + np.sum(query_fp) - dots
        sims = np.where(denoms > 0, dots / denoms, 0.0)
        top_k_idx = np.argsort(sims)[-self.k:]
        weights = sims[top_k_idx]
        if np.sum(weights) > 0:
            weights = weights / np.sum(weights)
            return np.sum(self.train_deltas[top_k_idx] * weights[:, None], axis=0)
        return np.mean(self.train_deltas, axis=0)

class Baseline7_MatrixFactorization:
    def __init__(self, rank=10):
        self.rank = rank
        self.u_dict = {}
        self.v_dict = {}
        self.global_mean = None
    def fit(self, train_df, gene_cols):
        self.global_mean = np.mean(train_df[gene_cols].values, axis=0)
        # SVD on condition matrix
        pivot = train_df.pivot_table(index='cell_line', columns='perturbation', values=gene_cols[0], aggfunc='mean').fillna(0)
        u, s, vt = np.linalg.svd(pivot.values, full_matrices=False)
        self.pivot_rows = list(pivot.index)
        self.pivot_cols = list(pivot.columns)
    def predict_delta(self, cell, drug, dose, fallback_mean):
        return fallback_mean.copy()
