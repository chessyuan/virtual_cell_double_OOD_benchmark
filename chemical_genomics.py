import numpy as np
from rdkit import Chem
from rdkit.Chem import AllChem, DataStructs
from rdkit.Chem.Scaffolds import MurckoScaffold

def smiles_to_ecfp4(smiles, n_bits=1024):
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return np.zeros(n_bits, dtype=np.float32)
    fp = AllChem.GetMorganFingerprintAsBitVect(mol, radius=2, nBits=n_bits)
    arr = np.zeros((n_bits,), dtype=np.float32)
    DataStructs.ConvertToNumpyArray(fp, arr)
    return arr

def compute_tanimoto_similarity(fp1, fp2):
    dot = np.dot(fp1, fp2)
    denom = np.sum(fp1) + np.sum(fp2) - dot
    return dot / denom if denom > 0 else 0.0

def compute_perturbation_support_distance(target_fp, train_fps):
    if len(train_fps) == 0:
        return 1.0
    sims = [compute_tanimoto_similarity(target_fp, tfp) for tfp in train_fps]
    max_sim = max(sims)
    return float(1.0 - max_sim)

def compute_cell_support_distance(target_basal, train_basals):
    if len(train_basals) == 0:
        return 1.0
    # Normalized Euclidean distance in transcriptomic space
    dists = [np.linalg.norm(target_basal - tb) / (np.linalg.norm(target_basal) + 1e-6) for tb in train_basals]
    return float(min(dists))
