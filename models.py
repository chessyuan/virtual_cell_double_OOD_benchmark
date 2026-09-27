import torch
import torch.nn as nn
import torch.optim as optim
import numpy as np

# 1. chemCPA / CPA: Compositional Perturbation Autoencoder
class ChemCPA_Module(nn.Module):
    def __init__(self, gene_dim=500, chem_dim=1024, latent_dim=64):
        super().__init__()
        self.cell_encoder = nn.Sequential(
            nn.Linear(gene_dim, 128),
            nn.LayerNorm(128),
            nn.ReLU(),
            nn.Linear(128, latent_dim)
        )
        self.drug_encoder = nn.Sequential(
            nn.Linear(chem_dim + 1, 128),
            nn.LayerNorm(128),
            nn.ReLU(),
            nn.Linear(128, latent_dim)
        )
        self.decoder = nn.Sequential(
            nn.Linear(latent_dim * 2, 256),
            nn.LayerNorm(256),
            nn.ReLU(),
            nn.Dropout(0.1),
            nn.Linear(256, gene_dim)
        )
        
    def forward(self, basal_expr, chem_fp, dose):
        z_cell = self.cell_encoder(basal_expr)
        drug_in = torch.cat([chem_fp, dose.unsqueeze(-1)], dim=-1)
        z_drug = self.drug_encoder(drug_in)
        z_comb = torch.cat([z_cell, z_drug], dim=-1)
        delta_pred = self.decoder(z_comb)
        return delta_pred

# 2. GEARS: Gene Regulatory Graph-Informed Neural Network
class GEARS_Module(nn.Module):
    def __init__(self, gene_dim=500, chem_dim=1024, hidden_dim=128):
        super().__init__()
        self.target_projector = nn.Sequential(
            nn.Linear(chem_dim, hidden_dim),
            nn.ReLU(),
            nn.Linear(hidden_dim, gene_dim)
        )
        self.regulatory_filter = nn.Sequential(
            nn.Linear(gene_dim, gene_dim),
            nn.Tanh()
        )
        self.dose_scale = nn.Sequential(
            nn.Linear(1, 32),
            nn.ReLU(),
            nn.Linear(32, 1),
            nn.Sigmoid()
        )
        
    def forward(self, basal_expr, chem_fp, dose):
        init_signal = self.target_projector(chem_fp)
        propagated = self.regulatory_filter(init_signal)
        scale = self.dose_scale(dose.unsqueeze(-1))
        delta_pred = propagated * scale
        return delta_pred

# 3. CellOT / CMonge: Neural Optimal Transport
class CellOT_Module(nn.Module):
    def __init__(self, gene_dim=500, chem_dim=1024, hidden_dim=128):
        super().__init__()
        self.transport_net = nn.Sequential(
            nn.Linear(gene_dim + chem_dim + 1, hidden_dim),
            nn.ReLU(),
            nn.Linear(hidden_dim, hidden_dim),
            nn.ReLU(),
            nn.Linear(hidden_dim, gene_dim)
        )
        
    def forward(self, basal_expr, chem_fp, dose):
        x_in = torch.cat([basal_expr, chem_fp, dose.unsqueeze(-1)], dim=-1)
        transport_map = self.transport_net(x_in)
        return transport_map

# 4. STATE: Stochastic Perturbation Autoencoder with Monte Carlo Uncertainty
class STATE_Module(nn.Module):
    def __init__(self, gene_dim=500, chem_dim=1024, latent_dim=64):
        super().__init__()
        self.encoder = nn.Sequential(
            nn.Linear(gene_dim + chem_dim + 1, 128),
            nn.ReLU(),
            nn.Linear(128, latent_dim * 2)
        )
        self.decoder = nn.Sequential(
            nn.Linear(latent_dim, 128),
            nn.ReLU(),
            nn.Dropout(0.2), # Dropout enabled during inference for MC uncertainty
            nn.Linear(128, gene_dim)
        )
        
    def forward(self, basal_expr, chem_fp, dose):
        x_in = torch.cat([basal_expr, chem_fp, dose.unsqueeze(-1)], dim=-1)
        stats = self.encoder(x_in)
        mu, logvar = torch.chunk(stats, 2, dim=-1)
        std = torch.exp(0.5 * logvar)
        eps = torch.randn_like(std)
        z = mu + eps * std
        delta = self.decoder(z)
        return delta, mu, logvar

# 5. MAP: Knowledge-Informed Virtual Cell
class MAP_Module(nn.Module):
    def __init__(self, gene_dim=500, chem_dim=1024, hidden_dim=128):
        super().__init__()
        self.cell_enc = nn.Linear(gene_dim, hidden_dim)
        self.drug_enc = nn.Linear(chem_dim + 1, hidden_dim)
        # Cross attention between drug embedding and cell context
        self.attn = nn.MultiheadAttention(embed_dim=hidden_dim, num_heads=4, batch_first=True)
        self.head = nn.Sequential(
            nn.Linear(hidden_dim, 256),
            nn.ReLU(),
            nn.Linear(256, gene_dim)
        )
        
    def forward(self, basal_expr, chem_fp, dose):
        c = self.cell_enc(basal_expr).unsqueeze(1)
        d_in = torch.cat([chem_fp, dose.unsqueeze(-1)], dim=-1)
        d = self.drug_enc(d_in).unsqueeze(1)
        attn_out, _ = self.attn(c, d, d)
        delta_pred = self.head(attn_out.squeeze(1))
        return delta_pred

# 6. scGPT / Foundation Model Perturbation Head
class scGPT_FM_Module(nn.Module):
    def __init__(self, gene_dim=500, chem_dim=1024, hidden_dim=256):
        super().__init__()
        # Simulated transformer-derived gene embeddings & projection head
        self.proj = nn.Sequential(
            nn.Linear(gene_dim + chem_dim + 1, hidden_dim),
            nn.LayerNorm(hidden_dim),
            nn.GELU(),
            nn.Linear(hidden_dim, gene_dim)
        )
        
    def forward(self, basal_expr, chem_fp, dose):
        x = torch.cat([basal_expr, chem_fp, dose.unsqueeze(-1)], dim=-1)
        return self.proj(x)
