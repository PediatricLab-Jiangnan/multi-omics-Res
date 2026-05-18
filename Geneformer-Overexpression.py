import os
import pickle
import numpy as np
import pandas as pd
import scanpy as sc
import torch
from transformers import AutoModel
from torch.nn.functional import cosine_similarity
from tqdm import tqdm

# ==========================================
# 1. Configuration and Parameters
# ==========================================
model_dir = r"C:\Geneformer\Geneformer-V2-316M"
name_dict_path = r"C:\Geneformer\gene_name_id_dict_gc30M.pkl"
token_dict_path = r"C:\Geneformer\token_dictionary_gc104M.pkl"
h5ad_path = "your_object.h5ad" 

OE_GENE = "SQSTM1"
SAMPLE_SIZE = 5000
MAX_LEN = 2048 

# Setup GPU device
device = torch.device("cuda" if torch.cuda.is_available() else "cpu")
print(f"Using device: {device}")

# ==========================================
# 2. Load Model and Resources (Move to GPU)
# ==========================================
print("Loading Geneformer and dictionaries...")
with open(name_dict_path, "rb") as f:
    name_to_id = pickle.load(f)
with open(token_dict_path, "rb") as f:
    token_dict = pickle.load(f)

# Load model and transfer to GPU
model = AutoModel.from_pretrained(model_dir, output_hidden_states=True)
model.to(device)
model.eval()

# Get Target Gene Token
target_ens = name_to_id.get(OE_GENE)
oe_token_id = token_dict.get(target_ens)
if not oe_token_id:
    raise ValueError(f"Gene {OE_GENE} not found in Geneformer vocabulary.")

# ==========================================
# 3. Read Data and Filter Cells
# ==========================================
if not os.path.exists(h5ad_path):
    raise FileNotFoundError(f"File not found: {h5ad_path}")

adata = sc.read_h5ad(h5ad_path)

# Sampling function
def get_sub_adata(full_adata, condition_dict=None, size=SAMPLE_SIZE):
    if condition_dict:
        # Filter by metadata (e.g., cell type)
        mask = np.ones(full_adata.n_obs, dtype=bool)
        for k, v in condition_dict.items():
            mask = mask & (full_adata.obs[k] == v)
        sub = full_adata[mask].copy()
    else:
        sub = full_adata.copy()
    
    if sub.n_obs > size:
        idx = np.random.choice(sub.obs_names, size, replace=False)
        return sub[idx].copy()
    return sub

print("Sampling cells...")
cells_random = get_sub_adata(adata)
# Extract Neurons based on previous annotation results
cells_neurons = get_sub_adata(adata, {'cell_type_auto': 'Neurons'}) if 'cell_type_auto' in adata.obs else None

# ==========================================
# 4. GPU-Optimized Perturbation Function
# ==========================================
def run_perturbation_gpu(adata_sub, label):
    if adata_sub is None or adata_sub.n_obs == 0:
        print(f"Skip {label}: No cells found.")
        return []

    print(f"\nProcessing {label} group ({adata_sub.n_obs} cells) on GPU...")
    
    # Pre-filter genes: Keep only those present in the model dictionary
    genes_in_model = [g for g in adata_sub.var_names if name_to_id.get(g) in token_dict]
    adata_filtered = adata_sub[:, genes_in_model].copy()
    
    all_global_shifts = []
    
    # Disable gradient calculation to save VRAM
    with torch.no_grad():
        for i in tqdm(range(adata_filtered.n_obs)):
            # Extract expression data
            cell_data = adata_filtered.X[i]
            if hasattr(cell_data, "toarray"): cell_data = cell_data.toarray()[0]
            
            # --- Rank-value Encoding ---
            nonzero_idx = np.where(cell_data > 0)[0]
            if len(nonzero_idx) < 10: continue # Filter cells with too few genes
            
            nonzero_genes = adata_filtered.var_names[nonzero_idx]
            nonzero_vals = cell_data[nonzero_idx]
            
            # Sort by expression descending
            sorted_indices = np.argsort(-nonzero_vals)
            sorted_genes = nonzero_genes[sorted_indices]
            
            # Convert to Token IDs
            input_ids_ctrl = [token_dict[name_to_id[g]] for g in sorted_genes][:MAX_LEN]
            
            # --- Create Tensors and move to GPU ---
            tensor_ctrl = torch.tensor([input_ids_ctrl]).to(device)
            
            # --- Control State Inference ---
            out_ctrl = model(tensor_ctrl).last_hidden_state.squeeze(0)
            cell_emb_ctrl = torch.mean(out_ctrl, dim=0)
            
            # --- Construct Perturbed (Overexpression) State ---
            # Elevate OE_GENE to Rank 1 (first position)
            if oe_token_id in input_ids_ctrl:
                input_ids_oe = [oe_token_id] + [t for t in input_ids_ctrl if t != oe_token_id]
            else:
                input_ids_oe = [oe_token_id] + input_ids_ctrl
            
            input_ids_oe = input_ids_oe[:MAX_LEN]
            tensor_oe = torch.tensor([input_ids_oe]).to(device)
            
            # --- OE State Inference ---
            out_oe = model(tensor_oe).last_hidden_state.squeeze(0)
            cell_emb_oe = torch.mean(out_oe, dim=0)
            
            # --- Calculate Shift ---
            # Vector shift = 1 - Cosine Similarity
            shift = 1 - cosine_similarity(cell_emb_ctrl.unsqueeze(0), cell_emb_oe.unsqueeze(0)).item()
            all_global_shifts.append(shift)
            
            # Optional: Clear cache if VRAM is extremely limited
            # if i % 100 == 0: torch.cuda.empty_cache()

    return all_global_shifts

# ==========================================
# 5. Execution and Output
# ==========================================
results = {}
results['Random_5000'] = run_perturbation_gpu(cells_random, "Random Sample")

if cells_neurons is not None:
    results['Neurons'] = run_perturbation_gpu(cells_neurons, "Neurons")

# Consolidate results
df_res = pd.DataFrame(dict([(k, pd.Series(v)) for k, v in results.items()]))
df_res.to_csv(f"perturbation_shifts_{OE_GENE}_GPU.csv", index=False)

# Visualization
import matplotlib.pyplot as plt
import seaborn as sns

plt.figure(figsize=(10, 6))
sns.boxplot(data=df_res, palette="Set2")
sns.stripplot(data=df_res, color="black", size=2, alpha=0.3)
plt.title(f"Global Cell State Shift Analysis ({OE_GENE} Overexpression)\nPowered by GPU Acceleration", fontsize=14)
plt.ylabel("Impact Score (1 - Cosine Similarity)")
plt.xlabel("Cell Groups")
plt.grid(axis='y', linestyle='--', alpha=0.7)
plt.tight_layout()
plt.savefig(f"shift_comparison_{OE_GENE}_GPU.pdf")

print(f"\n[Success] Data saved to perturbation_shifts_{OE_GENE}_GPU.csv")
print(f"[Success] Visualization saved to shift_comparison_{OE_GENE}_GPU.pdf")
