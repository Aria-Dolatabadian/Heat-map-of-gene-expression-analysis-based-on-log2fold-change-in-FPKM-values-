#FPKM: fragments per kilobase of transcript per million fragments mapped reads 
#Paper link: https://www.sciencedirect.com/science/article/pii/S2667064X26000709?via%3Dihub


import pandas as pd
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
from matplotlib.colors import LinearSegmentedColormap
import seaborn as sns

# ── 1. Load data ──────────────────────────────────────────────────────────────
INPUT_FILE = "gene_expression_input_data.xlsx"   # <- change path if needed

df_raw = pd.read_excel(INPUT_FILE, sheet_name="Expression Data")

# Metadata columns
gene_col      = "Gene_ID"
tf_col        = "TF_Family"
cluster_col   = "Cluster"
meta_cols     = [gene_col, tf_col, cluster_col]

# Expression matrix (remaining columns are samples)
sample_cols   = [c for c in df_raw.columns if c not in meta_cols]
df_expr       = df_raw.set_index(gene_col)[sample_cols].astype(float)

# ── 2. TF_Family colour mapping ───────────────────────────────────────────────
tf_family_colors = {
    "-":             "#A9A9A9",
    "Aa_trans":      "#FF69B4",
    "AAA":           "#90EE90",
    "NB-ARC":        "#FFD700",
    "Pkinase_Tyr":   "#FF8C00",
    "Transposase_21":"#9370DB",
    "ubiquitin":     "#4169E1",
    "UDPGT":         "#00CED1",
    "UvrD-helicase": "#2E8B57",
}

row_colors = (
    df_raw.set_index(gene_col)[tf_col]
          .map(tf_family_colors)
          .fillna("#A9A9A9")
          .rename("TF_Family")
)

# ── 3. Custom red-white-blue colormap ─────────────────────────────────────────
rwb = LinearSegmentedColormap.from_list(
    "rwb", ["#2166AC", "#92C5DE", "#FFFFFF", "#F4A582", "#D6604D", "#B2182B"]
)

# ── 4. Draw clustermap ────────────────────────────────────────────────────────
g = sns.clustermap(
    df_expr,
    cmap=rwb,
    vmin=-5, vmax=5,
    row_colors=row_colors,
    col_cluster=True,
    row_cluster=True,
    linewidths=0,
    figsize=(12, 18),
    dendrogram_ratio=(0.12, 0.08),
    colors_ratio=0.03,
    cbar_pos=(1.02, 0.35, 0.025, 0.25),
    yticklabels=True,
    xticklabels=True,
    method="average",
    metric="euclidean",
)

# ── 5. Style axes ─────────────────────────────────────────────────────────────
g.ax_heatmap.set_ylabel("Genes", fontsize=13, fontweight="bold")
g.ax_heatmap.set_xlabel("")
g.ax_heatmap.tick_params(axis="y", labelsize=6.5, right=True, left=False)
g.ax_heatmap.tick_params(axis="x", labelsize=10, rotation=45)
g.ax_heatmap.yaxis.set_label_position("left")

g.ax_col_dendrogram.set_title("Genotypes", fontsize=13, fontweight="bold", pad=6)

g.cax.set_ylabel("log$_2$ FPKM", fontsize=9)
g.cax.yaxis.set_label_position("right")
g.cax.tick_params(labelsize=8)

# ── 6. TF_Family legend ───────────────────────────────────────────────────────
# Only show families that actually appear in the data
present_tfs = df_raw[tf_col].unique()
legend_patches = [
    mpatches.Patch(color=col, label=name)
    for name, col in tf_family_colors.items()
    if name in present_tfs
]
g.ax_heatmap.legend(
    handles=legend_patches,
    title="TF_Family",
    title_fontsize=9,
    fontsize=8,
    loc="upper left",
    bbox_to_anchor=(1.18, 1.02),
    frameon=True,
    borderpad=0.8,
)

# ── 7. Save ───────────────────────────────────────────────────────────────────
out_path = "gene_expression_heatmap.png"
plt.savefig(out_path, dpi=180, bbox_inches="tight")
plt.close()
print(f"Heatmap saved → {out_path}")
print(f"  Genes: {df_expr.shape[0]}  |  Samples: {df_expr.shape[1]}")
