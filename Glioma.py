import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
from sklearn.decomposition import PCA
from sklearn.cluster import KMeans
from sklearn.preprocessing import StandardScaler

# ================================
# 1. Load GEO data
# ================================
df = pd.read_csv("GSE4290_series_matrix.txt", sep="\t", comment="!")
df.set_index("ID_REF", inplace=True)

# Select numeric columns
numeric_df = df.select_dtypes(include=["number"])





# ================================
# 2. Load annotation file
# ================================
anno_file = "GPL570-55999.txt"

header_line = None
with open(anno_file, encoding="utf-8") as f:
    for i, line in enumerate(f):
        if line.startswith("ID\t"):
            header_line = i
            break

if header_line is None:
    raise ValueError("⚠️ ID column not found")

anno = pd.read_csv(
    anno_file,
    sep="\t",
    skiprows=header_line,
    low_memory=False
)

# Detect gene column
gene_col = None
for col in anno.columns:
    if "Symbol" in col or "ILMN_Gene" in col:
        gene_col = col
        break

if gene_col is None:
    raise ValueError("⚠️ Gene column not found")

print("Using gene column:", gene_col)



# ================================
# 3. Merge expression + annotation
# ================================
df_reset = numeric_df.reset_index()

merged = df_reset.merge(
    anno[['ID', gene_col]],
    left_on='ID_REF',
    right_on='ID'
)

merged.set_index(gene_col, inplace=True)

merged_numeric = merged.select_dtypes(include=["number"])
merged_numeric.index = merged.index

merged = merged_numeric.groupby(merged_numeric.index).mean()

print("Final shape after grouping:", merged.shape)





# ================================
# 4. Extract phenotype (Tumor vs Normal)
# ================================
sample_ids = []
groups = []

with open("GSE4290_series_matrix.txt") as f:
    for line in f:
        if line.startswith("!Sample_geo_accession"):
            parts = line.strip().split("=")
            if len(parts) > 1:
                sample_ids.append(parts[1].strip())
        if line.startswith("!Sample_characteristics_ch1"):
            parts = line.strip().split("=")
            if len(parts) > 1:
                value = parts[1].strip()
                # تنظيف النصوص
                if "Tumor" in value:
                    groups.append("Tumor")
                elif "Normal" in value:
                    groups.append("Normal")
                else:
                    groups.append(value)

# بناء DataFrame متوافق
sample_info = pd.DataFrame({
    "Sample": sample_ids,
    "Group": groups
})

# فلترة بحيث تتطابق مع أعمدة merged
sample_info = sample_info[sample_info["Sample"].isin(merged.columns)]

print(sample_info["Group"].value_counts())



# ================================
# 4. MKI67 Tumor vs Normal (optional)
# ================================

mki67_values = merged.loc["MKI67"]

tumor_samples = sample_info[sample_info["Group"]=="Tumor"]["Sample"]
normal_samples = sample_info[sample_info["Group"]=="Normal"]["Sample"]

mki67_tumor = mki67_values[tumor_samples]
mki67_normal = mki67_values[normal_samples]

plt.figure(figsize=(6,5))
sns.boxplot(data=[mki67_tumor, mki67_normal])
plt.xticks([0,1], ["Tumor","Normal"])
plt.title("MKI67 Expression (Tumor vs Normal) in Glioma")
plt.savefig("Glioma images/MKI67_Tumor_vs_Normal.png")
plt.show()





# ================================
# 5. Proliferation genes analysis
# ================================
proliferation_genes = ["MKI67","PCNA","CCNB1","CDK1"]
proliferation_df = merged.loc[merged.index.intersection(proliferation_genes)]

plt.figure(figsize=(8,5))
sns.boxplot(data=proliferation_df.T)
plt.title("Proliferation Genes Expression in Glioma")
plt.savefig("Glioma images/Proliferation_Genes_Boxplot in Glioma.png")
plt.show()








# ================================
# 6. Correlation analysis
# ================================
cell_cycle_genes = [
    "MCM5", "PCNA", "TYMS", "FEN1", "MCM2", "MCM4",
    "RRM1", "MKI67", "TOP2A", "CCNB1", "CDC20",
    "CDK1", "UBE2C", "NDC80"
]

subset = merged.loc[merged.index.intersection(cell_cycle_genes)].T.dropna()

corr_matrix = subset.corr()

plt.figure(figsize=(10,8))
sns.heatmap(corr_matrix, cmap="coolwarm", center=0)
plt.title("Cell Cycle Correlation in Glioma")
plt.savefig("Glioma images/Correlation_Heatmap in Glioma.png")
plt.show()






# ================================
# 7. PCA + KMeans
# ================================
scaler = StandardScaler()
scaled = scaler.fit_transform(subset)

pca = PCA(n_components=2)
pca_result = pca.fit_transform(scaled)

print("Explained variance:", pca.explained_variance_ratio_)

kmeans = KMeans(n_clusters=3, random_state=42)
clusters = kmeans.fit_predict(pca_result)

plt.figure(figsize=(8,6))
plt.scatter(pca_result[:,0], pca_result[:,1], c=clusters, cmap="tab10")
plt.xlabel("PC1")
plt.ylabel("PC2")
plt.title("Cell-Cycle Genes PCA + KMeans in Glioma")
plt.colorbar(label="Cluster")
plt.savefig("Glioma images/PCA_KMeans in Glioma.png")
plt.show()






# ================================
# 8. 🔥 KEY STEP: MKI67 across clusters
# ================================
cluster_df = pd.DataFrame({
    "Cluster": clusters,
    "MKI67": subset["MKI67"].values
})

plt.figure(figsize=(6,5))
sns.boxplot(x="Cluster", y="MKI67", data=cluster_df)
plt.title("MKI67 Expression Across Clusters in Glioma")
plt.savefig("Glioma images/MKI67_Clusters in Glioma.png")
plt.show()