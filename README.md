# 🔬 Cross-Dataset Analysis of Tumor Proliferation Signatures

This project investigates tumor proliferation patterns using bulk RNA-seq gene expression data, with a focus on identifying intra-tumoral heterogeneity and comparing its structure across different cancer types.

---

## 📊 Datasets

- **Melanoma**
  - Source: NCBI GEO
  - Accession: GSE65904

- **Glioma**
  - Source: NCBI GEO
  - Accession: (add your dataset here)

- **Platform**
  - GPL10558 (Illumina HumanHT-12 V4.0)

Gene annotation was retrieved from the corresponding GPL platform file.

---

## ⚙️ Workflow Overview

### 1. Data Preprocessing
- Parsed GEO series matrix files
- Extracted numeric gene expression values
- Mapped probe IDs to gene symbols using GPL annotation
- Aggregated multiple probes per gene (mean expression)

---

### 2. Gene Selection
Focused on key proliferation and cell-cycle markers:

- **Core markers:**
  - MKI67, PCNA, CCNB1, CDK1

- **Extended cell-cycle set:**
  - MCM2, MCM4, MCM5, TOP2A, CDC20, UBE2C, NDC80, TYMS, FEN1, RRM1

---

### 3. Exploratory Analysis
- Boxplots and violin plots of gene expression
- Comparison of proliferation markers across samples
- Correlation heatmap of cell-cycle genes

---

### 4. Dimensionality Reduction & Clustering
- Standardization of expression values
- PCA (Principal Component Analysis)
- KMeans clustering to identify sample subgroups

---

### 5. Cluster-Level Interpretation
- Analysis of MKI67 expression across clusters
- Identification of distinct proliferative states

---

## 🔥 Key Findings

### 1. Proliferation is Not Uniform
Clustering based on cell-cycle genes reveals that:

> Tumor samples are not homogeneous — they exhibit distinct proliferative states.

---

### 2. Cross-Dataset Differences in Structure

Despite observing proliferation heterogeneity in both datasets, the **structure of this heterogeneity differs**:

- **Glioma**
  - Displays a **continuous gradient** of MKI67 expression
  - Suggests progressive variation in proliferative activity

- **Melanoma**
  - Exhibits **discrete clusters**
  - Indicates distinct proliferative states rather than a continuum

---

### 3. Key Insight

> Proliferation heterogeneity is conserved across tumor types, but its structural organization varies.

---

## 🧠 Biological Interpretation

- **Glioma**
  - May reflect gradual transitions between proliferative states
  - Suggests higher plasticity in tumor growth dynamics

- **Melanoma**
  - Likely composed of more defined subpopulations
  - Proliferation may behave as a state-switching process

---

## 📊 Outputs

- 📦 Gene expression boxplots and violin plots  
- 📊 Cell-cycle correlation heatmaps  
- 📉 PCA visualizations  
- 🔥 Cluster-specific MKI67 expression plots (per dataset)  

---

## 🛠️ Tools & Libraries

- Python  
- pandas  
- seaborn / matplotlib  
- scikit-learn  

---

## 🚀 Future Directions

- Validate findings across additional tumor datasets  
- Compare bulk-derived clusters with scRNA-seq cell states  
- Integrate immune-related genes to explore proliferation–immune interactions  
- Investigate links with clinical features (if available)  

---

## 📎 Reproducibility

### Requirements:
```bash
pip install pandas seaborn matplotlib scikit-learn
