# 🔬 Cross-Resolution Analysis of Tumor Proliferation Signatures

This project explores tumor proliferation patterns using bulk gene expression data from NCBI GEO, with a focus on uncovering intra-tumoral heterogeneity through computational analysis.

---

## 📊 Dataset

- **Source:** NCBI GEO  
- **Accession:** GSE71646  
- **Platform:** GPL10558 (Illumina HumanHT-12 V4.0 expression beadchip)

Gene annotation was obtained from the corresponding GPL platform file.

---

## ⚙️ Workflow Overview

1. **Data preprocessing**
   - Parsed GEO series matrix file
   - Extracted numeric expression values
   - Merged with GPL annotation (probe → gene mapping)
   - Aggregated multiple probes per gene (mean expression)

2. **Exploratory analysis**
   - Examined key proliferation markers:
     - *MKI67, PCNA, CCNB1, CDK1*
   - Compared expression distributions across samples

3. **Correlation analysis**
   - Constructed Pearson correlation matrix for cell-cycle genes
   - Observed strong co-expression across proliferation-related genes

4. **Dimensionality reduction & clustering**
   - Applied PCA on cell-cycle gene expression
   - Performed KMeans clustering to identify sample subgroups

5. **Cluster-level interpretation**
   - Evaluated MKI67 expression across clusters
   - Assessed variability in proliferative states

---

## 🔥 Key Findings

- **Expected:** Tumor samples show elevated proliferation marker expression (e.g., MKI67)

- **More importantly:**
  > Unsupervised clustering reveals **distinct proliferative subgroups within samples**, indicating that proliferation is not a uniform signal.

- **Interpretation:**
  - High MKI67 cluster → highly proliferative state  
  - Lower MKI67 clusters → reduced or alternative cell-cycle activity  
  - Suggests **intra-tumoral heterogeneity in proliferation**

---

## 🧠 Biological Insight

While bulk comparisons (tumor vs normal) capture global trends, clustering based on cell-cycle genes highlights **hidden structure within the data**, suggesting that:

> Proliferation is a heterogeneous and structured feature, not a single continuous gradient.

This aligns with observations from single-cell studies where tumor cell populations exist in distinct transcriptional states.

---

## 📁 Outputs

- 📦 Boxplots and violin plots for gene expression
- 📊 Correlation heatmap of cell-cycle genes
- 📉 PCA visualization of sample structure
- 🔥 Cluster-specific MKI67 expression analysis

---

## 🛠️ Tools & Libraries

- Python  
- pandas  
- seaborn / matplotlib  
- scikit-learn  

---

## 🚀 Future Directions

- Validate findings across independent GEO datasets  
- Compare bulk-derived clusters with scRNA-seq cell states  
- Expand gene sets (Hallmark cell cycle / proliferation signatures)  
- Integrate immune markers to explore proliferation–immune interactions  

---

## 📎 Reproducibility

To run the analysis:

1. Download the dataset from GEO:
   - GSE71646 series matrix file  
   - GPL10558 annotation file  

2. Place them in the `data/` directory  

3. Run:

```bash
python Plotting.py
