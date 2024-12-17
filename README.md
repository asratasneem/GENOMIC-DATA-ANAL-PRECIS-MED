## **1. Project Title**  
**Comparative Transcriptomic Profiling of COVID-19 and RSV: Insights from RNA-Seq and Weighted Gene Co-expression Network Analysis**

---

## **2. Project Description**  
This project focuses on **Precision Medicine** by analyzing transcriptomic data from patients with COVID-19 and Respiratory Syncytial Virus (RSV). Using RNA-Seq data and **Weighted Gene Co-expression Network Analysis (WGCNA)**, the project identifies:  
- **Key gene modules** associated with COVID-19 severity.  
- **Shared and distinct molecular pathways** between COVID-19 and RSV.  
- **Potential biomarkers and therapeutic targets** for severe respiratory viral infections.  

The findings provide insights into immune responses, platelet dysfunction, and antiviral pathways.

---

## **3. Dataset Details**  
- **GEO Accession**: [GSE152418](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE152418)  
- **Samples**:  
   - **32 COVID-19 patients**  
   - **34 Healthy controls**  
   - **2 Convalescent individuals**  
- **Data Type**: RNA-Seq gene counts and metadata.  

---

## **4. Objectives**  
1. Analyze transcriptomic profiles of COVID-19 and RSV using RNA-Seq.  
2. Identify gene co-expression modules correlated with COVID-19 severity.  
3. Determine shared and unique hub genes/pathways between COVID-19 and RSV.  
4. Highlight potential biomarkers and therapeutic targets.  

---

## **5. Workflow**

### **5.1 Data Acquisition and Preprocessing**  
1. **Data Download**:  
   - RNA-Seq samples downloaded using **SRA Toolkit** (fasterq-dump).  
   - Metadata file: `SraRunTable.csv`.  
2. **Quality Control**:  
   - **FastQC**: Quality assessment of raw reads.  
   - **Trim Galore**: Adapter trimming and low-quality base removal.  
   - Post-trimming QC confirmed with **FastQC**.  
3. **Alignment**:  
   - Tool: **HISAT2**  
   - Reads aligned to the **UCSC hg38 reference genome**.  
   - BAM files sorted and indexed using **SAMtools**.  
4. **Quantification**:  
   - **featureCounts**: Quantified gene expression levels to generate `gene_counts.txt`.  

---

### **5.2 Differential Expression Analysis (DEA)**  
- Tool: **DESeq2**  
- Input: `gene_counts.txt`  
- Output: Differentially expressed genes (DEGs).  

---

### **5.3 Weighted Gene Co-expression Network Analysis (WGCNA)**  
1. **Input**: Normalized gene counts (from DESeq2).  
2. **Steps**:  
   - Soft-threshold power selection for network construction.  
   - Detection of co-expression modules using hierarchical clustering.  
   - Correlation of modules with clinical traits (e.g., Severity, Gender).  
   - Identification of **hub genes** (Module Membership > 0.9).  
3. **Output**:  
   - Significant modules correlated with COVID-19 severity and healthy controls.  
   - List of hub genes for priority modules.  

---

## **6. Outputs**

### **Key Results**  
1. **Top Hub Genes**:  
   - Overall: `Top100_Overall_HubGenes.csv`  
   - RSV-specific: `top_100_genes_RSV.csv`  
2. **Shared Hub Genes**:  
   - Overlap: `overlap_genes.csv`.  

### **Visual Outputs**  
- **Heatmaps**:  
   - `heatmap_new_severity.png` – Module-trait correlations (Severity).  
   - `heatmap_new_gender.png` – Module-trait correlations (Gender).  
   - `heatmap_condition_new.png` – Module-trait correlations (Condition).  
- **Cluster Dendrogram**:  
   - `dendrogram.png` – Gene clustering with module colors.  
- **Soft-threshold Plots**:  
   - `Rplot02.pdf` – Soft-threshold model fit and connectivity.  

---

## **7. Code**  

### **R Script**  
- **File**: `WCGA.R`  
- **Purpose**:  
   - Differential expression analysis (DEA) with **DESeq2**.  
   - WGCNA: Co-expression module construction and hub gene identification.  
   - Module-trait correlation analysis and visualization.  

---

## **8. Dependencies**

### **Software**  
- **SRA Toolkit** (fasterq-dump)  
- **FastQC**: Quality control.  
- **Trim Galore**: Trimming low-quality reads.  
- **HISAT2**: Alignment to hg38.  
- **SAMtools**: BAM file processing.  
- **featureCounts**: Gene quantification.  

### **R Libraries**  
- `DESeq2`  
- `WGCNA`  
- `pheatmap`  
- `ggplot2`  

---

## **10. Execution Steps**

1. **Preprocessing and Alignment**:  
   Run Bash scripts sequentially:  
   ```bash
   bash download_sra_files.sh
   bash trim_galore.sh
   bash hisat2_alignment.sh
   bash feature_counts.sh
   ```

2. **Analysis**:  
   Execute the R script for DEA and WGCNA:  
   ```bash
   Rscript WCGA.R
   ```

---

