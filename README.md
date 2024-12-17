## **1. Project Title**  
**Comparative Transcriptomic Profiling of COVID-19 and RSV: Insights from RNA-Seq and Weighted Gene Co-expression Network Analysis**

---

## **2. Project Description**  
This project compares the transcriptomic profiles of COVID-19 and Respiratory Syncytial Virus (RSV) using **RNA-Seq** data. Key steps include **data preprocessing**, **alignment**, **gene quantification**, and **differential expression analysis (DEA)**, followed by **Weighted Gene Co-expression Network Analysis (WGCNA)** to identify hub genes and modules associated with disease states and clinical traits.

### **Key Findings**  
- Identified **disease-specific modules** linked to COVID-19 severity and immune response.  
- Shared hub genes (**OASL**, **TXN**, **RBCK1**) between COVID-19 and RSV.  
- Best-performing clustering pipeline selected based on **silhouette scores** and **stability metrics**.  

---

## **3. Dataset Details**  
- **GEO Accession**: [GSE152418](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE152418)  
- **Samples**:  
   - **32 COVID-19**  
   - **34 Healthy**  
   - **2 Convalescent**  
- **Data Type**: RNA-Seq SRR files (68 samples).  

---

## **4. Required Files**  
1. **Input Files**:  
   - `SraRunTable.csv`: Metadata for samples.  
   - Reference Genome: `Homo_sapiens_UCSC_hg38.tar.gz`.  
   - Annotation File: `genes.gtf`.  
2. **Preprocessed Outputs**:  
   - Gene Counts: `gene_counts.txt`.  
   - Normalized Counts: `normalized_gene_counts.csv`.  
3. **Results**:  
   - `Top100_Overall_HubGenes.csv`  
   - `top_100_genes_RSV.csv`  
   - `overlap_genes.csv`  

---

## **5. Tools and Libraries**

### **Programming Languages**  
- **Bash**: Data processing, alignment, and quantification.  
- **Python**: Pipeline evaluation, clustering analysis, and visualization.  

### **Software/Modules**  
| Tool                | Purpose                                  |  
|---------------------|------------------------------------------|  
| **sra-toolkit**     | Download and convert SRR files to FASTQ. |  
| **FastQC**          | Quality control for raw/trimmed reads.   |  
| **Trim Galore**     | Trimming low-quality reads.              |  
| **HISAT2**          | Read alignment to the reference genome.  |  
| **SAMtools**        | Processing, sorting, and indexing BAM files. |  
| **featureCounts**   | Gene expression quantification.          |  

### **Python Libraries**  
- `pandas`, `numpy`, `scikit-learn`, `matplotlib`, `seaborn`, `umap-learn`  

### **R Libraries**  
- **DESeq2**: Differential expression analysis.  
- **WGCNA**: Co-expression network analysis.  
- **pheatmap**: Visualization of heatmaps.  
- **ggplot2**: Data plotting.  

---

## **6. Workflow**

### **Step 1: Data Acquisition and Preprocessing**  
1. **Download Data**:  
   - Use `SRA Toolkit` to download SRR files.  
   - Script: `download_sra_files.sh`.  
2. **Convert SRA to FASTQ**:  
   - `fasterq-dump` for conversion and compression.  
   - Script: `fasterq_conversion.sh`.  
3. **Quality Control**:  
   - `FastQC` for raw reads and trimmed reads.  
   - Script: `fastqc_raw.sh` and `fastqc_trimmed.sh`.  

### **Step 2: Read Trimming**  
- Tool: **Trim Galore**  
- Parameters:  
   - `--length 30`, `--quality 20`, `--clip_R1 10`  
- Script: `trim_galore.sh`.

### **Step 3: Reference Genome Preparation**  
- Downloaded `Homo_sapiens_UCSC_hg38.tar.gz` and `genes.gtf` from UCSC.  
- Created a HISAT2 index.  
- Script: `hisat2_indexing.sh`.  

### **Step 4: Read Alignment**  
- Tool: **HISAT2**  
- Input: Trimmed FASTQ files.  
- Output: Aligned **BAM** files.  
- Script: `hisat2_alignment.sh`.

### **Step 5: Gene Expression Quantification**  
- Tool: **featureCounts**  
- Input: BAM files, `genes.gtf`.  
- Output: `gene_counts.txt` and summary file.  
- Script: `feature_counts.sh`.

---
### **5.1 Differential Expression Analysis (DEA)**  
- Tool: **DESeq2**  
- Input: `gene_counts.txt`  
- Output: Differentially expressed genes (DEGs).  

---

### **5.2 Weighted Gene Co-expression Network Analysis (WGCNA)**  
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

## **9. Results Summary**  
- **MEyellow** and **MEgrey** modules are strongly associated with COVID-19 severity.  
- Shared hub genes (**OASL**, **TXN**, and **RBCK1**) indicate conserved antiviral and immune response pathways between COVID-19 and RSV.  
- Limited overlap of hub genes highlights unique molecular pathways for each disease.

---
