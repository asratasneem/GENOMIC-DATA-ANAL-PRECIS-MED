# Install WGCNA if not installed
if (!requireNamespace("WGCNA", quietly = TRUE)) {
  install.packages("WGCNA")
}

# Load required packages
library(DESeq2)
library(WGCNA)
library(pheatmap)

# Enable multi-threading in WGCNA
options(WGCNA.allowWGCNAThreads = TRUE)


# Load gene count data
gene_counts <- read.table("gene_counts.txt", header = TRUE, row.names = 1, sep = "\t")

# Load sample traits metadata
sample_traits <- read.csv("SraRunTable.csv", row.names = 1)

# Ensure sample names match between counts and metadata
all(colnames(gene_counts) == rownames(sample_traits))

# Original column names
colnames(gene_counts)

# Clean column names to extract only sample identifiers (e.g., SRR12007878)
colnames(gene_counts) <- gsub(".*SRR([0-9]+).*", "SRR\\1", colnames(gene_counts))
# Remove the first 5 columns
gene_counts <- gene_counts[, -c(2:5)]
gene_counts <- gene_counts[, -c(1)]
# Ensure that gene_counts is numeric
gene_counts <- as.data.frame(lapply(gene_counts, as.numeric))
# Step 1: Remove the ".1" and other suffixes from column names
colnames(gene_counts) <- gsub("\\.\\d+", "", colnames(gene_counts))

# Check the updated column names
print(colnames(gene_counts))

# Step 2: Combine the duplicated columns by summing their values
merged_gene_counts <- sapply(unique(colnames(gene_counts)), function(col) {
  # Sum values of columns with the same name
  rowSums(gene_counts[, colnames(gene_counts) == col, drop = FALSE])
})
write.csv(merged_gene_counts, "merged_gene_counts.csv")
merged_gene_counts <- read.csv("merged_gene_counts.csv", row.names = 1)

# Convert the result back to a data frame
merged_gene_counts <- as.data.frame(merged_gene_counts)


# Step 4: Print the merged gene counts data
print(merged_gene_counts)

# Ensure sample names match between counts and metadata
all(colnames(merged_gene_counts) == rownames(sample_traits))

# Create DESeq2 dataset
dds <- DESeqDataSetFromMatrix(countData = merged_gene_counts,
                              colData = sample_traits,
                              design = ~ Condition)
# Keep genes with at least 10 counts across all samples
dds <- dds[rowSums(counts(dds)) > 10, ]

# Perform normalization
dds <- estimateSizeFactors(dds)

# Extract normalized counts
normalized_counts <- counts(dds, normalized = TRUE)

# Save normalized counts to a file
write.csv(normalized_counts, "normalized_counts.csv")

# Load normalized counts for WGCNA
exprData <- read.csv("normalized_counts.csv", row.names = 1)

# Transpose the data to ensure genes are rows and samples are columns
exprData <- t(exprData)

# Ensure metadata row names match expression data column names
sample_traits <- sample_traits[match(rownames(exprData), rownames(sample_traits)), ]

# Calculate variance for each gene
gene_variances <- apply(exprData, 1, var)

# Select the top variable genes (e.g., top 500)
top_genes <- order(gene_variances, decreasing = TRUE)[1:10000]

# Subset the expression data
exprData_variable_1 <- exprData[top_genes, ]
# Remove rows with NA values
exprData_variable_1 <- exprData_variable_1[complete.cases(exprData_variable_1), ]


# Candidate powers
powers <- c(c(1:10, seq(12, 50, by = 2)))


# Calculate soft-threshold
sft <- pickSoftThreshold(exprData_variable_1, powerVector = powers, networkType = "signed", verbose = 5)
library(ggplot2)

# Extract relevant data from the sft object
fitIndices <- sft$fitIndices
powers <- fitIndices[, 1]
scaleFreeTopology <- -sign(fitIndices[, 3]) * fitIndices[, 2]
meanConnectivity <- fitIndices[, 5]

# Create a data frame for ggplot
sft_df <- data.frame(
  Power = powers,
  ScaleFreeTopology = scaleFreeTopology,
  MeanConnectivity = meanConnectivity
)

# ggplot for Scale Free Topology Model Fit
p1 <- ggplot(sft_df, aes(x = Power, y = ScaleFreeTopology)) +
  geom_line(color = "blue") +
  geom_point(color = "red") +
  geom_text(aes(label = round(ScaleFreeTopology, 2)), hjust = -0.5, vjust = -0.5) +
  labs(
    title = "Scale Free Topology Model Fit",
    x = "Soft Threshold (Power)",
    y = "Model Fit (R^2)"
  ) +
  theme_minimal()

# ggplot for Mean Connectivity
p2 <- ggplot(sft_df, aes(x = Power, y = MeanConnectivity)) +
  geom_line(color = "darkgreen") +
  geom_point(color = "red") +
  geom_text(aes(label = round(MeanConnectivity, 2)), hjust = -0.5, vjust = -0.5) +
  labs(
    title = "Mean Connectivity",
    x = "Soft Threshold (Power)",
    y = "Mean Connectivity"
  ) +
  theme_minimal()

# Display both plots
library(gridExtra)
grid.arrange(p1, p2, nrow = 1)
while (dev.cur() > 1) dev.off()

softPower <- 14
enableWGCNAThreads()
exprData_variable_1 <- t(exprData_variable_1)

# Blockwise network construction and module detection
net <- blockwiseModules(
  datExpr = exprData_variable_1,        # Expression data (genes as rows, samples as columns)
  power = softPower,                    # Soft-threshold power
  TOMType = "signed",                   # Use a signed topological overlap matrix
  minModuleSize = 30,                   # Minimum module size
  reassignThreshold = 0,                # Threshold for reassigning genes to other modules
  mergeCutHeight = 0.25,                # Threshold for merging similar modules
  numericLabels = TRUE,                 # Use numeric labels instead of colors
  pamRespectsDendro = FALSE,            # Partitioning around medoids
  saveTOMs = TRUE,                      # Save topological overlap matrix (optional)
  saveTOMFileBase = "TOM",              # Base name for saved TOM files
  verbose = 3                           # Verbosity level
)
length(moduleColors) == nrow(exprData_variable_1)

#Save the module colors and labels
moduleColors <- labels2colors(net$colors)  # Convert numeric labels to color labels
table(moduleColors) 

# Save module assignments to a file
write.csv(data.frame(Gene = rownames(exprData_variable_1), Module = moduleColors), "Module_Assignments.csv")
# Transpose the exprData_variable_1 data
exprData_variable_1_transposed <- t(exprData_variable_1)

# Check dimensions
dim(exprData_variable_1)       # Should match the length of moduleColors
length(moduleColors)
nrow(exprData_variable_1)
df <- data.frame(
  Gene = rownames(exprData_variable_1),
  Module = moduleColors
)
write.csv(df, "Module_Assignments.csv")
# Make sure moduleColors is assigned to genes (columns)
if (length(moduleColors) == ncol(exprData_variable_1)) {
  # Create a data frame associating genes with their corresponding module colors
  module_associations <- data.frame(Gene = colnames(exprData_variable_1), Module = moduleColors)
  head(module_associations)
} else {
  stop("The length of moduleColors does not match the number of genes (columns) in exprData_variable_1")
}

# Set the output file as a PDF
pdf("dendrogram_with_module_colors_2.pdf", width = 8, height = 6)  # Adjust width and height as needed

# Plot dendrogram and module colors
plotDendroAndColors(
  dendro = net$dendrograms[[1]],
  colors = moduleColors[net$blockGenes[[1]]],
  groupLabels = "Module Colors",
  dendroLabels = FALSE,
  hang = 0.03,
  addGuide = TRUE,
  guideHang = 0.05
)

# Close the PDF device
dev.off()
# Get the count of genes per module
module_counts <- table(net$colors)

# Print the counts for each module
print(module_counts)
# Summary of network properties
print(net)




MEs <- moduleEigengenes(exprData_variable_1, moduleColors)$eigengenes
MEs <- orderMEs(MEs) 

rownames(MEs)

rownames(sample_traits)
print(colnames(sample_traits))
sample_traits <- sample_traits[rownames(MEs), ]

dim(MEs)
dim(sample_traits)
all(rownames(MEs) == rownames(sample_traits))
sample_traits <- sample_traits[rownames(MEs), ]
sample_traits <- sample_traits[rownames(MEs), ]
# Ensure that both sample_traits and MEs have the same row names
MEs <- MEs[match(rownames(sample_traits), rownames(MEs)), ]

# Convert categorical columns into factors and then into numeric
sample_traits$Condition <- as.factor(sample_traits$Condition)
sample_traits$Gender <- as.factor(sample_traits$Gender)
sample_traits$Severity <- as.factor(sample_traits$Severity)

# Convert factor levels into numeric values
sample_traits$Condition <- as.numeric(sample_traits$Condition)
sample_traits$Gender <- as.numeric(sample_traits$Gender)
sample_traits$Severity <- as.numeric(sample_traits$Severity)

sample_traits <- sample_traits[rownames(MEs), ]
moduleTraitCor <- cor(MEs, sample_traits, use = "pairwise.complete.obs")
moduleTraitPvalue <- corPvalueStudent(moduleTraitCor, nSamples = nrow(exprData_variable_1))



# Open a new plotting device
dev.new(width = 10, height = 8)

# Save the plot to a file (e.g., PNG)
png("heatmap_plot.png", width = 1300, height = 1000)



# Open a PDF file with smaller dimensions
pdf("heatmap_plot.pdf", width = 13, height = 12)  # Reduced dimensions

# Re-run the heatmap 
textMatrix <- paste(signif(moduleTraitCor, 2), "\n(", signif(moduleTraitPvalue, 1), ")", sep = "")
labeledHeatmap(
  Matrix = moduleTraitCor,
  xLabels = colnames(sample_traits),
  yLabels = names(MEs),
  ySymbols = names(MEs),
  colorLabels = FALSE,
  colors = blueWhiteRed(50),
  textMatrix = textMatrix,
  setStdMargins = TRUE,
  yLabelsPosition = "left",
  cex.text = 1,
    # Column label size (inside labeledHeatmap)
  cex.lab = 1.1,  
  cex.lab.y = 1.2,
  cex.lab.x = 1.5,# Axis labels size (inside labeledHeatmap)
  cex.legend = 1.2,  # Legend size (inside labeledHeatmap)
  zlim = c(-1, 1)
)

# Close the PDF device
dev.off()

# Close the file device
dev.off()







# Ensure that both sample_traits and MEs have the same row names
MEs <- MEs[match(rownames(sample_traits), rownames(MEs)), ]

# Convert the Condition column into a factor with the appropriate labels
sample_traits$Condition <- factor(sample_traits$Condition, 
                                  levels = c(1, 2, 3), 
                                  labels = c("Convalescent", "COVID-19", "Healthy"))

# Subset the data to focus on the Condition column
condition_data <- sample_traits[, "Condition", drop = FALSE]

# Ensure that the row names of condition_data match MEs
condition_data <- condition_data[rownames(MEs), ]

# Create separate columns for each condition
condition_levels <- levels(sample_traits$Condition)
condition_matrix <- as.data.frame(matrix(NA, nrow = nrow(MEs), ncol = length(condition_levels)))

# Assign column names as the condition levels (i.e., Healthy, COVID-19, Convalescent)
colnames(condition_matrix) <- condition_levels

# Fill the condition_matrix with numeric values based on condition
for (i in 1:length(condition_levels)) {
  condition_matrix[, i] <- ifelse(sample_traits$Condition == condition_levels[i], 1, 0)
}

# Calculate module-trait correlations for each condition
moduleTraitCor <- cor(MEs, condition_matrix, use = "pairwise.complete.obs")
moduleTraitPvalue <- corPvalueStudent(moduleTraitCor, nSamples = nrow(exprData_variable_1))

# Open a new plotting device and save the plot to PNG
png("heatmap_conditions_separate.png", width = 1300, height = 1000)

# Create a text matrix to display correlation and p-value values
textMatrix <- paste(
  signif(moduleTraitCor, 2), 
  "\n(", 
  signif(moduleTraitPvalue, 1), 
  ")", 
  sep = ""
)

# Plot the heatmap for different conditions as separate columns
labeledHeatmap(
  Matrix = moduleTraitCor,
  xLabels = condition_levels,  # Separate columns for each condition
  yLabels = names(MEs),
  ySymbols = names(MEs),
  colorLabels = FALSE,
  colors = blueWhiteRed(50),
  textMatrix = textMatrix,
  setStdMargins = TRUE,
  yLabelsPosition = "left",
  cex.text = 1,
  cex.lab = 1.1,     # Column label size
  cex.lab.y = 1.2,   # Row label size
  cex.lab.x = 1.5,   # Axis label size
  cex.legend = 1.2,  # Legend size
  zlim = c(-1, 1)
)

# Close the PNG device
dev.off()

# Open a PDF file and re-plot with smaller dimensions
pdf("heatmap_conditions_separate.pdf", width = 13, height = 12)

labeledHeatmap(
  Matrix = moduleTraitCor,
  xLabels = condition_levels,  # Separate columns for each condition
  yLabels = names(MEs),
  ySymbols = names(MEs),
  colorLabels = FALSE,
  colors = blueWhiteRed(50),
  textMatrix = textMatrix,
  setStdMargins = TRUE,
  yLabelsPosition = "left",
  cex.text = 1,
  cex.lab = 1.1,    
  cex.lab.y = 1.2,   
  cex.lab.x = 1.5,   
  cex.legend = 1.2,  
  zlim = c(-1, 1)
)

# Close the PDF device
dev.off()

# Ensure that both sample_traits and MEs have the same row names
MEs <- MEs[match(rownames(sample_traits), rownames(MEs)), ]

# Map severity numeric values to their respective categories
sample_traits$Severity <- factor(sample_traits$Severity, 
                                 levels = c(1, 4, 5, 3, 2), 
                                 labels = c("Convalescent", "Moderate", "Severe", "ICU", "Healthy"))

# Map gender numeric values to 'Female' and 'Male'
sample_traits$Gender <- factor(sample_traits$Gender, 
                               levels = c(1, 2), 
                               labels = c("Female", "Male"))

# Subset the data for severity and gender columns
severity_data <- sample_traits[, "Severity", drop = FALSE]
gender_data <- sample_traits[, "Gender", drop = FALSE]

# Ensure that the row names of severity and gender data match MEs
severity_data <- severity_data[rownames(MEs), ]
gender_data <- gender_data[rownames(MEs), ]

# Create separate matrices for each severity and gender group
severity_levels <- levels(sample_traits$Severity)
gender_levels <- levels(sample_traits$Gender)

# Create matrices to store severity and gender information
severity_matrix <- as.data.frame(matrix(NA, nrow = nrow(MEs), ncol = length(severity_levels)))
gender_matrix <- as.data.frame(matrix(NA, nrow = nrow(MEs), ncol = length(gender_levels)))

# Assign column names based on severity and gender levels
colnames(severity_matrix) <- severity_levels
colnames(gender_matrix) <- gender_levels

# Fill the severity matrix with 1 for matching condition and 0 otherwise
for (i in 1:length(severity_levels)) {
  severity_matrix[, i] <- ifelse(sample_traits$Severity == severity_levels[i], 1, 0)
}

# Fill the gender matrix with 1 for matching condition and 0 otherwise
for (i in 1:length(gender_levels)) {
  gender_matrix[, i] <- ifelse(sample_traits$Gender == gender_levels[i], 1, 0)
}

# Calculate module-trait correlations for severity and gender
moduleTraitCor_severity <- cor(MEs, severity_matrix, use = "pairwise.complete.obs")
moduleTraitPvalue_severity <- corPvalueStudent(moduleTraitCor_severity, nSamples = nrow(exprData_variable_1))

moduleTraitCor_gender <- cor(MEs, gender_matrix, use = "pairwise.complete.obs")
moduleTraitPvalue_gender <- corPvalueStudent(moduleTraitCor_gender, nSamples = nrow(exprData_variable_1))

# Open a new plotting device and save the plot to PNG for severity
png("heatmap_severity_separate.png", width = 1300, height = 1000)

# Create a text matrix to display correlation and p-value values for severity
textMatrix_severity <- paste(
  signif(moduleTraitCor_severity, 2), 
  "\n(", 
  signif(moduleTraitPvalue_severity, 1), 
  ")", 
  sep = ""
)

# Plot the heatmap for severity as separate columns
labeledHeatmap(
  Matrix = moduleTraitCor_severity,
  xLabels = severity_levels,  # Separate columns for each severity
  yLabels = names(MEs),
  ySymbols = names(MEs),
  colorLabels = FALSE,
  colors = blueWhiteRed(50),
  textMatrix = textMatrix_severity,
  setStdMargins = TRUE,
  yLabelsPosition = "left",
  cex.text = 1,
  cex.lab = 1.1,     # Column label size
  cex.lab.y = 1.2,   # Row label size
  cex.lab.x = 1.5,   # Axis label size
  cex.legend = 1.2,  # Legend size
  zlim = c(-1, 1)
)

# Close the PNG device
dev.off()

# Open a PDF file and re-plot with smaller dimensions for severity
pdf("heatmap_severity_separate.pdf", width = 13, height = 12)

labeledHeatmap(
  Matrix = moduleTraitCor_severity,
  xLabels = severity_levels,  # Separate columns for each severity
  yLabels = names(MEs),
  ySymbols = names(MEs),
  colorLabels = FALSE,
  colors = blueWhiteRed(50),
  textMatrix = textMatrix_severity,
  setStdMargins = TRUE,
  yLabelsPosition = "left",
  cex.text = 1,
  cex.lab = 1.1,    
  cex.lab.y = 1.2,   
  cex.lab.x = 1.5,   
  cex.legend = 1.2,  
  zlim = c(-1, 1)
)

# Close the PDF device
dev.off()

# Open a new plotting device and save the plot to PNG for gender
png("heatmap_gender_separate.png", width = 1300, height = 1000)

# Create a text matrix to display correlation and p-value values for gender
textMatrix_gender <- paste(
  signif(moduleTraitCor_gender, 2), 
  "\n(", 
  signif(moduleTraitPvalue_gender, 1), 
  ")", 
  sep = ""
)

# Plot the heatmap for gender as separate columns
labeledHeatmap(
  Matrix = moduleTraitCor_gender,
  xLabels = gender_levels,  # Separate columns for each gender
  yLabels = names(MEs),
  ySymbols = names(MEs),
  colorLabels = FALSE,
  colors = blueWhiteRed(50),
  textMatrix = textMatrix_gender,
  setStdMargins = TRUE,
  yLabelsPosition = "left",
  cex.text = 1,
  cex.lab = 1.1,     # Column label size
  cex.lab.y = 1.2,   # Row label size
  cex.lab.x = 1.5,   # Axis label size
  cex.legend = 1.2,  # Legend size
  zlim = c(-1, 1)
)

# Close the PNG device
dev.off()

# Open a PDF file and re-plot with smaller dimensions for gender
pdf("heatmap_gender_separate.pdf", width = 13, height = 12)

labeledHeatmap(
  Matrix = moduleTraitCor_gender,
  xLabels = gender_levels,  # Separate columns for each gender
  yLabels = names(MEs),
  ySymbols = names(MEs),
  colorLabels = FALSE,
  colors = blueWhiteRed(50),
  textMatrix = textMatrix_gender,
  setStdMargins = TRUE,
  yLabelsPosition = "left",
  cex.text = 1,
  cex.lab = 1.1,    
  cex.lab.y = 1.2,   
  cex.lab.x = 1.5,   
  cex.legend = 1.2,  
  zlim = c(-1, 1)
)

# Close the PDF device
dev.off()

# Plot the distribution of module membership for a specific module
module_num <- 1  # Replace with the module number of interest (e.g., module 1)
membership_values <- moduleMembership[, module_num]

# Increase margins and plot window size
par(mar = c(5, 5, 4, 2) + 0.1)  # Increase margins if necessary
dev.new(width = 8, height = 6) 
# Specify the output file (PNG)
png("module_membership_distribution.png", width = 800, height = 600)


# Histogram to visualize the distribution of module membership
hist(membership_values, breaks = 30, main = paste("Distribution of Module Membership for Module", module_num),
     xlab = "Module Membership", col = "lightblue", border = "black")

# Clear the current plot window
dev.off() 

threshold <- 0.90
# Calculate module membership for each gene
moduleMembership <- as.data.frame(cor(exprData_variable_1, MEs, use = "p"))
# Select hub genes for MEyellow
hubGenes_yellow <- rownames(moduleMembership)[moduleMembership$MEyellow > 0.93]

# Print the hub genes
cat("Hub genes in the yellow module:\n")
print(hubGenes_yellow)

# Select hub genes for MEred
hubGenes_red <- rownames(moduleMembership)[moduleMembership$MEred > threshold]

cat("\nHub genes in the red module:\n")
print(hubGenes_red)

# Save the top 10 genes into CSV files
write.csv(data.frame(Gene = hubGenes_yellow), "Top10_HubGenes_Yellow.csv", row.names = FALSE)
write.csv(data.frame(Gene = hubGenes_red), "Top10_HubGenes_Red.csv", row.names = FALSE)

threshold <- 0.90

# List of modules to extract hub genes from
modules_of_interest <- c("MEyellow", "MEgrey", "MEbrown", "MEtan", "MEgreenyellow", "MElightcyan", "MEblue", "MEpurple")

# Calculate module membership for each gene
moduleMembership <- as.data.frame(cor(exprData_variable_1, MEs, use = "p"))

# Initialize a vector to store all hub genes with their module membership values
allHubGenes <- data.frame(Gene = character(), Module = character(), Membership = numeric(), stringsAsFactors = FALSE)

# Loop through each module and extract hub genes
for (module in modules_of_interest) {
  # Select hub genes for the current module based on the threshold
  hubGenes <- rownames(moduleMembership)[moduleMembership[[module]] > threshold]
  
  # Check if there are hub genes for this module
  if (length(hubGenes) > 0) {
    # Add hub genes to the list along with their module and membership values
    temp <- data.frame(Gene = hubGenes, Module = module, Membership = moduleMembership[hubGenes, module], stringsAsFactors = FALSE)
    allHubGenes <- rbind(allHubGenes, temp)
  }
}

# Sort all hub genes by their membership value (in descending order)
sortedHubGenes <- allHubGenes[order(-allHubGenes$Membership), ]

# Select the top 20 hub genes overall
top_20_genes <- head(sortedHubGenes$Gene, 20)

# Print the top 20 hub genes
cat("Top 20 overall hub genes:\n")
print(top_20_genes)

# Save the top 20 hub genes into a CSV file
write.csv(data.frame(Gene = top_20_genes), "Top20_Overall_HubGenes.csv", row.names = FALSE)

