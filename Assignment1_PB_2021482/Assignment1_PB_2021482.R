#   Differential expression analysis with limma
library(GEOquery)
library(limma)
library(umap)

# load series and platform data from GEO

gset <- getGEO("GSE2034", GSEMatrix =TRUE, AnnotGPL=TRUE)
if (length(gset) > 1) idx <- grep("GPL96", attr(gset, "names")) else idx <- 1
gset <- gset[[idx]]

# make proper column names to match toptable 
fvarLabels(gset) <- make.names(fvarLabels(gset))

# group membership for all samples
gsms <- "undefined"
sml <- strsplit(gsms, split="")[[1]]

# log2 transformation
ex <- exprs(gset)
qx <- as.numeric(quantile(ex, c(0., 0.25, 0.5, 0.75, 0.99, 1.0), na.rm=T))
LogC <- (qx[5] > 100) ||
  (qx[6]-qx[1] > 50 && qx[2] > 0)
if (LogC) { ex[which(ex <= 0)] <- NaN
exprs(gset) <- log2(ex) }

# assign samples to groups and set up design matrix
gs <- factor(sml)
groups <- make.names(c("undefined"))
levels(gs) <- groups
gset$group <- gs
design <- model.matrix(~group + 0, gset)
colnames(design) <- levels(gs)

fit <- lmFit(gset, design)  # fit linear model

# set up contrasts of interest and recalculate model coefficients
cts <- paste(groups, c(tail(groups, -1), head(groups, 1)), sep="-")
cont.matrix <- makeContrasts(contrasts=cts, levels=design)
fit2 <- contrasts.fit(fit, cont.matrix)

# compute statistics and table of top significant genes
fit2 <- eBayes(fit2, 0.01)
tT <- topTable(fit2, adjust="fdr", sort.by="B", number=250)

tT <- subset(tT, select=c("ID","adj.P.Val","P.Value","F","GB_ACC","SPOT_ID","Gene.Symbol","Gene.symbol","Gene.title"))
write.table(tT, file=stdout(), row.names=F, sep="\t")

# Visualize and quality control test results.
# Build histogram of P-values for all genes. Normal test
# assumption is that most genes are not differentially expressed.
tT2 <- topTable(fit2, adjust="fdr", sort.by="B", number=Inf)
hist(tT2$adj.P.Val, col = "grey", border = "white", xlab = "P-adj",
     ylab = "Number of genes", main = "P-adj value distribution")

# summarize test results as "up", "down" or "not expressed"
dT <- decideTests(fit2, adjust.method="fdr", p.value=0.05)

# create Q-Q plot for t-statistic
t.good <- which(!is.na(fit2$F)) # filter out bad probes
qqt(fit2$t[t.good], fit2$df.total[t.good], main="Moderated t statistic")

# volcano plot (log P-value vs log fold change)
colnames(fit2) # list contrast names
ct <- 1        # choose contrast of interest
volcanoplot(fit2, coef=ct, main=colnames(fit2)[ct], pch=20,
            highlight=length(which(dT[,ct]!=0)), names=rep('+', nrow(fit2)))

# MD plot (log fold change vs mean log expression)
# highlight statistically significant (p-adj < 0.05) probes
plotMD(fit2, column=ct, status=dT[,ct], legend=F, pch=20, cex=1)
abline(h=0)

################################################################
# General expression data analysis
ex <- exprs(gset)

# box-and-whisker plot
ord <- order(gs)  # order samples by group
palette(c("#1B9E77", "#7570B3", "#E7298A", "#E6AB02", "#D95F02",
          "#66A61E", "#A6761D", "#B32424", "#B324B3", "#666666"))
par(mar=c(7,4,2,1))
title <- paste ("GSE2034", "/", annotation(gset), sep ="")
boxplot(ex[,ord], boxwex=0.6, notch=T, main=title, outline=FALSE, las=2, col=gs[ord])
legend("topleft", groups, fill=palette(), bty="n")

# expression value distribution
par(mar=c(4,4,2,1))
title <- paste ("GSE2034", "/", annotation(gset), " value distribution", sep ="")
plotDensities(ex, group=gs, main=title, legend ="topright")

# UMAP plot (dimensionality reduction)
ex <- na.omit(ex) # eliminate rows with NAs
ex <- ex[!duplicated(ex), ]  # remove duplicates
ump <- umap(t(ex), n_neighbors = 4, random_state = 123)
par(mar=c(3,3,2,6), xpd=TRUE)
plot(ump$layout, main="UMAP plot, nbrs=4", xlab="", ylab="", col=gs, pch=20, cex=1.5)
legend("topright", inset=c(-0.15,0), legend=levels(gs), pch=20,
       col=1:nlevels(gs), title="Group", pt.cex=1.5)
library("maptools")  # point labels without overlaps
pointLabel(ump$layout, labels = rownames(ump$layout), method="SANN", cex=0.6)

# mean-variance trend, helps to see if precision weights are needed
plotSA(fit2, main="Mean variance trend, GSE2034")



#PLease run from here above menton code is taaken as a reference GEO
#-----------------------------------Ques1--------------------------------------

# Install and load GEOquery package
install.packages("BiocManager")
BiocManager::install("GEOquery")
install.packages("GEOquery")
library(GEOquery)
library(limma)
library(ggplot2)

# Download microarray data from GEO
gse <- getGEO("GSE2034")

gse
# Extract expression matrix
exprs <- exprs(gse[[1]])



#-----------------------------------Ques2--------------------------------------
# Extract expression matrix
exprs_matrix <- exprs(gse[[1]])

# Perform preprocessing and EDA
exprs_norm <- normalizeQuantiles(exprs_matrix)
exprs_norm100 <- exprs_norm[1:100,1:100]
#boxplot(exprs_norm)

boxplot(exprs_norm100)
pdata <- pData(gse[[1]])
#pdata <- pdata[1:100,]
fdata <- fData(gse[[1]])
#fdata <- fdata[1:100,]

# List data attributes
# view data attributes of exprs_norm
attributes(exprs_norm)
attributes(pdata)
attributes(fdata)

#-----------------------------------Ques3--------------------------------------
# Perform log2 transformation
exprs_log <- log2(exprs_matrix)

# Perform log10 transformation
exprs_log10 <- log10(exprs_matrix)


#-----------------------------------Ques4--------------------------------------

# Split the data into two groups
# Split the data into two groups
# Calculate variance for each gene
gene_var <- apply(exprs_log, 1, var, na.rm = TRUE)

# Filter out genes with low variance
exprs_filtered <- exprs_log[gene_var > 0.1, ]

# Split the data into two groups
group1 <- pdata$characteristics_ch1.6 == "Clinical info: Final microarray diagnosis: ABC DLBCL"
group2 <- pdata$characteristics_ch1.6 == "Clinical info: Final microarray diagnosis: GCB DLBCL"
group1_data <- exprs_filtered[, group1]
group2_data <- exprs_filtered[, group2]

sum(group1)
sum(group2)

# Calculate variance of each gene across all samples
gene_var <- apply(exprs_log, 1, var, na.rm = TRUE)


# Remove rows with missing values
exprs_filtered <- exprs[complete.cases(exprs),]

# Perform t-test
t_test_result <- apply(exprs_filtered, 1, function(x) t.test(x[group1], x[group2], na.rm = TRUE))

# Filter out genes with low variance
#exprs_filtered <- exprs_log[gene_var > 0.1, ]

# Extract the log2 fold change for each gene
log_fc <- apply(exprs_filtered, 1, function(x) mean(x[group1]) - mean(x[group2]))

# Extract the p-value for each gene
p_val <- sapply(t_test_result, function(x) x$p.value)

# Correct the p-values for multiple comparisons using Holm's method
adj_p_val <- p.adjust(p_val, method = "holm")

# Create a data frame with the results
t_results <- data.frame(Gene = rownames(exprs_filtered), Log2FC = log_fc, P.Value = p_val, Adj.P.Value = adj_p_val)

# Plot the results using a volcano plot
ggplot(t_results, aes(x = Log2FC, y = -log10(P.Value))) +
  geom_point(aes(color = Adj.P.Value < 0.05), alpha = 0.5) +
  scale_color_manual(values = c("grey50", "red")) +
  geom_vline(xintercept = c(-1, 1), linetype = "dashed") +
  geom_hline(yintercept = -log10(0.05), linetype = "dashed") +
  labs(title = "Volcano plot of differential expression analysis",
       x = "Log2 Fold Change",
       y = "-log10(P-value)",
       color = "Adjusted P-value < 0.05") +
  theme_bw()


#-----------------------------------Ques5--------------------------------------
pdata <- pData(gse[[1]])
# Define the design matrix
design <- model.matrix(~ pdata$characteristics_ch1.6, data = pdata)

# Perform differential expression analysis
fit <- lmFit(exprs, design)
fit <- eBayes(fit)
results <- topTable(fit, coef=2, adjust="fdr", number=Inf)
top_table <- topTable(fit, coef=2, adjust="fdr", number=Inf)

# Draw a volcano plot
ggplot(results, aes(x = logFC, y = -log10(adj.P.Val))) +
  geom_point(aes(color = ifelse(abs(logFC) > 1 & adj.P.Val < 0.05, "red", "black"))) +
  scale_color_manual(values = c("black", "red")) +
  geom_hline(yintercept = -log10(0.05), linetype = "dashed") +
  labs(x = "Log Fold Change", y = "-log10(Adj. P-Value)", title = "Volcano Plot")

#-----------------------------------Ques6--------------------------------------

#Choose a significant cutoff based on log(FC) and p-values and justify why you chose those values as the cutoff.

# Define log fold change cutoff
lfc_cutoff <- 1

# Define adjusted p-value cutoff
pval_cutoff <- 0.05

# Identify significantly differentially expressed genes
sig_genes <- subset(top_table, abs(logFC) > lfc_cutoff & adj.P.Val < pval_cutoff)

# Get number of significantly differentially expressed genes
num_sig_genes <- nrow(sig_genes)

# Print summary information
cat(sprintf("%d genes were found to be significantly differentially expressed with a log fold change cutoff of %s and an adjusted p-value cutoff of %s.", num_sig_genes, lfc_cutoff, pval_cutoff))




#-----------------------------------Ques7--------------------------------------


BiocManager::install("clusterProfiler")
BiocManager::install("org.Hs.eg.db")
BiocManager::install("biomaRt")
# Load the required packages
library(biomaRt)
library(clusterProfiler)
library(org.Hs.eg.db)  # or any other appropriate annotation package for your organism

# Set the gene universe as all genes in the dataset
universe <- rownames(exprs_filtered)

# Set the gene set as the differentially expressed genes with adjusted p-value < 0.05
gene_set <- t_results$Gene[t_results$Adj.P.Value < 0.05]

# Set the universe of genes
universe <- names(org.Hs.eg.db)


library(GEOquery)


# Extract the expression matrix
expr_mat <- exprs(gse[[1]])

# Load required package
library(biomaRt)

# Specify the mart to use (e.g. for human)
mart <- useMart("ensembl", dataset = "hsapiens_gene_ensembl")
gene_list = rownames(expr_mat)


# Get the gene symbols from probe IDs
gene_list <- getBM(
  attributes = c("entrezgene_id", "affy_hg_u133_plus_2"),
  filters = "affy_hg_u133_plus_2",
  values = gene_list,
  mart = mart
)$entrezgene_id

# Run GO enrichment analysis with gene symbols
ego <- enrichGO(gene = gene_list,
                OrgDb = org.Hs.eg.db,
                keyType = "ENTREZID",
                ont = "BP",
                pvalueCutoff = 0.05,
                qvalueCutoff = 0.1,
                readable = TRUE)

# Show the enriched GO terms
ego@result





#-----------------------------------Ques8--------------------------------------
library(lattice)

# visulalisation

dotplot(ego)

barplot(ego)
# Create a data frame with GO term names and p-values
ego_df <- data.frame(GO.Term = ego@result$ID,
                     P.Value = ego@result$pvalue,
                     stringsAsFactors = FALSE)

# Create a dot plot of the p-values

dotplot(P.Value ~ GO.Term, data = ego_df, showCategory = 15)


barplot(P.Value ~ GO.Term, data = ego_df, showCategory = 15)

barplot(ego@result$Count[1:10], 
        horiz = TRUE, 
        las = 1, 
        main = "Enriched GO terms", 
        xlab = "Number of genes", 
        ylab = "GO term",
        xlim = c(0, max(ego@result$Count) + 10))

# import package
library("ggplot2")

# Extract the enriched GO terms and associated statistics from the enrichResult object
go_terms <- as.data.frame(ego@result[, c("Description", "GeneRatio", "pvalue", "qvalue")])
colnames(go_terms) <- c("GO Term", "Gene Ratio", "p-value", "qvalue")


install.packages("dplyr")
library(dplyr)
names(go_terms)[names(go_terms) == "p-value"] <- "pvalue"
go_terms_top <- go_terms %>%top_n(10, qvalue) 

# Create a dot plot
#ggplot(data = go_terms[1:20,], aes(x = qvalue, y = `GO Term`, size = as.numeric(`Gene Ratio`), color = -log10(`pvalue`))) +
  #geom_point() +
  #scale_color_gradient(low = "blue", high = "red") +
  #scale_size(range = c(5, 15)) +
  #labs(title = "Enriched GO Terms", x = "-log10(pvalue)", y = "GO Term", size = "Gene Ratio", color = "-log10(pvalue)") +
  #theme_minimal() +
  #theme(axis.text.y = element_text(size = 8))



#-----------------------------------Ques9--------------------------------------

if (!require("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

BiocManager::install("pathview")
BiocManager::install("clusterProfiler")

library(pathview)
library(clusterProfiler)



library(org.Hs.eg.db)

ensembl <- useMart("ensembl", dataset = "hsapiens_gene_ensembl")
# Convert the Affymetrix probe IDs to Ensembl IDs
ensembl_ids <- getBM(attributes = c("ensembl_gene_id"), 
                     filters = "affy_hg_u133_plus_2", 
                     values = gene_set, 
                     mart = ensembl)

# Convert the Ensembl IDs to gene symbols
gene_symbols <- getBM(attributes = c("ensembl_gene_id", "external_gene_name"), 
                      filters = "ensembl_gene_id", 
                      values = ensembl_ids$ensembl_gene_id, 
                      mart = ensembl)

# Reinstall the clusterProfiler package
BiocManager::install("clusterProfiler", dependencies = TRUE)

# Load the package
library(clusterProfiler)

BiocManager::install("gage")
library(KEGGgraph)
library(kegg)
library(gage)


# get the top pathways by adjusted p-value
top_pathways <- head(ego@result, n=10)

# extract the pathway names
pathway_names <- top_pathways$Description

# print the pathway names
cat("Top pathways:\n")
for (i in 1:length(pathway_names)) {
  cat(i, pathway_names[i], "\n")
}