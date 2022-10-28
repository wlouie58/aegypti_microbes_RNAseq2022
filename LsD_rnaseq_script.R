### RNA-Seq workflow - WL

### Pre-processing raw reads done on terminal, using bowtie2 and salmon

path <- "/Users/williamlouie/Dropbox/My Mac (Williams-MacBook-Pro.local)/Desktop/LsD_RNA-Seq"



library(tximport)
library(readr)
library(tximportData)

# First, identify the files with read counts, determined by salmon
# Read the folder called "salmon_quant", then name all subsequent quant files the name of their parent folder (sample ID)
dirs <- list.files("salmon_quant/")
quant_files <- list.files("salmon_quant/",pattern="quant.sf",recursive = TRUE,full.names = TRUE)
names(quant_files) <- dirs
quant_files

# Now we can inspect our count files from salmon
quants <- read_tsv(quant_files[1])
head(quants)

# The TPM are giving in the table, but we can demonstrate the calculation in a few steps
# divide the number of reads for each transcript by it’s length (reads per kilobase - RPK)
# sum the RPK values and divide by 1 million to get a scaling factor
# divide the RPK values by the scaling factor to get the TPM
rpk <- quants$NumReads / quants$EffectiveLength
scale_factor <- sum(rpk) / 1e6
tpm <- rpk / scale_factor



# Define transcript mapping
# need to input your .gtf file that has gene annotations
# *.gtf file didn't work... used the .gff file instead
gtf_file <- "VectorBase-54_AaegyptiLVP_AGWG.gff"
file.exists(gtf_file)
# Can also download .gtf file directly from the ENSEMBL, if needed
gtf_file <- "Aedes_gtf"
download.file("https://ftp.ncbi.nlm.nih.gov/genomes/all/annotation_releases/7159/101/GCF_002204515.2_AaegL5.0/GCF_002204515.2_AaegL5.0_genomic.gtf.gz", destfile = gtf_file)
txdb <- makeTxDbFromGFF(gtf_file, format="gff3", organism = "Aedes aegypti")
keytypes(txdb)
columns(txdb)
transcripts(txdb)
gene_names <- transcriptsBy(txdb, "gene")
length(gene_names)
exon_names <- transcriptsBy(txdb, "exon")
length(exon_names)
# To make gene names correspond to every transcript in the database
k <- keys(txdb, keytype="TXNAME")
cols <- c("GENEID", "CDSNAME", "CDSID", "EXONNAME", "EXONID")
tx_map <- select(txdb, keys = k, columns = "CDSNAME", keytype = "TXNAME")
head(tx_map)
# 111 exons couldn't be linked and were removed
# A lot of GENEIDs were gone, so I'm using a different method



# Using a .txt file from Biomart (ENSEMBL) as the annotation
# Create transcript database
library(GenomicFeatures)
library(readr)
library(dplyr)

tx2gene <- read_csv(file.path("anno_mart.txt"))
tx2gene

# Reorder the annotation coloumns so that transcript ID matches with the transcript ID in the quant files
# i.e. move to first column
tx2gene_order = tx2gene %>% select("Transcript stable ID", "Gene stable ID", "Protein stable ID", "Gene name", "GO term accession", "GO term name", "Gene description", "KEGG Pathway and Enzyme ID")
head(tx2gene_order)
colnames(tx2gene_order) <- c("Transcript", "Gene", "Protein_ID", "Gene_name", "GO_accession", "GO_term", "Gene_description", "KEGG_ID")
head(tx2gene_order)
# Use tximport
txi <- tximport(quant_files, type="salmon",tx2gene = tx2gene_order)
head(txi)
# These are your accessors for txi
names(txi)
# Intitial look at read counts
head(txi$counts)

# Extract transcript-level TPM values in the quant files, summarised to the gene-level
tpm <- txi$abundance
write.csv(tpm, file="tpm_values.csv",quote=FALSE)
# Also transform to fragments per kilobase per million mapped fragments (fpkm)
fpm <- fpm(dds)
write.csv(fpm, file="fpm_values.csv",quote=FALSE)
fpkm <- fpkm(dds)
write.csv(fpkm, file="fpkm_values.csv",quote=FALSE)




################################################################
# Now we can analyze the counts with edgeR
# Now we can analyze the counts with edgeR
# Now we can analyze the counts with edgeR

# Tutorial: https://ucdavis-bioinformatics-training.github.io/2018-June-RNA-Seq-Workshop/thursday/DE.html

library(edgeR)
library(Rgraphviz)
library(gplots)

# Import metadata
metadata <- read.table("LsD_metadata.txt", sep="\t", header=T, row.names=1)
View(metadata)
rownames(metadata) <- metadata$sample_id

# Counts table
counts <- txi$counts
head(counts)




############################################################################
# Now we can analyze the counts with limma + voom
# Now we can analyze the counts with limma + voom
# Now we can analyze the counts with limma + voom

library(limma)

# Create DGElist object
d0 <- DGEList(counts)
# Calculate normalization factors
d0 <- calcNormFactors(d0)
d0$samples

# Combine with metadata file
# NOTE: MAKE SURE THE SAMPLE ORDERS MATCH!

### Filter low-expressed genes
### Filter low-expressed genes
### Filter low-expressed genes
# at least 1cpm, at least 2 samples
drop <- which(apply(cpm(d0), 1, max) < 2)
d <- d0[-drop,] 
dim(d)
# Log-transform counts to cpm
logcpm <- cpm(d, prior.count=2, log=TRUE)
write.table(logcpm,"Normalized_counts.txt",sep="\t",quote=F)
# Store name identifiers for conditions
# take names for antibiotics variable and bloodmeal variable
# Take sample names
snames <- colnames(counts)
snames
antibiotics <- substr(snames, 1, 3) 
antibiotics
bloodmeal <- substr(snames, nchar(snames) - 1, nchar(snames) - 1)
bloodmeal
# Create groups that combines these variables
group <- interaction(antibiotics, bloodmeal)
group
# Plot by multidimensional scaling
plotMDS(d, col = as.numeric(group))

# Voom transformation and calculating variance weights
mm <- model.matrix(~ 0 + group)
# Let's voom
y <- voom(d, mm, plot = T)

# Looks good!
# What is voom doing?
# Counts are transformed to log2 counts per million reads (CPM), where “per million reads” is defined based on the normalization factors we calculated earlier
# A linear model is fitted to the log2 CPM for each gene, and the residuals are calculated
# A smoothed curve is fitted to the sqrt(residual standard deviation) by average expression (see red line in plot above)
# The smoothed curve is used to obtain weights for each gene and sample that are passed into limma along with the log2 CPMs


### Fit linear models in limma using weighted least squares for each gene
### Fit linear models in limma using weighted least squares for each gene
### Fit linear models in limma using weighted least squares for each gene
fit <- lmFit(y, mm)
head(coef(fit))

# Comparisons between groups (log fold-changes) are obtained as contrasts of these fitted linear models

#################

# Choose groups to compare
## Here we compare groupAbA.b to groupCtr.b
contr <- makeContrasts(groupAbA.b - groupCtr.b, levels = colnames(coef(fit)))
contr
# Estimate contrast for each gene
tmp <- contrasts.fit(fit, contr)
# Use Empirical Bayes smoothing of standard errors
tmp <- eBayes(tmp)

# WHICH GENES ARE DIFFERENTIALLY EXPRESSED?
# expressed as log2fold AbA.b / Ctr.b
top.table <- topTable(tmp, sort.by = "P", n = Inf)
head(top.table, 10)
# How many DE genes are there?
length(which(top.table$adj.P.Val < 0.05))
# How about including multiple test corrections?
# How many DE genes are there, after correcting for false discovery rate?
top.table <- topTable(tmp, adjust.method = "BH", sort.by = "P", n = Inf)
length(which(top.table$adj.P.Val < 0.05))
# Write to a top genes table
top.table$Gene <- rownames(top.table)
top.table <- top.table[,c("Gene", names(top.table)[1:6])]
top.table <- data.frame(top.table, tx2gene_order[match(top.table$Gene, tx2gene_order$Gene),], logcpm[match(top.table$Gene,rownames(logcpm)),])
head(top.table)
write.table(top.table, file = "AbA vs Ctr - before.txt", row.names = F, sep = "\t", quote = F)


plotMA(tmp, ylim=c(-5,5))




############### Repeat these lines for different groups
## NOTE: we can also use
# mm <- model.matrix(~x*y)
# for the model instead of making a group interaction first, as these do the same thing
# basically, it incorporates group x + group y + interaction effect of group x & y

## Here we repeat, but compare groupAbP.b to groupCtr.b
contr <- makeContrasts(groupAbP.b - groupCtr.b, levels = colnames(coef(fit)))
contr
tmp <- contrasts.fit(fit, contr)
tmp <- eBayes(tmp)
# WHICH GENES ARE DIFFERENTIALLY EXPRESSED?
# expressed as log2fold AbP.b / Ctr.b
top.table <- topTable(tmp, sort.by = "P", n = Inf)
head(top.table, 10)
length(which(top.table$adj.P.Val < 0.05))
top.table <- topTable(tmp, adjust.method = "BH", sort.by = "P", n = Inf)
length(which(top.table$adj.P.Val < 0.05))
# Write to a top genes table
top.table$Gene <- rownames(top.table)
top.table <- top.table[,c("Gene", names(top.table)[1:6])]
top.table <- data.frame(top.table, tx2gene_order[match(top.table$Gene, tx2gene_order$Gene),], logcpm[match(top.table$Gene,rownames(logcpm)),])
head(top.table)
write.table(top.table, file = "AbP vs Ctr - before.txt", row.names = F, sep = "\t", quote = F)

## Repeat, but compare groupAbA.b to groupAbP.b
contr <- makeContrasts(groupAbA.b - groupAbP.b, levels = colnames(coef(fit)))
contr
tmp <- contrasts.fit(fit, contr)
tmp <- eBayes(tmp)
# expressed as log2fold AbA.b / AbP.b
top.table <- topTable(tmp, sort.by = "P", n = Inf)
head(top.table, 10)
length(which(top.table$adj.P.Val < 0.05))
top.table$Gene <- rownames(top.table)
top.table <- top.table[,c("Gene", names(top.table)[1:6])]
top.table <- data.frame(top.table, tx2gene_order[match(top.table$Gene, tx2gene_order$Gene),], logcpm[match(top.table$Gene,rownames(logcpm)),])
head(top.table)
write.table(top.table, file = "AbA vs AbP - before.txt", row.names = F, sep = "\t", quote = F)

## Repeat, but compare groupCtr.z to groupCtr.m
contr <- makeContrasts(groupCtr.z - groupCtr.m, levels = colnames(coef(fit)))
contr
tmp <- contrasts.fit(fit, contr)
tmp <- eBayes(tmp)
# expressed as log2fold AbP.m / AbP.m
top.table <- topTable(tmp, sort.by = "P", n = Inf)
head(top.table, 10)
length(which(top.table$adj.P.Val < 0.05))
top.table$Gene <- rownames(top.table)
top.table <- top.table[,c("Gene", names(top.table)[1:6])]
top.table <- data.frame(top.table, tx2gene_order[match(top.table$Gene, tx2gene_order$Gene),], logcpm[match(top.table$Gene,rownames(logcpm)),])
head(top.table)
write.table(top.table, file = "ZIKV vs mock - Ctrl.txt", row.names = F, sep = "\t", quote = F)

# Save these variables from metadata
ab <- metadata$Antibiotic_treatment
bm <- metadata$bloodmeal
bmc <- metadata$Bloodmeal_content

# Voom transformation and calculating variance weights
mm <- model.matrix(~ab*bm)
colnames(mm)
# Let's voom
y <- voom(d, mm, plot = T)
# Looks good!
### Fit linear models in limma using weighted least squares for each gene
fit <- lmFit(y, mm)
head(coef(fit))

# Comparisons between groups (log fold-changes) are obtained as contrasts of these fitted linear models
# Choose groups to compare
## Here we compare groupAbA.b to groupCtr.b
contr <- makeContrasts(bloodmealm - bloodmealz, levels = colnames(coef(fit)))
contr
# Estimate contrast for each gene
tmp <- contrasts.fit(fit, contr)
# Use Empirical Bayes smoothing of standard errors
tmp <- eBayes(tmp)
# WHICH GENES ARE DIFFERENTIALLY EXPRESSED?
# expressed as log2fold AbA.b / Ctr.b
top.table <- topTable(tmp, sort.by = "P", n = Inf)
head(top.table, 10)
# How many DE genes are there?
length(which(top.table$adj.P.Val < 0.05))
# How about including multiple test corrections?
# How many DE genes are there, after correcting for false discovery rate?
top.table <- topTable(tmp, adjust.method = "BH", sort.by = "P", n = Inf)
length(which(top.table$adj.P.Val < 0.05))
# Write to a top genes table
top.table$Gene <- rownames(top.table)
top.table <- top.table[,c("Gene", names(top.table)[1:6])]
top.table <- data.frame(top.table, tx2gene_order[match(top.table$Gene, tx2gene_order$Gene),], logcpm[match(top.table$Gene,rownames(logcpm)),])
head(top.table)
write.table(top.table, file = "AbA vs Ctr - before.txt", row.names = F, sep = "\t", quote = F)









########################################################################
### ANALYZING WITH DESEQ2
### ANALYZING WITH DESEQ2
### ANALYZING WITH DESEQ2

library(DESeq2)
library(readr)
library(dplyr)
library(ggplot2)
library(PCAtools)
library(pheatmap)
library(RColorBrewer)
library(apeglm)
library(reshape2)
library(knitr)
library(ggpubr)

metadata <- read.table("LsD_metadata.txt", sep="\t", header=T, row.names=1)
View(metadata)
rownames(metadata) <- metadata$sample_id

# Check that the counts table and metadata match
counts <- txi$counts
head(counts)
mycounts <- as.data.frame(counts)
metadata <- as.data.frame(metadata)
class(mycounts)
class(metadata)
names(mycounts)[-1]
names(mycounts)[-1]==metadata$sample_id
all(names(mycounts)[-1]==metadata$sample_id)



# Make the dds
#   counts = counts table
#   d0 = normalization factors
#   Design: by bloodmeal content * antibiotic treatment
dds <- DESeqDataSetFromTximport(txi, 
                                colData = metadata,
                                design = ~Antibiotic_treatment + Bloodmeal_content + Antibiotic_treatment:Bloodmeal_content)
colData(dds)

#   Filter non-expressed genes (<5 ocunts)
is_expressed <- assay(dds) >= 10
head(is_expressed)
#   Look at the number of genes expressed in x number of samples
hist(rowSums(is_expressed), main="Number of samples a gene is expressed in", xlab="Sample Count")
#   Looks like most genes are expressed in every sample
#   Visualize library size
sum(assay(dds)[,1])
# Can also look by each sample
colSums(assay(dds))
#   Remove genes not expressed in at least 2 samples
keep <- rowSums(assay(dds) >= 10) >= 3
table(keep)
dds <- dds[keep,]

#   Reorder variables so you get the comparisons you want
dds$Bloodmeal_content = factor(dds$Bloodmeal_content, levels = c("before bloodfeed", "blood only", "blood + ZIKV"))
levels(dds$Bloodmeal_content)
dds$Antibiotic_treatment = factor(dds$Antibiotic_treatment, levels = c("AbxP", "AbxA", "Ctrl"))
levels(dds$Antibiotic_treatment)

dds$Antibiotic_treatment <- relevel(dds$Antibiotic_treatment, ref = "Ctrl")
dds$Bloodmeal_content <- relevel(dds$Bloodmeal_content, ref = "blood only")
# Run the analysis
dds <- DESeq(dds)
#   These are the resulting groups
resultsNames(dds)
#   Look at normalization factors
normalizationFactors(dds)





# Access results
# Alpha is FDR, default alpha=0.1
# WE GOT >10K GENES! TRY STRICTER FILTERING FOF LFCTHRESHOLD
# lfcThreshold can be changed from 0 to raise log-fold change

# With ZIKV treatment, what are the DE between AbP and Ctrl
# This looks at AbP vs. Ctrl + interaction effect ==> 63 genes!
res.zika_abp = results(dds, list( c("Antibiotic_treatment_pupa_vs_no_ab","Antibiotic_treatmentpupa.Bloodmeal_contentzikv") ))
ix = which.min(res.zika_abp$padj)
res.zika_abp <- res.zika_abp[order(res.zika_abp$padj),]
kable(res.zika_abp[1:5,-(3:4)])
summary(res.zika_abp)
#   How many are padj <0.05
res05 <- subset(res.zika_abp, padj < 0.05)
sum(res05$padj < 0.05)
res05 <- cbind(rownames(res05), res05)
rownames(res05) <- NULL
colnames(res05) <- c("Gene", "baseMean", "log2FoldChange","lfcSE", "stat", "pvalue", "padj")
res05 <- data.frame(res05, tx2gene_order[match(res05$Gene, tx2gene_order$Gene),], logcpm[match(res05$Gene, rownames(logcpm)),])
head(res05)
# Remove rows with NA
res05 <- res05[!grepl("NA", rownames(res05)),]
# Reorder by fold change
res05 <- res05[order(res05$log2FoldChange),]
head(res05)
# Export csv file
write.csv(as.data.frame(res05), file = "Sig - Zika response, AbP vs Ctrl.csv")

# With ZIKV treatment, what are the DE between AbA and Ctrl
#   This looks at AbA vs. Ctrl + interaction effect ==> 1978 genes!
res.zika_aba = results(dds, list( c("Antibiotic_treatment_adult_vs_no_ab","Antibiotic_treatmentadult.Bloodmeal_contentzikv") ))
summary(res.zika_aba)
#   How many are padj <0.05
res05 <- subset(res.zika_aba, padj < 0.05)
sum(res05$padj < 0.05)
res05 <- cbind(rownames(res05), res05)
rownames(res05) <- NULL
colnames(res05) <- c("Gene", "baseMean", "log2FoldChange","lfcSE", "stat", "pvalue", "padj")
res05 <- data.frame(res05, tx2gene_order[match(res05$Gene, tx2gene_order$Gene),], logcpm[match(res05$Gene, rownames(logcpm)),])
head(res05)
# Remove rows with NA
res05 <- res05[!grepl("NA", rownames(res05)),]
# Reorder by fold change
res05 <- res05[order(res05$log2FoldChange),]
head(res05)
# Export csv file
write.csv(as.data.frame(res05), file = "Sig - Zika response, AbA vs Ctrl.csv")






### PCA Plot
### PCA Plot
### PCA Plot

# Get log2 counts
vsd <- varianceStabilizingTransformation(dds, blind=TRUE)
#   Check distributions of samples using boxplots
boxplot(assay(vsd), xlab="", ylab="Log2 counts per million", las=2, main="Normalised Distributions")
#   Calculate sample distances
sampleDists <- dist(t(assay(vsd)))
sampleDists
sampleDistMatrix <- as.matrix( sampleDists )
rownames(sampleDistMatrix) <- paste(vsd$Bloodmeal_content, vsd$Antibiotic_treatment, sep = " - " )
colnames(sampleDistMatrix) <- NULL
colors <- colorRampPalette( rev(brewer.pal(9, "Blues")) )(255)
pheatmap(sampleDistMatrix,
         clustering_distance_rows = sampleDists,
         clustering_distance_cols = sampleDists,
         col = colors)
# Try out a PCA plot
pca_ab <- plotPCA(vsd, ntop = 500, intgroup = "Antibiotic_treatment")
pca_ab
pca_bmc <- plotPCA(vsd, ntop = 500, intgroup = "Bloodmeal_content")
pca_bmc
#   Can also generate from scratch, using values fed into ggplot
vst$Bloodmeal_content = factor(vst$Bloodmeal_content, levels = c("before bloodfeed", "blood only", "blood + ZIKV"))
levels(vst$Bloodmeal_content)
pca_data <- plotPCA(vsd, intgroup = c("Bloodmeal_content", "Antibiotic_treatment"), returnData = TRUE)
pca_data
percentVar <- round(100 * attr(pca_data, "percentVar"))
pca_all <- ggplot(pca_data, aes(x = PC1, y = PC2, shape = Bloodmeal_content, color = Antibiotic_treatment)) +
            geom_point(size =5) +
            xlab(paste0("PC1: ", percentVar[1], "%")) +
            ylab(paste0("PC2: ", percentVar[2], "%")) +
            coord_fixed() +
            theme(panel.background = element_rect(fill = "white", color = "black"), text = element_text(size=30)) +
            scale_color_manual(values = c( "darkgreen", "blue2", "darkorange3"))
pca_all
ggsave("PCA by Bm content and Abx treatment.pdf", plot = pca_all, width = 12, height = 12, dpi = 600)





### Heatmap
### Heatmap
### Heatmap

# Transform count data using the variance stablilizing transform
vsd
#   Convert the DESeq transformed object to a data frame
ds2vst <- assay(vsd)
ds2vst <- as.data.frame(ds2vst)
ds2vst$Gene <- rownames(ds2vst)
head(ds2vst)

## For AbP vs Ctrl
#   Extract only sig DE gene names
sigGenes <- res05$Gene
#   Filter deseq2VST by the sig genes
ds2vst <- ds2vst[ds2vst$Gene %in% sigGenes,]
#   Filter info from tx2gene_order by the sig genes, will add this to annotation
sigGO <- tx2gene_order[tx2gene_order$Gene %in% sigGenes,]
gene_name <- sigGO$Gene_name
#   Melt (long), overwrite our original data frame with the long format
ds2vst_long <- melt(ds2vst, id.vars=c("Gene"))
head(ds2vst_long)
ds2vst_long <- data.frame(ds2vst_long, sigGO[match(ds2vst_long$Gene, sigGO$Gene),])
# Plot heatmap
#   Annotate by GO term
htmp <- ggplot(ds2vst_long, aes(x=variable, y=Gene, fill=value)) + 
            geom_raster() + 
            facet_grid(rows = vars(GO_term), scales = "free", space = "free_y", switch = "y") +
            theme(axis.text.x = element_text(angle=65, hjust=1), 
        strip.placement = "outside", strip.text.y.left = element_text(angle = 0))
htmp
ggsave("Heatmap - AbP vs Ctrl.pdf", plot = htmp, width = 10, height = 15, dpi = 600)



### Compare AbA and AbP in context of ZIKV
# With ZIKV treatment, what are the DE between AbA and AbP
#   This looks at AbA vs AbP, accounting for ZIKV interaction ==> 175 genes!
res.zika_aba_abp = results(dds, list( c("Antibiotic_treatmentAbxA.Bloodmeal_contentblood...ZIKV","Antibiotic_treatmentAbxP.Bloodmeal_contentblood...ZIKV") ))
#   How many are padj <0.05
res05 <- subset(res.zika_aba_abp, padj < 0.05)
sum(res05$padj < 0.05)
res05 <- cbind(rownames(res05), res05)
rownames(res05) <- res05$Gene
colnames(res05) <- c("Gene", "baseMean", "log2FoldChange","lfcSE", "stat", "pvalue", "padj")
res05 <- data.frame(res05, tx2gene_order[match(res05$Gene, tx2gene_order$Gene),], logcpm[match(res05$Gene, rownames(logcpm)),])
head(res05)
# Remove rows with NA
res05 <- res05[!grepl("NA", rownames(res05)),]
# Reorder by fold change
res05 <- res05[order(res05$log2FoldChange),]
head(res05)
#   Filter only those with magnitude log2fold change higher than 4 (64 times)
res05_top <- subset(res05, abs(log2FoldChange) >= 4)
head(res05_top)
# Export csv file
write.csv(as.data.frame(res05), file = "Sig - Zika response, AbA vs AbP.csv")

# Transform count data using the variance stablilizing transform
vsd
#   Convert the DESeq transformed object to a data frame
ds2vst <- assay(vsd)
ds2vst <- as.data.frame(ds2vst)
ds2vst$Gene <- rownames(ds2vst)
head(ds2vst)

## For AbP vs AbP, ZIKV
# Use this
head(ds2vst)
#   Extract only sig DE gene names (from res.zika_aba_abp)
sigGenes <- res05_top$Gene
#   Filter deseq2VST by the sig genes
ds2vst <- ds2vst[ds2vst$Gene %in% sigGenes,]
head(ds2vst)
#   Average replicates
ds2vst$AbA_before <- rowMeans(ds2vst[1:3], na.rm=TRUE)
ds2vst$AbA_mock <- rowMeans(ds2vst[4:6], na.rm=TRUE)
ds2vst$AbA_zika <- rowMeans(ds2vst[7:9], na.rm=TRUE)
ds2vst$AbP_before <- rowMeans(ds2vst[10:12], na.rm=TRUE)
ds2vst$AbP_mock <- rowMeans(ds2vst[13:15], na.rm=TRUE)
ds2vst$AbP_zika <- rowMeans(ds2vst[16:18], na.rm=TRUE)
head(ds2vst)
ds2vst <- ds2vst[, -c(1:27)]
head(ds2vst)
#   Filter info from tx2gene_order by the sig genes, will add this to annotation
sigGO <- tx2gene_order[tx2gene_order$Gene %in% sigGenes,]
gene_name <- sigGO$Gene_name
#   Melt (long), overwrite our original data frame with the long format
ds2vst_long <- melt(ds2vst, id.vars=c("Gene"))
head(ds2vst_long)
ds2vst_long <- data.frame(ds2vst_long, sigGO[match(ds2vst_long$Gene, sigGO$Gene),])
# Plot heatmap
#   Annotate by GO term
htmp <- ggplot(ds2vst_long, aes(x=variable, y=Gene, fill=value)) + 
  geom_raster() + 
  facet_grid(rows = vars(GO_term), scales = "free", space = "free_y", switch = "y") +
  theme(axis.text.x = element_text(angle=65, hjust=1), 
        strip.placement = "outside", strip.text.y.left = element_text(angle = 0))
htmp
ggsave("Heatmap - AbA vs AbP.pdf", plot = htmp, width = 10, height = 15, dpi = 600)

# Add annotations (l2fc to heatmap)
head(res05_top)
head(ds2vst)
ds2vst_plot <- ds2vst[, -1]
head(ds2vst_plot)
ds2vst_plot <- as.matrix(ds2vst_plot)

annotation_row = data.frame(res05_top$log2FoldChange)
colnames(annotation_row) = "AbA/AbP (Log2-fold)"
rownames(annotation_row) = res05_top$Gene
head(annotation_row)

labs.row <- paste(res05_top$Gene, res05_top$GO_term)

myColor <- list("AbA/AbP (Log2-fold)" = brewer.pal(11, "RdBu"))

phtmp <- pheatmap(ds2vst_plot, cluster_cols = F, cluster_rows = F, color = colorRampPalette(c("black", "green"))(20),
                  annotation_row = annotation_row, annotation_colors = myColor, annotation_names_row = F,
                  labels_row = labs.row, fontsize = 25)
phtmp
ggsave("Heatmap - AbA vs AbP.pdf", plot = phtmp, width = 25, height = 25, dpi = 600)





### Filter by immune genes
### Filter by immune genes
### Filter by immune genes

immune_list <- read.csv("Immune genes list.csv", header=T)
head(immune_list)



# Using res05 from AbA vs Ctrl
res05 <- subset(res.zika_aba, padj < 0.05)
sum(res05$padj < 0.05)
res05 <- cbind(rownames(res05), res05)
rownames(res05) <- NULL
colnames(res05) <- c("Gene", "baseMean", "log2FoldChange","lfcSE", "stat", "pvalue", "padj")
res05 <- data.frame(res05, tx2gene_order[match(res05$Gene, tx2gene_order$Gene),], logcpm[match(res05$Gene, rownames(logcpm)),])
head(res05)
res05 <- res05[!grepl("NA", rownames(res05)),]
res05 <- res05[order(res05$log2FoldChange),]
head(res05)
# Filter out immune genes, if present
immune_gene_list <- immune_list$Gene
immune_genes <- subset(res05, Gene %in% immune_gene_list)
immune_genes
# Paste my custom annotations
immune_genes$name <- immune_list$Name[match(immune_genes$Gene, immune_list$Gene)]
immune_genes$pathway <- immune_list$Pathway[match(immune_genes$Gene, immune_list$Gene)]
immune_genes$description <- immune_list$Description[match(immune_genes$Gene, immune_list$Gene)]

ggplot(immune_genes, aes(x=name, y=log2FoldChange, size = padj)) +
  geom_point(position=position_jitter(width=.1, height=0)) +
  scale_size(trans = 'reverse') +
  facet_grid(~pathway, scales = "free_x", space = "free_x") + 
  theme(panel.background = element_rect(fill = "white", color = "black"),
          text = element_text(size=20), axis.text.x = element_text(angle = 50, hjust = 1))

# Filter ALL immune genes, not just res05
immune.zika_aba <- subset(res.zika_aba, rownames(res.zika_aba) %in% immune_list$Gene)
immune.zika_aba <- cbind(rownames(immune.zika_aba), immune.zika_aba)
rownames(immune.zika_aba) <- NULL
colnames(immune.zika_aba) <- c("Gene", "baseMean", "log2FoldChange","lfcSE", "stat", "pvalue", "padj")
immune.zika_aba <- data.frame(immune.zika_aba, tx2gene_order[match(immune.zika_aba$Gene, tx2gene_order$Gene),], logcpm[match(immune.zika_aba$Gene, rownames(logcpm)),])
head(immune.zika_aba)
immune.zika_aba$Name <- immune_list$Name[match(immune.zika_aba$Gene, immune_list$Gene)]
immune.zika_aba$pathway <- immune_list$Pathway[match(immune.zika_aba$Gene, immune_list$Gene)]
immune.zika_aba$description <- immune_list$Description[match(immune.zika_aba$Gene, immune_list$Gene)]
head(immune.zika_aba)
# Plot
immune_scatter.aba <- ggplot(immune.zika_aba, aes(x=Name, y=log2FoldChange, color = ifelse(padj < 0.05, "Sig", "NonSig"))) +
                  geom_point(size = 5, position=position_jitter(width=.1, height=0)) +
                  scale_color_manual(name="padj", values = c("black", "red")) + 
                  facet_grid(~pathway, scales = "free_x", space = "free_x") + 
                  ylim(-4, 4) +
                  geom_hline(yintercept = 0, linetype = "dashed", color = "black") +
                  ggtitle("AbxA") +
                  ylab("Log2-Fold Change") +
                  xlab("Gene name") +
                  theme(panel.background = element_rect(fill = "white", color = "black"),
                  text = element_text(size=30), axis.title.x = element_blank(), axis.text.x = element_text(angle = 90, hjust = 1))
immune_scatter.aba
ggsave("Immune lfc - AbA.pdf", plot = immune_scatter.aba, width = 30, height = 15, dpi = 600)



# Do the same, but for res05 from res.zika_abp
res05 <- subset(res.zika_abp, padj < 0.05)
sum(res05$padj < 0.05)
res05 <- cbind(rownames(res05), res05)
rownames(res05) <- NULL
colnames(res05) <- c("Gene", "baseMean", "log2FoldChange","lfcSE", "stat", "pvalue", "padj")
res05 <- data.frame(res05, tx2gene_order[match(res05$Gene, tx2gene_order$Gene),], logcpm[match(res05$Gene, rownames(logcpm)),])
head(res05)
res05 <- res05[!grepl("NA", rownames(res05)),]
res05 <- res05[order(res05$log2FoldChange),]
head(res05)
# Filter out immune genes, if present
immune_gene_list <- immune_list$Gene
immune_genes <- subset(res05, Gene %in% immune_gene_list)
immune_genes
# Paste my custom annotations
immune_genes$name <- immune_list$Name[match(immune_genes$Gene, immune_list$Gene)]
immune_genes$pathway <- immune_list$Pathway[match(immune_genes$Gene, immune_list$Gene)]
immune_genes$description <- immune_list$Description[match(immune_genes$Gene, immune_list$Gene)]

ggplot(immune_genes, aes(x=name, y=log2FoldChange, size = padj)) +
  geom_point(position=position_jitter(width=.1, height=0)) +
  scale_size(trans = 'reverse') +
  facet_grid(~pathway, scales = "free_x", space = "free_x") + 
  theme(panel.background = element_rect(fill = "white", color = "black"),
        text = element_text(size=20), axis.text.x = element_text(angle = 50, hjust = 1))

# Filter ALL immune genes, not just res05
immune.zika_abp <- subset(res.zika_abp, rownames(res.zika_abp) %in% immune_list$Gene)
immune.zika_abp <- cbind(rownames(immune.zika_abp), immune.zika_abp)
rownames(immune.zika_abp) <- NULL
colnames(immune.zika_abp) <- c("Gene", "baseMean", "log2FoldChange","lfcSE", "stat", "pvalue", "padj")
immune.zika_abp <- data.frame(immune.zika_abp, tx2gene_order[match(immune.zika_abp$Gene, tx2gene_order$Gene),], logcpm[match(immune.zika_abp$Gene, rownames(logcpm)),])
head(immune.zika_abp)
immune.zika_abp$Name <- immune_list$Name[match(immune.zika_abp$Gene, immune_list$Gene)]
immune.zika_abp$pathway <- immune_list$Pathway[match(immune.zika_abp$Gene, immune_list$Gene)]
immune.zika_abp$description <- immune_list$Description[match(immune.zika_abp$Gene, immune_list$Gene)]
head(immune.zika_abp)
# Plot
immune_scatter.abp <- ggplot(immune.zika_abp, aes(x=Name, y=log2FoldChange, color = ifelse(padj < 0.05, "Sig", "NonSig"))) +
  geom_point(size = 5, position=position_jitter(width=.1, height=0)) +
  scale_color_manual(name="padj", values = c("black", "red")) + 
  facet_grid(~pathway, scales = "free_x", space = "free_x") + 
  ylim(-4, 4) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "black") +
  ggtitle("AbxP") +
  ylab("Log2-Fold Change") +
  xlab("Gene name") +
  theme(panel.background = element_rect(fill = "white", color = "black"),
        text = element_text(size=30), axis.title.x = element_blank(), axis.text.x = element_text(angle = 90, hjust = 1))
immune_scatter.abp
ggsave("Immune lfc - AbP.pdf", plot = immune_scatter.abp, width = 30, height = 15, dpi = 600)



# Merge lfc AbA and AbP into a panel
immune_scatter.merge <- ggarrange(immune_scatter.aba, immune_scatter.abp, align = "v", common.legend = TRUE, legend = "right",
                                  ncol=1, nrow=2)
ggsave("SFig 2 - Immune lfc AbA AbP panels.pdf", plot = immune_scatter.merge, width = 25, height = 25, dpi = 600)



# Do the same, but for res05 from res.zika_aba_abp
head(res.zika_aba_abp)
#   Filtered
res05 <- subset(res.zika_aba_abp, padj < 0.05)
sum(res05$padj < 0.05)
immune_genes <- subset(res05, Gene %in% immune_gene_list)
head(immune_genes)
immune_genes$pathway <- immune_list$Pathway[match(immune_genes$Gene, immune_list$Gene)]
ggplot(immune_genes, aes(x=Gene, y=log2FoldChange, color=pathway, size = padj)) +
  geom_point(position=position_jitter(width=.1,height=0)) +
  scale_size(trans = 'reverse')
#   Unfiltered
immune_genes <- subset(res.zika_aba_abp, Gene %in% immune_gene_list)
head(immune_genes)
immune_genes$Name <- immune_list$Name[match(immune_genes$Gene, immune_list$Gene)]
immune_genes$pathway <- immune_list$Pathway[match(immune_genes$Gene, immune_list$Gene)]
#   Remove NAs from padj
immune_genes <- subset(immune_genes, !is.na(padj))
head(immune_genes)

immune_scatter<- ggplot(immune_genes, aes(x=Name, y=log2FoldChange, color = ifelse(padj < 0.05, "Sig", "NonSig"))) +
                        geom_point(size = 8) +
                        scale_color_manual(name="padj", values = c("black", "red")) + 
                        geom_hline(yintercept = 0, linetype = "dashed", color = "black") +
                        facet_grid(~pathway, scales = "free_x", space = "free_x") + 
                        ylab("Log2-Fold Change") +
                        xlab("Gene name") +
                        theme(panel.background = element_rect(fill = "white", color = "black"),
                              text = element_text(size = 30),
                              axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1))
immune_scatter
ggsave("Fig 2E - Immune lfc - AbA vs AbP.pdf", plot = immune_scatter, width = 30, height = 15, dpi = 600)




### Make Venn diagram of AbP, AbA, and Ctrl DEG (circle is  Zikv/blood)
# We already have res.zika_abp and res.zika_aba
#   These are the DEGs between ab treatments + zikv status
#   Now we just need ctrl - zikv vs blood
# Contrast ZIKV vs blood, for Ctrl
#   1260 genes! 1176 after filtering
#   Redo dds
dds <- DESeqDataSetFromTximport(txi, 
                                colData = metadata,
                                design = ~Bloodmeal_content + Antibiotic_treatment + Bloodmeal_content:Antibiotic_treatment)
colData(dds)
is_expressed <- assay(dds) >= 10
head(is_expressed)
sum(assay(dds)[,1])
colSums(assay(dds))
keep <- rowSums(assay(dds) >= 10) >= 3
table(keep)
dds <- dds[keep,]
dds$Bloodmeal_content = factor(dds$Bloodmeal_content, levels = c("before bloodfeed", "blood only", "blood + ZIKV"))
levels(dds$Bloodmeal_content)
dds$Antibiotic_treatment = factor(dds$Antibiotic_treatment, levels = c("AbxP", "AbxA", "Ctrl"))
levels(dds$Antibiotic_treatment)

dds$Antibiotic_treatment <- relevel(dds$Antibiotic_treatment, ref = "AbxA")
dds$Bloodmeal_content <- relevel(dds$Bloodmeal_content, ref = "blood only")
dds <- DESeq(dds)
resultsNames(dds)

### VENN DIAGRAM - BEFORE VS AFTER BLOODFEED, ZIKV VS BLOOD
venn_ctrl.bef.aft = results(dds, name = "Bloodmeal_content_before.bloodfeed_vs_blood.only", lfcThreshold = 0.05, tidy=TRUE)
res05 <- subset(venn_ctrl.bef.aft, padj < 0.05)
rownames(res05) <- NULL
colnames(res05) <- c("Gene", "baseMean", "log2FoldChange","lfcSE", "stat", "pvalue", "padj")
res05 <- data.frame(res05, tx2gene_order[match(res05$Gene, tx2gene_order$Gene),], logcpm[match(res05$Gene, rownames(logcpm)),])
res05 <- res05[!grepl("NA", rownames(res05)),]
res05 <- res05[order(res05$log2FoldChange),]
sum(res05$padj < 0.05)
venn_ctrl.bef.aft_res <- res05
venn_ctrl.zikv = results(dds, name = "Bloodmeal_content_blood...ZIKV_vs_blood.only",
                        lfcThreshold = 0.05, tidy=TRUE)
res05 <- subset(venn_ctrl.zikv, padj < 0.05)
rownames(res05) <- NULL
colnames(res05) <- c("Gene", "baseMean", "log2FoldChange","lfcSE", "stat", "pvalue", "padj")
res05 <- data.frame(res05, tx2gene_order[match(res05$Gene, tx2gene_order$Gene),], logcpm[match(res05$Gene, rownames(logcpm)),])
res05 <- res05[!grepl("NA", rownames(res05)),]
res05 <- res05[order(res05$log2FoldChange),]
sum(res05$padj < 0.05)
venn_ctrl.zikv_res <- res05
venn_aba.bef.aft = results(dds, name = "Bloodmeal_content_before.bloodfeed_vs_blood.only"  ,
                           lfcThreshold = 0.05, tidy=TRUE)
res05 <- subset(venn_aba.bef.aft, padj < 0.05)
rownames(res05) <- NULL
colnames(res05) <- c("Gene", "baseMean", "log2FoldChange","lfcSE", "stat", "pvalue", "padj")
res05 <- data.frame(res05, tx2gene_order[match(res05$Gene, tx2gene_order$Gene),], logcpm[match(res05$Gene, rownames(logcpm)),])
res05 <- res05[!grepl("NA", rownames(res05)),]
res05 <- res05[order(res05$log2FoldChange),]
sum(res05$padj < 0.05)
venn_aba.bef.aft_res <- res05
venn_aba.zikv = results(dds, name = "Bloodmeal_content_blood...ZIKV_vs_blood.only",
                            lfcThreshold = 0.05, tidy=TRUE)
res05 <- subset(venn_aba.zikv, padj < 0.05)
rownames(res05) <- NULL
colnames(res05) <- c("Gene", "baseMean", "log2FoldChange","lfcSE", "stat", "pvalue", "padj")
res05 <- data.frame(res05, tx2gene_order[match(res05$Gene, tx2gene_order$Gene),], logcpm[match(res05$Gene, rownames(logcpm)),])
res05 <- res05[!grepl("NA", rownames(res05)),]
res05 <- res05[order(res05$log2FoldChange),]
sum(res05$padj < 0.05)
venn_aba.zikv_res <- res05
venn_abp.bef.aft = results(dds, name = "Bloodmeal_content_before.bloodfeed_vs_blood.only",
                        lfcThreshold = 0.05, tidy=TRUE)
res05 <- subset(venn_abp.bef.aft, padj < 0.05)
rownames(res05) <- NULL
colnames(res05) <- c("Gene", "baseMean", "log2FoldChange","lfcSE", "stat", "pvalue", "padj")
res05 <- data.frame(res05, tx2gene_order[match(res05$Gene, tx2gene_order$Gene),], logcpm[match(res05$Gene, rownames(logcpm)),])
res05 <- res05[!grepl("NA", rownames(res05)),]
res05 <- res05[order(res05$log2FoldChange),]
sum(res05$padj < 0.05)
venn_abp.zikv_res <- res05
venn_abp.zikv = results(dds, name = "Bloodmeal_content_blood...ZIKV_vs_blood.only",
                        lfcThreshold = 0.05, tidy=TRUE)
res05 <- subset(venn_abp.zikv, padj < 0.05)
rownames(res05) <- NULL
colnames(res05) <- c("Gene", "baseMean", "log2FoldChange","lfcSE", "stat", "pvalue", "padj")
res05 <- data.frame(res05, tx2gene_order[match(res05$Gene, tx2gene_order$Gene),], logcpm[match(res05$Gene, rownames(logcpm)),])
res05 <- res05[!grepl("NA", rownames(res05)),]
res05 <- res05[order(res05$log2FoldChange),]
sum(res05$padj < 0.05)
venn_abp.zikv_res <- res05


# Plot venn diagram for all antibiotic treatments, BEFORE VS AFTER BLOODFEED
venn_data <- data.frame(AbxA = venn_aba.bef.aft$padj < 0.05,
                        AbxP = venn_abp.bef.aft$padj < 0.05,
                        Ctrl = venn_ctrl.bef.aft$padj < 0.05)
venn <- vennDiagram(venn_data, circle.col = c("orange","blue", "green"), cex=2.5)
dev.copy(pdf,'Venn, for bloodfeeding DEGs.pdf', width=10, height=8)
dev.off()

venn_data <- data.frame(AbxA = venn_aba.zikv$padj < 0.05,
                        AbxP = venn_abp.zikv$padj < 0.05,
                        Ctrl = venn_ctrl.zikv$padj < 0.05)
venn <- vennDiagram(venn_data, circle.col = c("orange","blue", "green"), cex=2.5)
dev.copy(pdf,'Venn, for ZIKV DEGs.pdf', width=10, height=8)
dev.off()





# Plot venn diagram for all antibiotic treatments
venn_data <- data.frame(AbA_vs_Ctrl = res.aba.ctrl$padj < 0.05,
                        AbP_vs_Ctrl = res.abp.ctrl$padj < 0.05)
venn <- vennDiagram(venn_data, circle.col = c("#994F00","#006CD1"), cex=2.5)
dev.copy(pdf,'Venn abx treatment.pdf', width=10, height=8)
dev.off()






######################################################
# KEGG Pathways

if (!requireNamespace("BiocManager", quietly=TRUE))
  install.packages("BiocManager")
BiocManager::install(c("pathview", "gage", "gageData", "GenomicAlignments",
                       "TxDb.Hsapiens.UCSC.hg19.knownGene"))

# Install db package for "fly"
if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
BiocManager::install("org.Dm.eg.db")

library(topGO)
library(KEGGREST)

# Use annotation file, tx2gene_order
# take out genes without GO terms
tx2gene_order <- tx2gene_order[which(tx2gene_order$GO_accession != ""),]
head(tx2gene_order)
# Make reference, assigning GO accession numbers to a gene
gene2GO <- tapply(tx2gene_order$GO_accession, tx2gene_order$Gene, function(x)x)
head(gene2GO)

# Read DE file
DE_b.nb <- read.csv("Sig - Zika response, AbA.csv", stringsAsFactors = F)
head(DE_b.nb)
# Define gene list as 1's if adjP < cutoff, 0, otherwise
geneList_b.nb <- ifelse(DE_b.nb$padj < 0.05, 1, 0)
names(geneList_b.nb) <- unlist(lapply(strsplit(DE_b.nb$Gene, split = ".", fixed = T), function(x)x[1]))
head(geneList_b.nb)
# Creat GO object
GOdata_b.nb <- new("topGOdata",
              ontology = "BP",
              allGenes = geneList_b.nb,
              geneSelectionFun = function(x)(x == 1),
              annot = annFUN.gene2GO, gene2GO = gene2GO)
# Run Fisher's exact test
resultFisher <- runTest(GOdata_b.nb, algorithm = "elim", statistic = "fisher")
# Get results
# Ok if we get all "expected by chance" since our file has only the significant DE genes, not like in example
tab_b.nb <- GenTable(GOdata_b.nb, raw.p.value = resultFisher, topNodes = length(resultFisher@score),
                numChar = 120)
head(tab_b.nb)
write_csv(tab_b.nb, file = "Sig - Bloodmeal vs no blood - GO.csv")

# Read DE file
DE_b.z <- read.csv("Sig - ZIKV vs blood only.csv", stringsAsFactors = F)
head(DE_b.z)
# Define gene list as 1's if adjP < cutoff, 0, otherwise
geneList_b.z <- ifelse(DE_b.z$padj < as0.05, 1, 0)
names(geneList_b.z) <- unlist(lapply(strsplit(DE_b.z$Gene, split = ".", fixed = T), function(x)x[1]))
head(geneList_b.z)
# Creat GO object
GOdata_b.z <- new("topGOdata",
                   ontology = "BP",
                   allGenes = geneList_b.z,
                   geneSelectionFun = function(x)(x == 1),
                   annot = annFUN.gene2GO, gene2GO = gene2GO)
# Run Fisher's exact test
resultFisher <- runTest(GOdata_b.z, algorithm = "elim", statistic = "fisher")
# Get results
# Ok if we get all "expected by chance" since our file has only the significant DE genes, not like in example
tab_b.z <- GenTable(GOdata_b.z, raw.p.value = resultFisher, topNodes = length(resultFisher@score),
                     numChar = 120)
head(tab_b.z)
write_csv(tab_b.z, file = "Sig - ZIKV vs blood only - GO.csv")



# Extract Ae. aegypti info from KEGG website
# https://www.genome.jp/entry/T01053

kegg_abp <- read.csv("Sig - Zika response, AbP.csv", header = T)
kegg_aba <- read.csv("Sig - Zika response, AbA.csv", header = T)
kegg_ctrl <- read.csv("Sig - Zika response, Ctrl.csv", header = T)

# Count number of unieq KEGG IDs
unique(kegg_abp$KEGG_ID)
length(unique(kegg_abp$KEGG_ID))
kegg_abp_pie <- 
  kegg_abp %>%
  group_by(KEGG_ID) %>%
  summarise(kegg_number = n())
kegg_abp_pie

unique(kegg_aba$KEGG_ID)
length(unique(kegg_aba$KEGG_ID))
kegg_aba_pie <- 
  kegg_aba %>%
  group_by(KEGG_ID) %>%
  summarise(kegg_number = n())
kegg_aba_pie

unique(kegg_ctrl$KEGG_ID)
length(unique(kegg_ctrl$KEGG_ID))
kegg_ctrl_pie <- 
  kegg_ctrl %>%
  group_by(KEGG_ID) %>%
  summarise(kegg_number = n())
kegg_ctrl_pie

pie_labels <- c("511" = "Glycan Degradation", "1100" = "Metabolic Processes", "4714" = "Thermogenesis", "NA" = "NA")
pie_ctrl <- pie(kegg_ctrl_pie$kegg_number, labels = pie_labels)

unique(kegg_abp$KEGG_ID)
length(unique(kegg_abp$KEGG_ID))
kegg_abp_pie <- 
  kegg_abp %>%
  group_by(KEGG_ID) %>%
  summarise(kegg_number = n())
kegg_abp_pie

unique(kegg_aba$KEGG_ID)
length(unique(kegg_aba$KEGG_ID))
kegg_aba_pie <- 
  kegg_aba %>%
  group_by(KEGG_ID) %>%
  summarise(kegg_number = n())
kegg_aba_pie


