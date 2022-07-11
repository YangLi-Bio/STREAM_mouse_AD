###############################################################
#                                                             #
# Perform eGRN analyses associated with time points in AD     #
#                                                             #
###############################################################


# Libraries
library(Seurat)
library(dplyr)
library(pbmcapply)
library(ggvenn)


# Global parameters
R.dir <- "/fs/ess/PCON0022/liyang/STREAM/Case_2_AD/Rfile/"
table.dir <- "/fs/ess/PCON0022/liyang/STREAM/Case_2_AD/Tables/"
tool.dir <- "/fs/ess/PCON0022/liyang/r_utilities/functions/"
image.dir <- "/fs/ess/PCON0022/liyang/STREAM/Case_2_AD/Images/"
orig.dir <- "/fs/ess/PCON0022/liyang/STREAM/Case_2_AD/10X_data/"


# Source files
source(paste0(tool.dir, "visual_tools.R"))


# Load Seurat object
obj <- qs::qread(paste0(R.dir, "Obj_clustered.qsave"))
obj@assays
head(colnames(obj))

# Load and meta data
meta.df <- read.csv(paste0(orig.dir, "Multiome_RNA_ATAC_Mouse_Brain_Alzheimers_AppNote_aggr.csv"))
dim(meta.df)
head(meta.df)
sample.ids <- unique(strsplit(colnames(obj), split = "-") %>% 
                       sapply(., function(x) {return(x[2])}))
sample.ids
which(as.numeric(sample.ids) > nrow(meta.df))


# Extract the time points
meta.df$Time.point
time.vec <- sapply(colnames(obj), function(x) {
  meta.df$Time.point[as.numeric(strsplit(x, split = "-")[[1]][2])]
})
obj <- AddMetaData(object = obj, metadata = time.vec, col.name = "Time")
obj$Time
unique(obj$Time)
obj@reductions$wnn.umap
qs::qsave(obj, paste0(R.dir, "Obj_time_points.qsave"))


# Generate the UMAP for time points
obj$Time <- factor(x = obj$Time, levels = c("2.5 months", "5.7 months", 
                                            "13+ months"))
Idents(obj) <- obj$Time
p.time.umap <- get_UMAP(object = obj, txt = "Time",
                     legend.position = "right", 
                     margin.quad = c(0.1, -5.0, 0.1, 0.0))
p.time.umap
qs::qsave(p.time.umap, paste0(R.dir, "UMAP_time_points.qsave"))


# Divide eGRNs into the three time points
eGRN.list <- qs::qread(paste0(R.dir, "eGRNs_list.qsave"))
length(eGRN.list)
core.eGRNs <- lapply(eGRN.list, function(x) {
  list(TF = x$TF, genes = x$genes[x$gene.status], 
       peaks = x$peaks, cells = x$cells, atac.ratio = x$atac.ratio, 
       score = x$score)
})
qs::qsave(core.eGRNs, paste0(R.dir, "Core_eGRNs.qsave"))


# basic information of core eGRNs
core.eGRN.dt <- rbindlist(lapply(core.eGRNs, function(x) {
  list(Genes = length(x$genes), CREs = length(x$peaks), Cells = length(x$cells))
}))
core.eGRN.dt
apply(core.eGRN.dt, 2, range)
qs::qsave(core.eGRN.dt, paste0(R.dir, "Core_eGRN_info.qsave"))


# Generate boxplots for the core eGRNs
p.boxplot.eGRN.info <- get_boxplot(obj = core.eGRN.dt$Genes, jitter.flag = F, 
         y.lab = "Genes", color = "black", fill = "#BE2A3E", 
         width = 200, height = 800) | 
  get_boxplot(obj = core.eGRN.dt$CREs, jitter.flag = F, 
              y.lab = "CREs", color = "black", fill = "#BE2A3E", 
              width = 200, height = 800) | 
  get_boxplot(obj = core.eGRN.dt$Cells, jitter.flag = F, 
              y.lab = "Cells", color = "black", fill = "#BE2A3E", 
              width = 200, height = 800)
qs::qsave(p.boxplot.eGRN.info, paste0(R.dir, "Boxplot_eGRN_genes_CREs_cells.qsave"))


# Calculate ratios of eGRN cells against the three time points
time.cells <- split(colnames(obj), obj$Time)
length(time.cells)
sapply(time.cells, length)
time.cells
eGRN.time.dt <- Reduce("rbind", lapply(core.eGRNs, function(x) {
  sapply(time.cells, function(y) {
    length(intersect(x$cells, y))
  })
}))
dim(eGRN.time.dt)
rownames(eGRN.time.dt) <- seq_along(core.eGRNs)
eGRN.time.dt
qs::qsave(eGRN.time.dt, paste0(R.dir, "eGRN_cells_in_time.qsave"))


# Generate stacked barplots of eGRN cells across time points
table(obj$Time)
eGRN.time.summ <- as(eGRN.time.dt, "sparseMatrix") %>% summary
head(eGRN.time.summ)
stacked.bar.df <- Reduce("rbind", pbmclapply(1:nrow(eGRN.time.summ), function(i) {
  c(rownames(eGRN.time.dt)[eGRN.time.summ[i, 1]], 
    colnames(eGRN.time.dt)[eGRN.time.summ[i, 2]], 
    eGRN.time.summ[i, 3])
}, mc.cores = detectCores())) %>% as.data.frame
rownames(stacked.bar.df) <- 1:nrow(stacked.bar.df)
colnames(stacked.bar.df) <- c("eGRN", "Time", "Ratio")
stacked.bar.df$eGRN <- factor(stacked.bar.df$eGRN, levels = seq_along(core.eGRNs))
stacked.bar.df$Time <- factor(stacked.bar.df$Time, levels = c("2.5 months", "5.7 months", 
                                                         "13+ months"))
stacked.bar.df$Ratio <- as.numeric(stacked.bar.df$Ratio)
qs::qsave(stacked.bar.df, paste0(R.dir, "Stacked_barplot_eGRN_time.qsave"))
p.stacked <- ggplot(stacked.bar.df, aes(fill = Time, y = Ratio, x = eGRN)) + 
  geom_bar(position = "fill", stat = "identity") + 
  scale_fill_manual(values = c("red", "green", "blue"))
p.stacked
qs::qsave(p.stacked, paste0(R.dir, "Stacked_barplots_eGRN_time.qsave"))


# Rule out the 5.7 months
# 2.5 months, 5.7 months, and 13+ months occupy 35.2%, 11.0%, and 53.8% barcodes
# Therefore, we only ocus on the 1.5 and 13+ month conditions
table(obj$Time) / sum(table(obj$Time))
stacked.bars.no.5.7 <- stacked.bar.df
stacked.bars.no.5.7 <- stacked.bars.no.5.7[stacked.bars.no.5.7$Time != "5.7 months",]
dim(stacked.bar.df)
dim(stacked.bars.no.5.7)
p.stacked.no.5.7 <- ggplot(stacked.bars.no.5.7, aes(fill = Time, y = Ratio, x = eGRN)) + 
  geom_bar(position = "fill", stat = "identity") + 
  scale_fill_manual(values = c("red", "blue"))
p.stacked.no.5.7
qs::qsave(p.stacked.no.5.7, paste0(R.dir, "Stacked_barplots_eGRN_time_no_5.7_months.qsave"))


# Partition the eGRNs into three categories: 
# 1. Active in 50%+ cells of 2.5 months
# 2. Active in 50%+ cells of 13+ months
# 3. Other cells
eGRNs.splited <- split(stacked.bars.no.5.7, f = stacked.bars.no.5.7$eGRN)
eGRNs.times <- sapply(eGRNs.splited, function(x) {
  x$Ratio / sum(x$Ratio)
}) %>% t
colnames(eGRNs.times) <- c("2.5 months", "13+ months")
head(eGRNs.times)
dim(eGRNs.times)
eGRNs.2.5 <- core.eGRNs[eGRNs.times[, 1] > 0.50]
length(eGRNs.2.5)
eGRNs.13 <- core.eGRNs[eGRNs.times[, 2] > 0.50]
length(eGRNs.13)
TFs.2.5 <- sapply(eGRNs.2.5, "[[", "TF") %>% table %>% names
TFs.13 <- sapply(eGRNs.13, "[[", "TF") %>% table %>% names
qs::qsave(eGRNs.2.5, paste0(R.dir, "eGRNs_2.5_months.qsave"))
qs::qsave(eGRNs.13, paste0(R.dir, "eGRNs_13_months.qsave"))


# Get the intersected cells where both eGRNs.2.5 and eGRNs.13 are active

p.venn.TFs <- ggvenn(list("2.5 months" = TFs.2.5, "13+ months" = TFs.13), 
       columns = c("2.5 months", "13+ months"), 
       # show_elements = T,
       fill_color = c("#FF0000", "#0000FF"), 
       text_size = 6, 
       stroke_size = 0.3)
qs::qsave(p.venn.TFs, paste0(R.dir, "Venn_diagrams_TFs_2.5_13.qsave"))


# Get the heatmap of TF expression, chromatin accessibility, and target expression across
# 2.5 months and 13+ months

# We cannot observe obvious difference between expression of target genes between 
# 2.5 months v.s. 13+ months. Hence, we will not include this heatmap.

# We cannot observe obvious difference in expression in eGRN cells between 2.5 months v.s. 
# 13+ months


# Count the number of target genes in eGRNs
genes.2.5 <- Reduce("union", sapply(eGRNs.2.5, "[[", "genes"))
length(genes.2.5)
TFs.2.5 <- Reduce("union", sapply(eGRNs.2.5, "[[", "TF")) %>% tolower %>% capitalize
length(TFs.2.5)
length(intersect(genes.2.5, TFs.2.5))
genes.13 <- Reduce("union", sapply(eGRNs.13, "[[", "genes"))
length(genes.13)
TFs.13 <- Reduce("union", sapply(eGRNs.13, "[[", "TF")) %>% tolower %>% capitalize
length(TFs.13)
length(intersect(genes.13, TFs.13))


# Create assay for the RNA expression
eGRNs.two.times <- c(eGRNs.2.5, eGRNs.13)
length(eGRNs.two.times)
rna.m <- GetAssayData(obj, slot = "data", assay = "RNA")
qs::qsave(eGRNs.two.times, paste0(R.dir, "eGRNs_two_times.qsave"))
eGRN.exp <- Reduce("rbind", pbmclapply(eGRNs.two.times, function(x) {
  v <- apply(rna.m[x$genes, , drop = F], 2, mean)
  v[setdiff(colnames(rna.m), x$cells)] <- 0
  v
}), mc.cores = detectCores())
rownames(eGRN.exp) <- seq_along(eGRNs.two.times)
dim(eGRN.exp)
head(eGRN.exp[1:5, 1:5])
sapply(eGRNs.two.times, "[[", "TF")
obj[["eGRN.exp"]] <- CreateAssayObject(counts = eGRN.exp)
obj[['eGRN.exp']][1:5, 1:5]


# Create assay for the accessibility
atac.m <- GetAssayData(obj, slot = "data", assay = "ATAC")
eGRN.acc <- Reduce("rbind", pbmclapply(eGRNs.two.times, function(x) {
  v <- apply(atac.m[x$peaks, , drop = F], 2, mean)
  v[setdiff(colnames(atac.m), x$cells)] <- 0
  v
}, mc.cores = detectCores()))
dim(eGRN.acc)
rownames(eGRN.acc) <- seq_along(eGRNs.two.times)
eGRN.acc[1:5, 1:5]
obj[["eGRN.acc"]] <- CreateAssayObject(counts = eGRN.acc)
obj[["eGRN.acc"]][1:5, 1:5]
qs::qsave(obj, paste0(R.dir, "Obj_eGRN_assays.qsave"))


# There are too few cells in FeaturePlots.
# Hence, I will not use FeaturePlots.

# Get DEGs and DARs belonging to eGRNs
eGRN.genes <- Reduce("union", lapply(eGRNs.two.times, "[[", "genes"))
length(eGRN.genes)
eGRN.peaks <- Reduce("union", lapply(eGRNs.two.times, "[[", "peaks"))
length(eGRN.peaks)
Idents(obj) <- obj$Time
DefaultAssay(obj) <- "RNA"
DEGs <- FindMarkers(object = obj, features = eGRN.genes, 
                    ident.1 = "2.5 months", ident.2 = "13+ months")
qs::qsave(DEGs, paste0(R.dir, "DEGs.qsave"))
DefaultAssay(obj) <- "ATAC"
# DARs <- FindMarkers(object = obj, features = eGRN.peaks, 
#                     ident.1 = "2.5 months", ident.2 = "13+ months")
# qs::qsave(DARs, paste0(R.dir, "DARs.qsave"))


# Filter DEGs
DEGs <- DEGs[DEGs$p_val_adj < 0.05,]
dim(DEGs)
# Only 181 DEGs are retained for 2.5 months and 13+ months.
# Therefore, I will discard this plot.


# AD-associated proteins: JUND, MAP1B, FOS, MFGE8, JUNB, and JUN
# Ref: Xu J, Zhang P, Huang Y, et al. Multimodal single-cell/nucleus RNA sequencing data analysis uncovers molecular networks between disease-associated microglia and astrocytes with implications for drug repurposing in Alzheimer's disease[J]. Genome research, 2021, 31(10): 1900-1912.
TFs <- sapply(eGRNs.two.times, "[[", "TF") %>% tolower %>% capitalize # we have Jund, Jun, and Fos
TFs
table(TFs)

# Nfe2: Osama A, Zhang J, Yao J, Yao X, Fang J. Nrf2: a dark horse in Alzheimer's disease treatment. Ageing Res Rev. 2020 Dec;64:101206. doi: 10.1016/j.arr.2020.101206. Epub 2020 Nov 2. PMID: 33144124.
# Plag1: Liu C, Chyr J, Zhao W, Xu Y, Ji Z, Tan H, Soto C, Zhou X; Alzheimer’s Disease Neuroimaging Initiative. Genome-Wide Association and Mechanistic Studies Indicate That Immune Response Contributes to Alzheimer's Disease Development. Front Genet. 2018 Sep 24;9:410. doi: 10.3389/fgene.2018.00410. PMID: 30319691; PMCID: PMC6166008.

# Esr2: Ulhaq ZS, Garcia CP. Estrogen receptor beta (ESR2) gene polymorphism and susceptibility to dementia. Acta Neurol Belg. 2021 Oct;121(5):1281-1293. doi: 10.1007/s13760-020-01360-z. Epub 2020 Apr 25. PMID: 32335869.

# Nr2c2: Dharshini S A P, Taguchi Y H, Gromiha M M. Exploring the selective vulnerability in Alzheimer disease using tissue specific variant analysis[J]. Genomics, 2019, 111(4): 936-949.

# Elk4: Liu Q, Zhu L, Liu X, et al. TRA2A-induced upregulation of LINC00662 regulates blood-brain barrier permeability by affecting ELK4 mRNA stability in Alzheimer’s microenvironment[J]. RNA biology, 2020, 17(9): 1293-1308.

# Nfatc2: Manocha GD, Ghatak A, Puig KL, Kraner SD, Norris CM, Combs CK. NFATc2 Modulates Microglial Activation in the AβPP/PS1 Mouse Model of Alzheimer's Disease. J Alzheimers Dis. 2017;58(3):775-787. doi: 10.3233/JAD-151203. PMID: 28505967; PMCID: PMC6265241.

# NKX3-1: Khayer, N., Jalessi, M., Jahanbakhshi, A. et al. Nkx3-1 and Fech genes might be switch genes involved in pituitary non-functioning adenoma invasiveness. Sci Rep 11, 20943 (2021). https://doi.org/10.1038/s41598-021-00431-2

# AR: Ferrari R, Dawoodi S, Raju M, Thumma A, Hynan LS, Maasumi SH, Reisch JS, O'Bryant S, Jenkins M, Barber R, Momeni P. Androgen receptor gene and sex-specific Alzheimer's disease. Neurobiol Aging. 2013 Aug;34(8):2077.e19-20. doi: 10.1016/j.neurobiolaging.2013.02.017. Epub 2013 Mar 29. PMID: 23545426; PMCID: PMC4012749.


Jund.genes <- Reduce("union", lapply(eGRNs.two.times[c(5, 7)], "[[", "genes"))
length(Jund.genes)
Jund.cells <- Reduce("union", lapply(eGRNs.two.times[c(5, 7)], "[[", "cells"))
length(Jund.cells)


# Rearrange genes and cells to showcase the matrix
Jund.5 <- eGRNs.two.times[[5]]
Jund.7 <- eGRNs.two.times[[7]]

# Genes
genes.1 <- setdiff(Jund.5$genes, Jund.7$genes)
genes.2 <- intersect(Jund.5$genes, Jund.7$genes)
genes.3 <- setdiff(Jund.7$genes, Jund.5$genes)

# Peaks
peaks.1 <- setdiff(Jund.5$peaks, Jund.7$peaks)
peaks.2 <- intersect(Jund.5$peaks, Jund.7$peaks)
peaks.3 <- setdiff(Jund.7$peaks, Jund.5$peaks)

# Cells
cells.1 <- setdiff(Jund.5$cells, Jund.7$cells)
cells.2 <- intersect(Jund.5$cells, Jund.7$cells)
cells.3 <- setdiff(Jund.7$cells, Jund.5$cells)


# Build matrix for RNA
Jund.rna <- rna.m[c(genes.1, genes.2, genes.3), 
                     c(cells.1, cells.2, cells.3)]
Jund.rna <- Jund.rna / apply(Jund.rna, 1, max)
dim(Jund.rna)
class(Jund.rna)
Jund.exp.summ <- summary(Jund.rna)
Jund.exp.tile <- rbindlist(apply(Jund.exp.summ, 1, function(x) {
  list(rownames(Jund.rna)[x[1]], colnames(Jund.rna)[x[2]], 
       x[3])
}))
colnames(Jund.exp.tile) <- c("Gene", "Cell", "Expression")
Jund.exp.tile$Gene <- factor(Jund.exp.tile$Gene, levels = rownames(Jund.rna))
Jund.exp.tile$Cell <- factor(Jund.exp.tile$Cell, levels = colnames(Jund.rna))
p.Jund.exp <- ggplot(Jund.exp.tile, aes(x = Cell, y = Gene, fill = Expression)) + 
  geom_tile(color = "white") + 
  scale_fill_gradient(low = "#FFDEDE", high = "#FF6464") + coord_fixed() + 
  theme(axis.text.x = element_blank(), axis.ticks.x = element_blank())
p.Jund.exp
qs::qsave(p.Jund.exp, paste0(R.dir, "Heatmap_Jund_expression.qsave"))


# Build matrix for ATAC
# Jund.atac <- atac.m[c(peaks.1, peaks.2, peaks.3),
#                   c(cells.1, cells.2, cells.3)] > 0 %>% as.numeric
# # Jund.atac <- Jund.atac / apply(Jund.atac, 1, max)
# dim(Jund.atac)
# class(Jund.atac)
# range(Jund.atac)
# Jund.acc.summ <- summary(Jund.atac)
# Jund.acc.summ$x <- as.numeric(Jund.acc.summ$x)
# Jund.acc.tile <- rbindlist(apply(Jund.acc.summ, 1, function(x) {
#   list(rownames(Jund.atac)[x[1]], colnames(Jund.atac)[x[2]],
#        x[3])
# }))
# colnames(Jund.acc.tile) <- c("CRE", "Cell", "Accessibility")
# Jund.acc.tile$CRE <- factor(Jund.acc.tile$CRE, levels = rownames(Jund.atac))
# Jund.acc.tile$Cell <- factor(Jund.acc.tile$Cell, levels = colnames(Jund.atac))
# p.Jund.acc <- ggplot(Jund.acc.tile, aes(x = Cell, y = CRE, fill = Accessibility)) +
#   geom_tile(color = "white") +
#   scale_fill_gradient(low = "#FFDEDE", high = "#FF6464") + coord_fixed() +
#   theme(axis.text.x = element_blank(), axis.ticks.x = element_blank())
# p.Jund.acc

# We cannot observe obvious difference in accessibility between eGRNs.
# Hence, we will not include this plot.
# Now, there are three panels.


# Comparative analyses
# 1. KEGG pathway analyses and comparison
# 2. eGRNs regulated by the same TFs, e.g., Ar, Fos, Nfe2, Elk4, Esr2, and Zbtb33
# 3. Difference in expression and accessibility, TF expression
# 4. Emphasize the CREs (composition and accessibility) in the eGRNs overrepresented in one of the eGRNs
# 5. Inspect the pathways specific to the key gene or their associated CREs