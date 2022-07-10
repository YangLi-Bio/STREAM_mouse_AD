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

# Use a volcano plot or MA plot to visualize the difference in expression, accessibility, and TF
# between 2.5 months v.s. 13+ months
rna.m <- GetAssayData(obj, "data", "RNA")
atac.m <- GetAssayData(obj, "data", "ATAC")
cells.2.5 <- colnames(obj)[which(obj$Time == "2.5 months")]
length(cells.2.5)
head(cells.2.5)
cells.13 <- colnames(obj)[which(obj$Time == "13+ months")]
length(cells.13)
head(cells.13)
eGRNs.timed <- c(eGRNs.2.5, eGRNs.13)
cells.eGRN <- Reduce("union", lapply(eGRNs.timed, "[[", "cells"))
