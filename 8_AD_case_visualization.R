##############################################################
#                                                            #
#                  Make figures for Case AD                  #
#                                                            #
##############################################################


# Global parameters
tool.dir <- "/fs/ess/PCON0022/liyang/r_utilities/functions/"
R.dir <- "/fs/ess/PCON0022/liyang/STREAM/Case_2_AD/Rfile/"
table.dir <- "/fs/ess/PCON0022/liyang/STREAM/Case_2_AD/Tables/"


# Source codes
source(paste0(tool.dir, "transcriptome_tools.R"))
source(paste0(tool.dir, "cistrome_tools.R"))
source(paste0(tool.dir, "cistrome_tools.R"))


# Libraries
library(Seurat)
library(Signac)
library(dplyr)
library(qs)
library(ggvenn)
library(Hmisc)


# Load eGRNs
eGRNs.2.5 <- qs::qread(paste0(R.dir, "eGRNs_2.5_months.qsave"))
eGRNs.13 <- qs::qread(paste0(R.dir, "eGRNs_13+_months.qsave"))


# Panel A: UMAP plot of three time points
p.umap <- qs::qread(paste0(R.dir, "UMAP_time_points.qsave"))
p.umap


# Panel B: Venn diagram for TFs between 2.5 months v.s. 13+ months
# Path: Panel_B_UMAP_time_points
p.Venn.TFs
TFs.2.5
TFs.13
length(intersect(TFs.2.5, TFs.13))
setdiff(TFs.2.5, TFs.13)
setdiff(TFs.13, TFs.2.5)


# Panel C: eGRNs and pathways of Jund
links.Jund.2.5
length(links.Jund.2.5)
links.Jund.2.5$gene
intersect(links.Jund.2.5$gene, AD.genes)
GRangesToString(links.Jund.2.5[links.Jund.2.5$gene %in% AD.genes])


# Prepare Cytoscape files
prepare_Cytoscape(gr = links.Jund.2.5, org = "mm10", 
                  obj = obj.2.5, cells = Cores.2.5[[11]]$cells,
                  path = paste0(table.dir, "Jund_2.5_eGRN/"), 
                  prefix = "Jund",
                  binary.accessibility = F)


# Barplots of pathways


# UCSC Genome Browser

