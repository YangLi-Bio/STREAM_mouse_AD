###############################################################
#                                                             #
#             Annotate cell types using markers               #
#                                                             #
###############################################################


# Libraries
library(Seurat)
library(Signac)
library(readxl)
library(pbmcapply)
library(pbapply)
library(dplyr)


# Global parameters
R.dir <- "/fs/ess/PCON0022/liyang/STREAM/Case_2_AD/Rfile/"
table.dir <- "/fs/ess/PCON0022/liyang/STREAM/Case_2_AD/Tables/"
tool.dir <- "/fs/ess/PCON0022/liyang/r_utilities/functions/"
image.dir <- "/fs/ess/PCON0022/liyang/STREAM/Case_2_AD/Images/"


# Source codes
source(paste0(tool.dir, "visual_tools.R"))


# Load Seurat object
obj <- qs::qread(paste0(R.dir, "Obj_clustered.qsave"))
obj
obj@assays
levels(obj$seurat_clusters)


# Load the marker gene table
marker.df <- read_excel(paste0(table.dir, "Brain_cell_type_markers.xlsx"), 
                        sheet = 1)
class(marker.df)
dim(marker.df)
marker.list <- apply(marker.df, 2, function(x) {
  capitalize(tolower(x[!is.na(x)]))
})
marker.list
range(sapply(marker.list, length))
qs::qsave(marker.list, paste0(R.dir, "Marker_list.qsave"))


# Build marker-cluster expression matrix
rna.m <- GetAssayData(object = obj, assay = "RNA", slot = "data")
setdiff(Reduce("union", marker.list), rownames(rna.m))
marker.m <- rna.m[Reduce("union", marker.list),]
marker.cluster.m <- apply(marker.m, 1, function(r) {
  sapply(levels(obj$seurat_clusters), function(x) {
    mean(r[which(obj$seurat_clusters == x)])
  })
}) %>% t
rownames(marker.cluster.m)
colnames(marker.cluster.m)
marker.cluster.m[1:4, 1:4]
qs::qsave(marker.cluster.m, paste0(R.dir, "Marker_cluster_matrix.qsave"))


# Build marker-cluster activity matrix
act.m <- GeneActivity(
  object = obj,
  assay = "ATAC",
  features = Reduce("union", marker.list),
  extend.upstream = 2000,
  extend.downstream = 0,
  biotypes = "protein_coding",
  max.width = 5e+05,
  process_n = 2000,
  gene.id = FALSE,
  verbose = TRUE
)
dim(act.m)
setdiff(Reduce("union", marker.list), rownames(act.m))
marker.act.m <- act.m[Reduce("union", marker.list),]
act.cluster.m <- apply(marker.act.m, 1, function(r) {
  sapply(levels(obj$seurat_clusters), function(x) {
    mean(r[which(obj$seurat_clusters == x)])
  })
}) %>% t
rownames(act.cluster.m)
colnames(act.cluster.m)
act.cluster.m[1:4, 1:4]
dim(act.cluster.m)
qs::qsave(act.cluster.m, paste0(R.dir, "Activity_cluster_matrix.qsave"))


# generate marker-cluster expression matrix
cell.type.list <- do.call("c", lapply(seq_along(marker.list), function(i) {
  rep(names(marker.list)[i], length(marker.list[[i]]))
}))
length(cell.type.list) == nrow(marker.cluster.m)
row.split <- sapply(cell.type.list, function(x) {
  substr(x, 1, 3)
})
p.marker.exp.heatmap <- get_complex_heatmap(obj = marker.cluster.m, top.fill = NULL, 
                    top.label = NULL, left.fill = NULL, 
                    left.label = names(marker.list), norm = T, 
                    col.split = NULL, row.split = row.split)
qs::qsave(p.marker.exp.heatmap, paste0(R.dir, "Marker_expression_heatmap.qsave"))


# generate marker-cluster activity matrix
p.activity.heatmap <- get_complex_heatmap(obj = act.cluster.m, top.fill = NULL, 
                                            top.label = NULL, left.fill = NULL, 
                                            left.label = names(marker.list), norm = T, 
                                            col.split = NULL, row.split = row.split)
qs::qsave(p.activity.heatmap, paste0(R.dir, "Activity_expression_heatmap.qsave"))


# Annotate cell types
p.exp.act <- p.marker.exp.heatmap | p.activity.heatmap
p.exp.act
qs::qsave(p.exp.act, paste0(R.dir, "Expression_activity_heatmaps.qsave"))
