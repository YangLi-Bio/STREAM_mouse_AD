###############################################################
#                                                             #
# Perform eGRN analyses associated with time points in AD     #
#                                                             #
###############################################################


# Libraries
library(Seurat)
library(dplyr)
library(pbmcapply)


# Global parameters
R.dir <- "/fs/ess/PCON0022/liyang/STREAM/Case_2_AD/Rfile/"
table.dir <- "/fs/ess/PCON0022/liyang/STREAM/Case_2_AD/Tables/"
tool.dir <- "/fs/ess/PCON0022/liyang/r_utilities/functions/"
image.dir <- "/fs/ess/PCON0022/liyang/STREAM/Case_2_AD/Images/"
orig.dir <- "/fs/ess/PCON0022/liyang/STREAM/Case_2_AD/10X_data/"


# Load Seurat object
obj <- qs::qread(paste0(R.dir, "Obj_QC.qsave"))
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
qs::qsave(obj, paste0(R.dir, "Obj_time_points.qsave"))
