###############################################################
#                                                             #
#          Convert 10X Multiome data into Seurat object       #
#                                                             #
###############################################################


# Global variables
tool.dir <- "/fs/ess/PCON0022/liyang/r_utilities/functions/"
code.dir <- "/fs/ess/PCON0022/liyang/STREAM/Case_2_AD/Codes/"
R.dir <- "/fs/ess/PCON0022/liyang/STREAM/Case_2_AD/Rfile/"
image.dir <- "/fs/ess/PCON0022/liyang/STREAM/Case_2_AD/Images/"


# Local variables
orig.dir <- "/fs/ess/PCON0022/liyang/STREAM/Case_2_AD/10X_data/"
h5.file <- "Multiome_RNA_ATAC_Mouse_Brain_Alzheimers_AppNote_filtered_feature_bc_matrix.h5"
frag.file <- "Multiome_RNA_ATAC_Mouse_Brain_Alzheimers_AppNote_atac_fragments.tsv.gz"


# # Preprocess the 10X Genomics Multiome data
# source(paste0(tool.dir, "IO_tools.R"))
# obj <- load_10X_Multiome(h5.file = paste0(orig.dir, h5.file), 
#                          frag.file = paste0(orig.dir, frag.file),  
#                          org = "mm10", image.path = paste0(image.dir, 
#                                                            "Violin_plots_QC.png"))
# 
# 
# # Save the Seurat object and violin plots
# qs::qsave(obj, paste0(R.dir, "Obj_preprocessed.qsave"))


# Generate a small object for parameter tuning of STREAM
obj <- qs::qread(paste0(R.dir, "Obj_preprocessed.qsave"))
source(paste0(tool.dir, "multiome_tools.R"))
obj.sampled <- sample_data(obj = obj, assays = c("RNA", "ATAC"),
                        n.assays = c(1000, 1000), n.cells = 1000, 
                        method = "variable")
qs::qsave(obj.sampled, paste0(R.dir, "Obj_sampled.qsave"))
obj.sampled@assays$RNA
obj.sampled@assays$ATAC


# Quality control (QC)
obj.subsetted <- subset(x = obj, 
                        subset = nCount_ATAC < 100000 &
                          nCount_RNA < 25000 &
                          nCount_ATAC > 1000 &
                          nCount_RNA > 1000 &
                          nucleosome_signal < 1 &
                          TSS.enrichment > 1)
ncol(obj)
ncol(obj.subsetted)
qs::qsave(obj.subsetted, paste0(R.dir, "Obj_QC.qsave"))
