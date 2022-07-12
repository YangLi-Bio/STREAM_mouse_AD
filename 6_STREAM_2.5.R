###############################################################
#                                                             #
#               Test STREAM on a small dataset                #
#                                                             #
###############################################################


# Global parameters
tool.dir <- "/fs/ess/PCON0022/liyang/r_utilities/functions/"
code.dir <- "/fs/ess/PCON0022/liyang/STREAM/Case_2_AD/Codes/"
R.dir <- "/fs/ess/PCON0022/liyang/STREAM/Case_2_AD/Rfile/"
out.dir <- "/fs/ess/scratch/PCON0022/liyang/stream/Case_2_AD/Rfile/Rfile_2.5/"


# Run STREAM
source("/fs/ess/PCON0022/liyang/STREAM/Codes/stream/R/stream.R")
obj <- qs::qread(paste0(R.dir, "Obj_2.5_months.qsave"))
eGRNs.sampled <- run_stream(obj = obj, var.genes = 3000, top.peaks = 3000,
                       min.cells = 10, out.dir = out.dir, org = "mm10",
                       top.ngenes = 15, c.cutoff = 1.0, n.blocks = 500,
                       seed.ratio = 0.30, cicero.covar = 0.00,
                       signac.score = 0.00, min.eGRNs = 100,
                       peak.assay = "ATAC", sim.mode = "both",
                       cover.blocks = 10, KL = 6, expand.dist = 250000, 
                       expand.cutoff = 0.70)
qs::qsave(eGRNs.sampled, paste0(out.dir, "eGRNs_list.qsave"))