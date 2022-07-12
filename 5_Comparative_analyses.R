##############################################################
#                                                            #
# Comparative analyses between eGRNs inferred in 2.5 months  #
# v.s. 13+ months                                            #
#                                                            #
##############################################################


# Global parameters
tool.dir <- "/fs/ess/PCON0022/liyang/r_utilities/functions/"
R.dir <- "/fs/ess/PCON0022/liyang/STREAM/Case_2_AD/Rfile/"


# Source codes
source(paste0(tool.dir, "transcriptome_tools.R"))
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
length(eGRNs.2.5)
length(eGRNs.13)
TFs.2.5 <- Reduce("union", lapply(eGRNs.2.5, "[[", "TF"))
TFs.2.5
length(TFs.2.5)
TFs.13 <- Reduce("union", lapply(eGRNs.13, "[[", "TF"))
TFs.13
length(TFs.13)
TFs.intersected <- intersect(TFs.2.5, TFs.13)
length(TFs.intersected)
p.Venn.TFs <- ggvenn(list("2.5 months" = TFs.2.5, "13+ months" = TFs.13), 
                     fill_color = c("#ff0000", "#0000FF"), stroke_size = 0.3, 
                     set_name_size = 5, fill_alpha = 1.0)
p.Venn.TFs
qs::qsave(p.Venn.TFs, paste0(R.dir, "Venn_TFs_intersected.qsave"))


# Important TFs

# AD-associated proteins: JUND, MAP1B, FOS, MFGE8, JUNB, and JUN
# Ref: Xu J, Zhang P, Huang Y, et al. Multimodal single-cell/nucleus RNA sequencing data analysis uncovers molecular networks between disease-associated microglia and astrocytes with implications for drug repurposing in Alzheimer's disease[J]. Genome research, 2021, 31(10): 1900-1912.

# Nfe2: Osama A, Zhang J, Yao J, Yao X, Fang J. Nrf2: a dark horse in Alzheimer's disease treatment. Ageing Res Rev. 2020 Dec;64:101206. doi: 10.1016/j.arr.2020.101206. Epub 2020 Nov 2. PMID: 33144124.

# Plag1: Liu C, Chyr J, Zhao W, Xu Y, Ji Z, Tan H, Soto C, Zhou X; Alzheimer’s Disease Neuroimaging Initiative. Genome-Wide Association and Mechanistic Studies Indicate That Immune Response Contributes to Alzheimer's Disease Development. Front Genet. 2018 Sep 24;9:410. doi: 10.3389/fgene.2018.00410. PMID: 30319691; PMCID: PMC6166008.

# Esr2: Ulhaq ZS, Garcia CP. Estrogen receptor beta (ESR2) gene polymorphism and susceptibility to dementia. Acta Neurol Belg. 2021 Oct;121(5):1281-1293. doi: 10.1007/s13760-020-01360-z. Epub 2020 Apr 25. PMID: 32335869.

# Nr2c2: Dharshini S A P, Taguchi Y H, Gromiha M M. Exploring the selective vulnerability in Alzheimer disease using tissue specific variant analysis[J]. Genomics, 2019, 111(4): 936-949.

# Elk4: Liu Q, Zhu L, Liu X, et al. TRA2A-induced upregulation of LINC00662 regulates blood-brain barrier permeability by affecting ELK4 mRNA stability in Alzheimer’s microenvironment[J]. RNA biology, 2020, 17(9): 1293-1308.

# Nfatc2: Manocha GD, Ghatak A, Puig KL, Kraner SD, Norris CM, Combs CK. NFATc2 Modulates Microglial Activation in the AβPP/PS1 Mouse Model of Alzheimer's Disease. J Alzheimers Dis. 2017;58(3):775-787. doi: 10.3233/JAD-151203. PMID: 28505967; PMCID: PMC6265241.

# NKX3-1: Khayer, N., Jalessi, M., Jahanbakhshi, A. et al. Nkx3-1 and Fech genes might be switch genes involved in pituitary non-functioning adenoma invasiveness. Sci Rep 11, 20943 (2021). https://doi.org/10.1038/s41598-021-00431-2

# Ar: Ferrari R, Dawoodi S, Raju M, Thumma A, Hynan LS, Maasumi SH, Reisch JS, O'Bryant S, Jenkins M, Barber R, Momeni P. Androgen receptor gene and sex-specific Alzheimer's disease. Neurobiol Aging. 2013 Aug;34(8):2077.e19-20. doi: 10.1016/j.neurobiolaging.2013.02.017. Epub 2013 Mar 29. PMID: 23545426; PMCID: PMC4012749.


# Get Core eGRNs
Cores.2.5 <- lapply(eGRNs.2.5, function(x) {
  xx <- x
  xx$genes <- x$genes[x$gene.status]
  xx$gene.status <- NULL
  xx$TF <- capitalize(tolower(x$TF))
  if (x$TF == "MAF::NFE2") {
    xx$TF <- "Nfe2"
  }
  if (x$TF == "STAT1::STAT2") {
    xx$TF <- "Stat2"
  }
  if (x$TF == "NFIC::TLX1") {
    xx$TF <- "Tlx1"
  }
  xx
})
Cores.2.5[[1]]

Cores.13 <- lapply(eGRNs.13, function(x) {
  xx <- x
  xx$genes <- x$genes[x$gene.status]
  xx$gene.status <- NULL
  xx$TF <- capitalize(tolower(x$TF))
  if (x$TF == "MAF::NFE2") {
    xx$TF <- "Nfe2"
  }
  if (x$TF == "STAT1::STAT2") {
    xx$TF <- "Stat2"
  }
  if (x$TF == "NFIC::TLX1") {
    xx$TF <- "Tlx1"
  }
  xx
})
Cores.13[[1]]

freq.2.5 <- sapply(Cores.2.5, "[[", "TF")
table(freq.2.5)
TFs.2.5 <- names(table(freq.2.5))

freq.13 <- sapply(Cores.13, "[[", "TF")
table(freq.13)
TFs.13 <- names(table(freq.13))

qs::qsave(Cores.2.5, paste0(R.dir, "Core_eGRNs_2.5.qsave"))
qs::qsave(Cores.13, paste0(R.dir, "Core_eGRNs_13.qsave"))


# Select eGRNs for comparison
TFs.2.5
TFs.13
ref.TFs <- c("Jund", "Fos", "Jun", "Nfe2", "Plag1", "Esr2", "Elk4", 
             "Nkx3-1", "Ar")
intersect(TFs.2.5, ref.TFs)
intersect(TFs.13, ref.TFs)
setdiff(TFs.2.5, TFs.13)
setdiff(TFs.13, TFs.2.5)

# We will explore a pair of eGRNs of Jund, Jun, Fos in 2.5 months v.s. 13+ months
# Retain other Cytoscape networks as Supplementary Figures

# Other TFs (Ar, Elk4, Esr2, and Nfe2) may be put to Supplemtary Figures
# Nfatc2 and Plag1


# Pathway enrichment analyses
gene.list.2.5 <- lapply(Cores.2.5, "[[", "genes")
gene.list.13 <- lapply(Cores.13, "[[", "genes")
pathways.2.5 <- run_GO_and_KEGG(genes.ll = gene.list.2.5, 
                                dbs = "KEGG", org = "mouse")$KEGG_2019_Mouse
pathways.13 <- run_GO_and_KEGG(genes.ll = gene.list.13, 
                               dbs = "KEGG", org = "mouse")$KEGG_2019_Mouse
qs::qsave(pathways.2.5, paste0(R.dir, "Pathways_2.5.qsave"))
qs::qsave(pathways.13, paste0(R.dir, "Pathways_13.qsave"))


# Check pathways of Jund
pathways.2.5 <- pathways.2.5[pathways.2.5$Adjusted.P.value < 0.05,]
pathways.13 <- pathways.13[pathways.13$Adjusted.P.value < 0.05,]

Jund.2.5.ids <- which(freq.2.5 == "Jund")
which(grepl("Alzheimer", pathways.2.5[pathways.2.5$Id %in% Jund.2.5.ids, 2]))
which(grepl("rain", pathways.2.5[pathways.2.5$Id %in% Jund.2.5.ids, 2]))
which(grepl("ipid", pathways.2.5[pathways.2.5$Id %in% Jund.2.5.ids, 2]))
# Make it!

Jund.13.ids <- which(freq.13 == "Jund")
which(grepl("Alzheimer", pathways.13[pathways.13$Id %in% Jund.13.ids, 2]))
which(grepl("rain", pathways.13[pathways.13$Id %in% Jund.13.ids, 2]))
# None!


# Check pathways of Jun
Jun.2.5.ids <- which(freq.2.5 == "Jun")
which(grepl("Alzheimer", pathways.2.5[pathways.2.5$Id %in% Jun.2.5.ids, 2]))
which(grepl("rain", pathways.2.5[pathways.2.5$Id %in% Jun.2.5.ids, 2]))
# None!

Jun.13.ids <- which(freq.13 == "Jun")
which(grepl("Alzheimer", pathways.13[pathways.13$Id %in% Jun.2.5.ids, 2]))
which(grepl("rain", pathways.13[pathways.13$Id %in% Jun.2.5.ids, 2]))
# None!


# Check pathways of Fos
Fos.2.5.ids <- which(freq.2.5 == "Fos")
which(grepl("Alzheimer", pathways.2.5[pathways.2.5$Id %in% Fos.2.5.ids, 2]))
which(grepl("rain", pathways.2.5[pathways.2.5$Id %in% Fos.2.5.ids, 2]))
# Make it!

Fos.13.ids <- which(freq.13 == "Fos")
which(grepl("Alzheimer", pathways.13[pathways.13$Id %in% Fos.13.ids, 2]))
which(grepl("rain", pathways.13[pathways.13$Id %in% Fos.13.ids, 2]))
# None!


# Check pathways of Ar
Ar.2.5.ids <- which(freq.2.5 == "Ar")
which(grepl("Alzheimer", pathways.2.5[pathways.2.5$Id %in% Ar.2.5.ids, 2]))
which(grepl("rain", pathways.2.5[pathways.2.5$Id %in% Ar.2.5.ids, 2]))
# None!

Ar.13.ids <- which(freq.13 == "Ar")
which(grepl("Alzheimer", pathways.13[pathways.13$Id %in% Ar.13.ids, 2]))
which(grepl("rain", pathways.13[pathways.13$Id %in% Ar.13.ids, 2]))
# None!


# Check pathways of Elk4
Elk4.2.5.ids <- which(freq.2.5 == "Elk4")
which(grepl("Alzheimer", pathways.2.5[pathways.2.5$Id %in% Elk4.2.5.ids, 2]))
which(grepl("rain", pathways.2.5[pathways.2.5$Id %in% Elk4.2.5.ids, 2]))
# Bingo!

Elk4.13.ids <- which(freq.13 == "Elk4")
which(grepl("Alzheimer", pathways.13[pathways.13$Id %in% Elk4.13.ids, 2]))
which(grepl("rain", pathways.13[pathways.13$Id %in% Elk4.13.ids, 2]))
# None!


# Check pathways of Esr2
Esr2.2.5.ids <- which(freq.2.5 == "Esr2")
which(grepl("Alzheimer", pathways.2.5[pathways.2.5$Id %in% Esr2.2.5.ids, 2]))
which(grepl("rain", pathways.2.5[pathways.2.5$Id %in% Esr2.2.5.ids, 2]))
# None!

Esr2.13.ids <- which(freq.13 == "Esr2")
which(grepl("Alzheimer", pathways.13[pathways.13$Id %in% Esr2.13.ids, 2]))
which(grepl("rain", pathways.13[pathways.13$Id %in% Esr2.13.ids, 2]))
# None!


# Check pathways of Nfe2
Nfe2.2.5.ids <- which(freq.2.5 == "Nfe2")
which(grepl("Alzheimer", pathways.2.5[pathways.2.5$Id %in% Nfe2.2.5.ids, 2]))
which(grepl("rain", pathways.2.5[pathways.2.5$Id %in% Nfe2.2.5.ids, 2]))
# Bingo!

Nfe2.13.ids <- which(freq.13 == "Nfe2")
which(grepl("Alzheimer", pathways.13[pathways.13$Id %in% Nfe2.13.ids, 2]))
which(grepl("rain", pathways.13[pathways.13$Id %in% Nfe2.13.ids, 2]))
# Bingo!


# Check pathways of Nfatc2 (difference)
Nfatc2.2.5.ids <- which(freq.2.5 == "Nfatc2")
which(grepl("Alzheimer", pathways.2.5[pathways.2.5$Id %in% Nfatc2.2.5.ids, 2]))
which(grepl("rain", pathways.2.5[pathways.2.5$Id %in% Nfatc2.2.5.ids, 2]))
# None!

Nfatc2.13.ids <- which(freq.13 == "Nfatc2")
which(grepl("Alzheimer", pathways.13[pathways.13$Id %in% Nfatc2.13.ids, 2]))
which(grepl("rain", pathways.13[pathways.13$Id %in% Nfatc2.13.ids, 2]))
# None!


# Check pathways of Plag1 (difference)
Plag1.2.5.ids <- which(freq.2.5 == "Plag1")
which(grepl("Alzheimer", pathways.2.5[pathways.2.5$Id %in% Plag1.2.5.ids, 2]))
which(grepl("rain", pathways.2.5[pathways.2.5$Id %in% Plag1.2.5.ids, 2]))
# None!

Plag1.13.ids <- which(freq.13 == "Plag1")
which(grepl("Alzheimer", pathways.13[pathways.13$Id %in% Plag1.13.ids, 2]))
which(grepl("rain", pathways.13[pathways.13$Id %in% Plag1.13.ids, 2]))
# None!


# Check pathways of Nkx3-1 (difference)
Nkx3.2.5.ids <- which(freq.2.5 == "Nkx3-1")
which(grepl("Alzheimer", pathways.2.5[pathways.2.5$Id %in% Nkx3.2.5.ids, 2]))
which(grepl("rain", pathways.2.5[pathways.2.5$Id %in% Nkx3.2.5.ids, 2]))
# None!

Nkx3.13.ids <- which(freq.13 == "Nkx3-1")
which(grepl("Alzheimer", pathways.13[pathways.13$Id %in% Nkx3.13.ids, 2]))
which(grepl("rain", pathways.13[pathways.13$Id %in% Nkx3.13.ids, 2]))
# None!


# Cytoscape analysis for Jund
Jund.2.5.ids
dim(pathways.2.5)
dim(pathways.13)
Jund.2.5.pathways <- pathways.2.5[pathways.2.5$Id == "11",]
dim(Jund.2.5.pathways)
unique(Jund.2.5.pathways$Term)
jund.2.5.selected.paths <- 
  unique(Jund.2.5.pathways$Term)[c(3, 8, 9, 13, 37, 39, 41, 52, 57, 58, 65)]
Jund.2.5.pathways[Jund.2.5.pathways$Term %in% jund.2.5.selected.paths,]$P.value

Jund.13.pathways <- pathways.13[pathways.13$Id %in% Jund.13.ids]
dim(Jund.13.pathways)
unique(Jund.13.pathways$Term)

# Gene comparison
intersect(Cores.2.5[[11]]$genes, Cores.13[[15]]$genes)
setdiff(Cores.2.5[[11]]$genes, Cores.13[[15]]$genes)
setdiff(Cores.13[[15]]$genes, Cores.2.5[[11]]$genes)

# Peak comparison
intersect(Cores.2.5[[11]]$peaks, Cores.13[[15]]$peaks)
setdiff(Cores.2.5[[11]]$peaks, Cores.13[[15]]$peaks)
setdiff(Cores.13[[15]]$peaks, Cores.2.5[[11]]$peaks)

# Important genes in AD pathway
AD.genes <- strsplit(Jund.2.5.pathways["1324", "Genes"], split = ";")[[1]] %>% 
  tolower %>% capitalize
AD.genes
intersect(AD.genes, Cores.2.5[[11]]$genes)
intersect(AD.genes, Cores.13[[15]]$genes)


# Linked peaks of Jund eGRN of 2.5 months
links.Jund.2.5 <- link_peaks_to_genes(peak.obj = Cores.2.5[[11]]$peaks, 
                                      gene.obj = Cores.2.5[[11]]$genes, 
                    org = "mm10", distance = "gene")
links.Jund.2.5$gene
links.Jund.2.5[links.Jund.2.5$gene %in% AD.genes]


# Load the objects
obj.2.5 <- qs::qread(paste0(R.dir, "Obj_2.5_months.qsave"))
obj.13 <- qs::qread(paste0(R.dir, "Obj_13+_months.qsave"))
rna.2.5 <- GetAssayData(obj.2.5, "data", "RNA")
rna.13 <- GetAssayData(obj.13, "data", "RNA")
atac.2.5 <- GetAssayData(obj.2.5, "data", "ATAC")
atac.13 <- GetAssayData(obj.13, "data", "ATAC")
rna.2.5[AD.genes, Cores.2.5[[11]]$cells]
atac.2.5[Cores.2.5[[11]]$peaks, Cores.2.5[[11]]$cells]
GRangesToString(links.Jund.2.5) %>% unique

# Showcase the UCSC Genome Browser