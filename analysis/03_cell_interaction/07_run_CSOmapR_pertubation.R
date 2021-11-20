suppressPackageStartupMessages({
  library(scater)
  library(scran)
  library(CSOmapR)
  library(logger)
  library(tidyverse)
})
rm(list = ls())
######################
# Smartseq2 CD45 All
######################
sce <- readRDS("../../data/expression/sce/sce_Smartseq2_scHCC-CD45_featureCounts_qc_clustered_analysed.rds")
assay(sce, "tpm") <- calculateTPM(sce, rowData(sce)$Length)

tpmdata <- tpm(sce)
tpmdata <- as.data.frame(as.matrix(tpmdata))

caf_cells <- colnames(sce[, grepl("^CAF", sce$celltype)])

#####
#perform CSOmap

log_info("CSOmap started ...")
TPM <- tpmdata
LR <- read.table('/raid1/zhangqiming/01Project_HCC/02Analysis/CSOmap/CSOmap.R/data/demo/LR_pairs.txt',
                 header = F)
colnames(LR) <- c('lignad','receptor','weight')
LR[,1:2] <- sapply(LR[,1:2], as.character)

affinityMat = getAffinityMat(tpmdata, LR, verbose = T)
log_info("TSNE step")

TPM['HGF',caf_cells] <- 0
tpmdata['HGF',caf_cells] <- 0

coords_res = runExactTSNE_R(
  X = affinityMat,
  no_dims = 3,
  max_iter = 1000,
  verbose = T
)

coords = coords_res$Y
rownames(coords) <- colnames(TPM)
colnames(coords) <- c('x', 'y', 'z')
coords_tbl = bind_cols(cellName = rownames(coords), as.data.frame(coords))



#signif_results = getSignificance(coords, labels = cellinfo_tbl$labels, verbose = T)
#contribution_list = getContribution(TPM, LR, signif_results$detailed_connections)


coords_outdir = paste0("./out/CSOmapR", "/coordinates_20210610_pertubateHGF_in_CAF.txt")
log_info("Writing coordintates to ", coords_outdir)
write_tsv(coords_tbl, path = coords_outdir)




sce <- readRDS("../../data/expression/sce/sce_Smartseq2_scHCC-CD45_featureCounts_qc_clustered_analysed.rds")
assay(sce, "tpm") <- calculateTPM(sce, rowData(sce)$Length)

tpmdata <- tpm(sce)
tpmdata <- as.data.frame(as.matrix(tpmdata))

caf_cells <- colnames(sce[, grepl("^CAF", sce$celltype)])

#####
#perform CSOmap

log_info("CSOmap started ...")
TPM <- tpmdata
LR <- read.table('/raid1/zhangqiming/01Project_HCC/02Analysis/CSOmap/CSOmap.R/data/demo/LR_pairs.txt',
                 header = F)
colnames(LR) <- c('lignad','receptor','weight')
LR[,1:2] <- sapply(LR[,1:2], as.character)

affinityMat = getAffinityMat(tpmdata, LR, verbose = T)
log_info("TSNE step")

TPM['TGFB1',caf_cells] <- 0
tpmdata['TGFB1',caf_cells] <- 0

coords_res = runExactTSNE_R(
  X = affinityMat,
  no_dims = 3,
  max_iter = 1000,
  verbose = T
)

coords = coords_res$Y
rownames(coords) <- colnames(TPM)
colnames(coords) <- c('x', 'y', 'z')
coords_tbl = bind_cols(cellName = rownames(coords), as.data.frame(coords))



#signif_results = getSignificance(coords, labels = cellinfo_tbl$labels, verbose = T)
#contribution_list = getContribution(TPM, LR, signif_results$detailed_connections)


coords_outdir = paste0("./out/CSOmapR", "/coordinates_20210610_pertubateTGFB1_in_CAF.txt")
log_info("Writing coordintates to ", coords_outdir)
write_tsv(coords_tbl, path = coords_outdir)





log_info("CSOmap started ...")
TPM <- tpmdata
LR <- read.table('/raid1/zhangqiming/01Project_HCC/02Analysis/CSOmap/CSOmap.R/data/demo/LR_pairs.txt',
                 header = F)
colnames(LR) <- c('lignad','receptor','weight')
LR[,1:2] <- sapply(LR[,1:2], as.character)

affinityMat = getAffinityMat(tpmdata, LR, verbose = T)
log_info("TSNE step")

TPM['FGB',caf_cells] <- 0
tpmdata['FGB',caf_cells] <- 0

coords_res = runExactTSNE_R(
  X = affinityMat,
  no_dims = 3,
  max_iter = 1000,
  verbose = T
)

coords = coords_res$Y
rownames(coords) <- colnames(TPM)
colnames(coords) <- c('x', 'y', 'z')
coords_tbl = bind_cols(cellName = rownames(coords), as.data.frame(coords))



#signif_results = getSignificance(coords, labels = cellinfo_tbl$labels, verbose = T)
#contribution_list = getContribution(TPM, LR, signif_results$detailed_connections)


coords_outdir = paste0("./out/CSOmapR", "/coordinates_20210610_pertubateFGB_in_CAF.txt")
log_info("Writing coordintates to ", coords_outdir)
write_tsv(coords_tbl, path = coords_outdir)


log_info("CSOmap started ...")
TPM <- tpmdata
LR <- read.table('/raid1/zhangqiming/01Project_HCC/02Analysis/CSOmap/CSOmap.R/data/demo/LR_pairs.txt',
                 header = F)
colnames(LR) <- c('lignad','receptor','weight')
LR[,1:2] <- sapply(LR[,1:2], as.character)

affinityMat = getAffinityMat(tpmdata, LR, verbose = T)
log_info("TSNE step")

TPM['FGG',caf_cells] <- 0
tpmdata['FGG',caf_cells] <- 0

coords_res = runExactTSNE_R(
  X = affinityMat,
  no_dims = 3,
  max_iter = 1000,
  verbose = T
)

coords = coords_res$Y
rownames(coords) <- colnames(TPM)
colnames(coords) <- c('x', 'y', 'z')
coords_tbl = bind_cols(cellName = rownames(coords), as.data.frame(coords))



#signif_results = getSignificance(coords, labels = cellinfo_tbl$labels, verbose = T)
#contribution_list = getContribution(TPM, LR, signif_results$detailed_connections)


coords_outdir = paste0("./out/CSOmapR", "/coordinates_20210610_pertubateFGG_in_CAF.txt")
log_info("Writing coordintates to ", coords_outdir)
write_tsv(coords_tbl, path = coords_outdir)


log_info("CSOmap started ...")
TPM <- tpmdata
LR <- read.table('/raid1/zhangqiming/01Project_HCC/02Analysis/CSOmap/CSOmap.R/data/demo/LR_pairs.txt',
                 header = F)
colnames(LR) <- c('lignad','receptor','weight')
LR[,1:2] <- sapply(LR[,1:2], as.character)

affinityMat = getAffinityMat(tpmdata, LR, verbose = T)
log_info("TSNE step")

TPM['HLA-B',caf_cells] <- 0
tpmdata['HLA-B',caf_cells] <- 0

coords_res = runExactTSNE_R(
  X = affinityMat,
  no_dims = 3,
  max_iter = 1000,
  verbose = T
)

coords = coords_res$Y
rownames(coords) <- colnames(TPM)
colnames(coords) <- c('x', 'y', 'z')
coords_tbl = bind_cols(cellName = rownames(coords), as.data.frame(coords))



#signif_results = getSignificance(coords, labels = cellinfo_tbl$labels, verbose = T)
#contribution_list = getContribution(TPM, LR, signif_results$detailed_connections)


coords_outdir = paste0("./out/CSOmapR", "/coordinates_20210610_pertubateHLA-B_in_CAF.txt")
log_info("Writing coordintates to ", coords_outdir)
write_tsv(coords_tbl, path = coords_outdir)


log_info("CSOmap started ...")
TPM <- tpmdata
LR <- read.table('/raid1/zhangqiming/01Project_HCC/02Analysis/CSOmap/CSOmap.R/data/demo/LR_pairs.txt',
                 header = F)
colnames(LR) <- c('lignad','receptor','weight')
LR[,1:2] <- sapply(LR[,1:2], as.character)

affinityMat = getAffinityMat(tpmdata, LR, verbose = T)
log_info("TSNE step")

TPM['FGA',caf_cells] <- 0
tpmdata['FGA',caf_cells] <- 0

coords_res = runExactTSNE_R(
  X = affinityMat,
  no_dims = 3,
  max_iter = 1000,
  verbose = T
)

coords = coords_res$Y
rownames(coords) <- colnames(TPM)
colnames(coords) <- c('x', 'y', 'z')
coords_tbl = bind_cols(cellName = rownames(coords), as.data.frame(coords))



#signif_results = getSignificance(coords, labels = cellinfo_tbl$labels, verbose = T)
#contribution_list = getContribution(TPM, LR, signif_results$detailed_connections)


coords_outdir = paste0("./out/CSOmapR", "/coordinates_20210610_pertubateFGA_in_CAF.txt")
log_info("Writing coordintates to ", coords_outdir)
write_tsv(coords_tbl, path = coords_outdir)