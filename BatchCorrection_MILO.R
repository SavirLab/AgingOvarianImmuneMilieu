library(dplyr)
library(Seurat)
library(SeuratData)
library(patchwork)
library(edgeR)
library(limma)
library(scater)
library(miloR)
library(statmod)
library(MultinomialCI)

# Loading the data
data_dir_old <- '/Users/eliel/Library/CloudStorage/OneDrive-Technion/R Code/10XData/01_Ovary_Old_2_Oct_19'
data_dir_young <- '/Users/eliel/Library/CloudStorage/OneDrive-Technion/R Code/10XData/02_Ovary_Young_26_Nov_19/filtered'
expression_matrix_old <- Read10X(data.dir = data_dir_old)
old = CreateSeuratObject(counts = expression_matrix_old, project = "Old")
expression_matrix_young <- Read10X(data.dir = data_dir_young)
young = CreateSeuratObject(counts = expression_matrix_young, project = "young")

# Merge old and young
all_data <- merge(young, y = old, add.cell.ids = c("Young", "Old"), project = "All10X")

### All data analysis ###
all_data[["percent.mt"]] <- PercentageFeatureSet(all_data, pattern = "^mt-")
VlnPlot(all_data, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
all_data <- subset(all_data, subset = nFeature_RNA > 200 & nFeature_RNA < 2500 & percent.mt < 10)
all_data <- NormalizeData(all_data, normalization.method = "LogNormalize", scale.factor = 10000)
all_data <- FindVariableFeatures(all_data, selection.method = "vst", nfeatures = 2000)
all.genes <- rownames(all_data)
all_data <- ScaleData(all_data, features = all.genes)
all_data <- RunPCA(all_data, features = VariableFeatures(object = all_data))
all_data <- FindNeighbors(all_data, dims = 1:15)
all_data <- FindClusters(all_data, resolution = 0.55)

all_data <- RunTSNE(all_data, dims = 1:15)
all_data <- RunUMAP(all_data, dims = 1:15)
DimPlot(all_data, reduction = "tsne",split.by = "orig.ident",label = "TRUE")
DimPlot(all_data, reduction = "umap")

###############################################


# Write to csv 
data_to_write_out <- as.data.frame(as.matrix(all_data@reductions$tsne))
fwrite(x = data_to_write_out, file = "outfile.csv")
write.csv(x=all_data@reductions$tsne@cell.embeddings,file ='/Users/eliel/Library/CloudStorage/OneDrive-Technion/R Code/EA/tsnecorr_wo_batch.csv')
write.csv(x=all_data@meta.data$seurat_clusters,file ='/Users/eliel/Library/CloudStorage/OneDrive-Technion/R Code/EA/tsnecorr_wo_batch_clusters.csv')

# only good clusters
# clustred_data <- subset(all_data,idents=c(13,14) , invert = TRUE)
# DimPlot(clustred_data, reduction = "tsne",split.by = "orig.ident",label = "TRUE")


######################### Batch correction #########################

age.list <- SplitObject(all_data, split.by = "orig.ident")
features <- SelectIntegrationFeatures(object.list = age.list)
immune.anchors <- FindIntegrationAnchors(object.list = age.list, anchor.features = features)
immune.combined <- IntegrateData(anchorset = immune.anchors)
immune.combined <- ScaleData(immune.combined, features = all.genes)
immune.combined <- RunPCA(immune.combined, features = VariableFeatures(object = all_data))
immune.combined <- FindNeighbors(immune.combined, dims = 1:15)
immune.combined <- FindClusters(immune.combined, resolution = 0.55)
immune.combined <- RunTSNE(immune.combined, dims = 1:15)
                           
DimPlot(immune.combined, reduction = "tsne",split.by = "orig.ident",label = "TRUE")

# write.csv(x=immune.combined@reductions$tsne@cell.embeddings,file ='/Users/eliel/Library/CloudStorage/OneDrive-Technion/R Code/EA/tsnecorr_w_batch.csv')
# write.csv(x=immune.combined@meta.data$seurat_clusters,file ='/Users/eliel/Library/CloudStorage/OneDrive-Technion/R Code/EA/tsnecorr_w_batch_clusters.csv')
# write.csv(x=immune.combined@meta.data$orig.ident,file ='/Users/eliel/Library/CloudStorage/OneDrive-Technion/R Code/EA/age_factor.csv')
# write.csv(x=immune.combined@meta.data$seurat_clusters,file ='/Users/eliel/Library/CloudStorage/OneDrive-Technion/R Code/EA/tsnecorr_w_batch_clusters.csv')

immune.combined_new_cluster <- immune.combined

# load new clustering idx file and update values in immune seurat object
new_clustering_value <- read.csv('/Users/eliel/Library/CloudStorage/OneDrive-Technion/R Code/EA/new_clusters_idx_with_cord.csv')
for (idx in 1:length(new_clustering_value$idx_aligned)) {
  immune.combined_new_cluster@meta.data$seurat_clusters[idx] = new_clustering_value$idx_aligned[idx]
  immune.combined_new_cluster@active.ident[idx] = new_clustering_value$idx_aligned[idx]
}
immune.combined_new_cluster <- subset(immune.combined_new_cluster, ident=c(15) , invert = TRUE)

DimPlot(immune.combined_new_cluster, reduction = "tsne",split.by = "orig.ident",label = "TRUE")
DimPlot(immune.combined_new_cluster, reduction = "tsne",label = "TRUE")
###########################################################################
# DotPlot for all clusters (equivalent to Fig. 1C)
immune.combined_new_cluster_tby <- all_data
new_clustering_value <- read.csv('C:/Users/SavirLab/Technion/Yoni Savir - TalBenYakov/R Code/EA/new_clusters_idx_with_cord.csv')

for (idx in 1:length(new_clustering_value$idx_aligned)) {
  all_data@meta.data$seurat_clusters[idx] = new_clustering_value$idx_aligned[idx]
  all_data@active.ident[idx] = new_clustering_value$idx_aligned[idx]
}
all_data <- subset(immune.combined_new_cluster_tby, ident=c(14) , invert = TRUE)

DimPlot(all_data, reduction = "tsne",split.by = "orig.ident",label = "FALSE")
DimPlot(all_data, reduction = "tsne",label = "FALSE")

new.cluster.ids <- c("ILC1", "DNT", "CD8+CD4 T", "NKT", "Neutrophils", "NK",
                     "Macrophages-1", "B cells", "Dendritic cells-1","Macrophages-2",
                     "Dendritic cells-2","ILC2","ILC3","X","X")
names(new.cluster.ids) <- levels(immune.combined_new_cluster_tby)
immune.combined_new_cluster_tby <- RenameIdents(immune.combined_new_cluster_tby, new.cluster.ids)
cd_genes <- c("S100a8","Itgam","Adgre1","Itgax","Klrb1c",
              "Ccl5","Itga1","Gata3","Il17rb","Ly6c2","Cd3e",
              "Trbc2","Cd8b1","Cd4","Cd28","Tmem176b","Il7r",
              "Tcrg-C2","Il2ra","Cd19")
DotPlot(object = immune.combined_new_cluster_tby, features = cd_genes)

###########################################################################


# Perform Milo on the immune.combined_new_cluster (after batch correction) - EA
immune.combined_bc_sce <- as.SingleCellExperiment(immune.combined_new_cluster) # make a sc object
immune.combined_bc_milo <- Milo(immune.combined_bc_sce)

k = 13
d = 30
immune.combined_bc_milo <- buildGraph(immune.combined_bc_milo, k = k, d = d) #k=13, d=30: n=330
immune.combined_bc_milo <- makeNhoods(immune.combined_bc_milo, prop = 0.1, k = k, d = d, refined = TRUE)
plotNhoodSizeHist(immune.combined_bc_milo)

immune.combined_bc_milo@colData$Samples <- c(rep(1, 1657), rep(2, 1658), rep(3, 2738), rep(4, 2738))
immune.combined_bc_milo <- countCells(immune.combined_bc_milo, meta.data = as.data.frame(colData(immune.combined_bc_milo)), sample="Samples")
immune.combined_bc_milo <- calcNhoodDistance(immune.combined_bc_milo, d=60)
immune.combined_bc_milo <- buildNhoodGraph(immune.combined_bc_milo)
head(nhoodCounts(immune.combined_bc_milo))

age_design <- data.frame(colData(immune.combined_bc_milo))[c("Samples", "orig.ident")]
age_design <- distinct(age_design)
rownames(age_design) <- age_design$Samples

da_results <- testNhoods(immune.combined_bc_milo, design = ~ orig.ident, design.df = age_design)
head(da_results)

# plot milo da results
tsne_plot <-plotReducedDim(immune.combined_bc_sce, colour_by="orig.ident", dimred = "TSNE") 

ggplot(da_results, aes(logFC, -log10(SpatialFDR))) + 
  geom_point() +
  geom_hline(yintercept = 1) ## Mark significance threshold (10% FDR)

nh_graph_pl <- plotNhoodGraphDA(immune.combined_bc_milo, da_results, layout="TSNE",alpha=0.1) 

tsne_plot + nh_graph_pl +
  plot_layout(guides="collect")

# add cluster data
da_results <- annotateNhoods(immune.combined_bc_milo, da_results, coldata_col = "ident")
head(da_results)
da_results <- groupNhoods(immune.combined_bc_milo, da_results, max.lfc.delta = 5)
plotNhoodGroups(immune.combined_bc_milo, da_results, layout="TSNE")
plotDAbeeswarm(da_results, group.by = "ident")
da_results$logFC <- (da_results$logFC * (-1) - log(5476/3315, 2))
ident_name <- da_results$ident
ident_name <- plyr::mapvalues(ident_name, from=c("1", "2", "3", "4", "5", "6", "7", "8", "9", "10", "11", "12", "13"), 
                        to=c("Neutrophils", "Macrophages", "cDC2", "cDC1", "NK cells", "ILC1", "ILC2", "ILC3", "NKT cells",
                             "CD8+ T", "CD4+ T", "CD8- CD4- cells", "B cells"))
da_results$ident_name <- ident_name
plotDAbeeswarm(da_results, group.by = "ident_name")
#####################


write.csv(x=da_results,file ='/Users/eliel/Library/CloudStorage/OneDrive-Technion/R Code/EA/da_results.csv')
write.csv(x=immune.combined_bc_milo@nhoods, file='/Users/eliel/Library/CloudStorage/OneDrive-Technion/R Code/EA/bc_milo.csv')

# do some changes in da_results
# first we want to inverse the LFC, for examples ident '12' should be with positive LFC
# second we want to normalize the number 

da_results_fc_correction <- da_results
da_results_fc_correction$logFC <- (da_results_fc_correction$logFC * (-1) )
plotDAbeeswarm(da_results_fc_correction, group.by = "ident")

#load new da after normalize data
da_results_norm <- read.csv('/Users/eliel/Library/CloudStorage/OneDrive-Technion/R Code/EA/da_results_cell_591nhoods_after_correction_v2.csv')
da_results_norm <- annotateNhoods(immune.combined_bc_milo, da_results_norm, coldata_col = "ident")
plotDAbeeswarm(da_results_norm, group.by = "ident")

ident_name <- da_results_norm$ident
ident_name <- plyr::mapvalues(ident_name, from=c("1", "2", "3", "4", "5", "6", "7", "8", "9", "10", "11", "12", "13"), 
                        to=c("Neutrophils", "Macrophages", "cDC2", "cDC1", "NK cells", "ILC1", "ILC2", "ILC3", 
                             "NKT cells", "CD8+ T", "CD4+ T", "CD8- CD4- cells", "B cells"))
da_results_norm$ident_name <- ident_name
plotDAbeeswarm(da_results_norm, group.by = "ident_name")
################################################################################


# multiomial error
m = multinomialCI(c(1308,1990,462,124,134,181,258,227,146,152,47,71), 0.05)


# Perform Milo on the clustred data - OLD

clustred_data_sce <- as.SingleCellExperiment(clustred_data) # make a sc object


clustred_data_milo <- Milo(clustred_data_sce) # make a milo obkect



clustred_data_milo <- buildGraph(clustred_data_milo, k = 13, d = 30)

clustred_data_milo <- makeNhoods(clustred_data_milo, prop = 0.1, k = 13, d=30, refined = TRUE)

plotNhoodSizeHist(clustred_data_milo)

clustred_data_milo@colData$Samples <- c(rep(1, 1657), rep(2, 1658), rep(3, 2738), rep(4, 2738))
clustred_data_milo <- countCells(clustred_data_milo, meta.data = as.data.frame(colData(clustred_data_milo)), sample="Samples")

clustred_data_milo <- calcNhoodDistance(clustred_data_milo, d=30, reduced.dim = "pca.corrected")

head(nhoodCounts(clustred_data_milo))

###############
age_design <- data.frame(colData(clustred_data_milo))[c("Samples", "orig.ident")]
age_design <- distinct(age_design)
rownames(age_design) <- age_design$Samples

da_results <- testNhoods(clustred_data_milo, design = ~ orig.ident, design.df = age_design)

clustred_data_milo <- buildNhoodGraph(clustred_data_milo)

head(da_results)

# plot milo da results
tsne_plot <-plotReducedDim(clustred_data_sce, colour_by="orig.ident", dimred = "TSNE") 

ggplot(da_results, aes(logFC, -log10(SpatialFDR))) + 
  geom_point() +
  geom_hline(yintercept = 1) ## Mark significance threshold (10% FDR)

nh_graph_pl <- plotNhoodGraphDA(clustred_data_milo, da_results, layout="TSNE",alpha=0.1) 


tsne_plot + nh_graph_pl +
  plot_layout(guides="collect")

# add cluster data

da_results <- annotateNhoods(clustred_data_milo, da_results, coldata_col = "ident")

head(da_results)

da_results <- groupNhoods(clustred_data_milo, da_results, max.lfc.delta = 5)

plotNhoodGroups(clustred_data_milo, da_results, layout="TSNE")

plotDAbeeswarm(da_results, group.by = "ident")

write.csv(x=da_results,file ='C:/Users/SavirLab/OneDrive - Technion/TalBenYakov/R Code/da_results.csv')

# multiomial error

m = multinomialCI(c(1308,1990,462,124,134,181,258,227,146,152,47,71), 0.05)

