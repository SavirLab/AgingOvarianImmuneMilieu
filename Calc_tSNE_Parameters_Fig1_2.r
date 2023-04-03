library(Seurat)
data_dir_old <- 'E:/Tal/10XData/01_Ovary_Old_2_Oct_19'
data_dir_young <- 'E:/Tal/10XData/02_Ovary_Young_26_Nov_19/filtered'
expression_matrix_old <- Read10X(data.dir = data_dir_old)
old = CreateSeuratObject(counts = expression_matrix_old, project = "Old")
expression_matrix_young <- Read10X(data.dir = data_dir_young)
young = CreateSeuratObject(counts = expression_matrix_young, project = "young")
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
DimPlot(all_data, reduction = "tsne", label = "True")



##### Dot Plot for Fig 1 - supp 3
immune.combined_new_cluster_tby <- all_data

cd_genes <- c("S100a8","Itgam","Adgre1","Itgax","Klrb1c",
              "Ccl5","Itga1","Gata3","Il17rb","Ly6c2","Cd3e",
              "Trbc2","Cd8b1","Cd4","Cd28","Tmem176b","Il7r",
              "Tcrg-C2","Il2ra","Cd19")
DotPlot(object = immune.combined_new_cluster_tby, features = cd_genes)

