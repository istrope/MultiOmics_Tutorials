library(Seurat)
library(dplyr)

# Assuming 'data' is your raw count matrix of scRNA-seq data
fariba <- readRDS('seurat_10x.rds')

# Add percent of mitochondrial and ribosomal genes as metadata columns
mito_genes <- grep(pattern = "^MT-", x = rownames(fariba@assays$RNA@counts), value = TRUE)
ribo_genes <- grep(pattern = "^RP[SL]", x = rownames(fariba@assays$RNA@counts), value = TRUE)

fariba[["percent.mt"]] <- PercentageFeatureSet(fariba, features = mito_genes)
fariba[["percent.rb"]] <- PercentageFeatureSet(fariba, features = ribo_genes)

print('Filtering Samples')
# Filtering cells based on provided criteria
fariba <- subset(fariba, subset = nCount_RNA > 500 & nFeature_RNA > 200 & nCount_RNA < 25000 & nFeature_RNA < 6000 & percent.mt < 12 & percent.rb < 50)

print('SCTTransform')
# Normalization and variance stabilization
fariba <- SCTransform(fariba, vars.to.regress = c("percent.mt", "percent.rb"), verbose = FALSE,variable.features.n = 2000)

print('Dimensionality Reduction')
# Dimensional reduction and clustering
fariba <- RunPCA(fariba, verbose = FALSE,features = VariableFeatures(object = fariba))
fariba <- RunUMAP(fariba, dims = 1:15)


library(ggplot2)
# Assuming you have already performed PCA, UMAP, and clustering as in the previous examples
umapplot <- DimPlot(fariba, reduction = "umap",group.by=c('orig.ident'), label = TRUE)
print('Saving Plots and object')
# Save the UMAP plot to a file
ggsave(filename = "UMAPPlot.png", plot = umapplot, width = 8, height = 6, dpi = 300)
saveRDS(fariba,'fariba.rds')
