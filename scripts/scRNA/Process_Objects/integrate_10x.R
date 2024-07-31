# Load libraries
library(Seurat)
library(dplyr)
library(Matrix)

# Function to create a Seurat object from a directory
create_seurat_object <- function(directory) {
  project_name <- basename(directory)  
  # Read the data
  print(paste0('Creating Seurat Object for: ',project_name))

 #READ 10X Data
  data <- Read10X(data.dir = directory)
  colnames(data) <- gsub('-1',paste0('-',directory),colnames(data))
  # Create a Seurat object with the directory name as project name
  seurat_object <- CreateSeuratObject(counts = data,project = project_name)
  print(paste0('Finished on Sample: ',project_name))
  # Add features as row names
  return(seurat_object)
}

# Directory containing subdirectories with the data files
root_directory <- "."
subdirectories <- list.dirs(root_directory, full.names = TRUE, recursive = FALSE)

# Create Seurat objects and merge them
seurat_objects <- lapply(subdirectories, create_seurat_object)

print('Filtering')
seurat_objects <- lapply(seurat_objects,function(x){
					x[['percent.mt']] <- PercentageFeatureSet(x,pattern='^MT-')
					x[['percent.rb']] <- PercentageFeatureSet(x,pattern = "^RP[SL]")
					x <- subset(x,subset = nCount_RNA > 500 & nFeature_RNA > 200 & nCount_RNA < 25000 & nFeature_RNA < 6000 & percent.mt < 12 & percent.rb < 50)})
print('SCTRANSFORM')
seurat_objects <- lapply(seurat_objects, SCTransform, verbose = FALSE)
# Finding integration anchors
print('Finding Anchors')
features <- SelectIntegrationFeatures(object.list = seurat_objects, nfeatures = 3000)
seurat_objects <- PrepSCTIntegration(object.list = seurat_objects, anchor.features = features)
anchors <- FindIntegrationAnchors(
  object.list = seurat_objects,
  normalization.method = "SCT",  # Ensure SCT normalization method is specified
  anchor.features = features,        # Number of features to use
  dims = 1:20,                   # Dimensions to use
  k.filter = 30,                 # Filtering parameter
  k.score = 30                   # Scoring parameter
)
print('Integrating Data')
# Integrating data
integrated_data <- IntegrateData(
  anchorset = anchors,
  dims = 1:20
)
print('Scaling and Dim Reductin')
integrated_data <- ScaleData(integrated_data,verbose=FALSE)
# Run PCA
integrated_data <- RunPCA(integrated_data, verbose = FALSE)

# Find neighbors
integrated_data <- FindNeighbors(integrated_data, dims = 1:20)

# Run UMAP
integrated_data <- RunUMAP(integrated_data, dims = 1:20)
# Optional: Save the final Seurat object
saveRDS(integrated_data, file = "integrated_data.rds")

