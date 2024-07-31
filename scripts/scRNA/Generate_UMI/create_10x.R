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
final_seurat_object <- Reduce(function(x, y) merge(x, y), seurat_objects)

# Optional: Save the final Seurat object
saveRDS(final_seurat_object, file = "seurat_10x.rds")

