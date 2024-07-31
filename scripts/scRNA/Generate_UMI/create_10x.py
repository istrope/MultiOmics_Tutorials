import os
import scanpy as sc
from scipy.io import mmread
from scipy.sparse import csr_matrix
import pandas as pd
import anndata as ad
import re
def create_anndata_object(directory):
    project_name = os.path.basename(directory)
    print(f'Creating AnnData object for: {project_name}')

    # Read the 10X data
    X = mmread(os.path.join(directory,'matrix.mtx.gz'))
    X = csr_matrix(X)
    X = X.transpose()
    var = pd.read_csv(os.path.join(directory,'features.tsv.gz'),sep='\t',header=None)
    var.columns = ['gene_ids','feature_name','type']
    var.index = var['gene_ids']
    obs = pd.read_csv(os.path.join(directory,'barcodes.tsv.gz'),sep='\t',header=None)
    obs.columns = ['cell_id']
    obs.index = obs['cell_id']

    adata = ad.AnnData(X = X, obs = obs,var = var)
    # Assign the project name
    adata.obs['project'] = project_name
    adata.obs.index = [re.sub('-1','',x) for x in adata.obs.index.astype(str)]
    adata.obs.index = adata.obs.index.astype(str) + '-' + adata.obs['project'].astype(str)
    file = project_name + '.h5ad'
    print(file)
    adata.write_h5ad(file)
    print(f'Finished on Sample: {project_name}')
    
    return adata

# Directory containing subdirectories with the data files
root_directory = "."
subdirectories = [os.path.join(root_directory, sub_dir) for sub_dir in os.listdir(root_directory) 
                  if os.path.isdir(os.path.join(root_directory, sub_dir))]

# Create AnnData objects and concatenate them
[create_anndata_object(directory) for directory in subdirectories]


final_anndata_object = anndata_objects[0].concatenate(*anndata_objects[1:], batch_key="batch")

# Optional: Save the final AnnData object
final_anndata_object.write("anndata_10x.h5ad")
print("Final AnnData object saved as anndata_10x.h5ad")

