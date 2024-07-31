library(Seurat)
library(infercnv)
library(dplyr)
library(tidyr)

print('reading data')
fariba <- readRDS('fariba.rds')

print('Editing Metadata')
metadata <- fariba@meta.data

annotations <- write.table(data.frame(cbind(rownames(metadata),metadata$scpred_no_rejection)),file='inferCNV/annotations.txt',sep='\t',row.names=F,col.names=F)

print('Create Infercnv Object:')
data.infercnv <- CreateInfercnvObject(raw_counts_matrix = as.matrix(fariba@assays$RNA@counts),
				      annotations_file='inferCNV/annotations.txt',
				      delim='\t',
				      gene_order_file='inferCNV/gene_pos.txt',
				      ref_group_names=c('vascular','lymphatic','pericytes','bcells','myeloid','tcells'))

data.infercnv <- infercnv::run(data.infercnv,
			       cutoff=0.1,
			       out_dir = 'inferCNV',
			       cluster_by_groups = F,
			       plot_steps=F,
			       max_centered_threshold = 3,
			       sd_amplifier = 1.3,
			       analysis_mode='samples',
			       cluster_references = T,
			       denoise=T,
			       no_plot=T,
			       HMM=F)
