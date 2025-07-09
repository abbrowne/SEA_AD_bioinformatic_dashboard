#!/usr/bin/env Rscript
# R script to convert aggregated CSV files to Seurat objects
# Usage: Rscript convert_to_seurat.R <input_dir> <output_dir>

library(Seurat)
library(reticulate)
library(anndata)
library(Signac)
library(EnsDb.Hsapiens.v86)
library(ggplot2)
library(cowplot)
library(biovizBase)


rna_data <- read_h5ad("test_results/aggregated_multiomic/aggregated_rna_data.h5ad")
temp_counts <- t(as.matrix(rna_data$X))  # or "counts" if that exists
temp_meta <- as.data.frame(as.matrix(rna_data$obs))      # cell-level metadata
rna_seurat <- CreateSeuratObject(counts = temp_counts, meta.data = temp_meta)
saveRDS(rna_seurat,file="E:/Projects/Brain/SEA_AD_analysis/test_results/seurat_objects/aggregated_rna_seurat.RDS")
rm(list=c("temp_counts","temp_meta","rna_data","rna_seurat"))

atac_data <- read_h5ad("test_results/aggregated_multiomic/aggregated_atac_data.h5ad")
temp_counts <- t(as.matrix(atac_data$X))
temp_meta <- as.data.frame(as.matrix(atac_data$obs))      # cell-level metadata
atac_seurat <- CreateSeuratObject(counts = temp_counts, meta.data = temp_meta)

peak_names <- rownames(temp_counts)  # Ensure format is like "chr1:1000-1500"
gr <- StringToGRanges(peak_names, sep = c(":", "-"))
chrom_assay <- CreateChromatinAssay(
  counts = temp_counts,
  ranges = gr,
  genome = 'hg38',
  sep = c(":", "-"),
  annotation = NULL  # can add later using EnsDb
)
atac_seurat[["ATAC"]] <- chrom_assay
DefaultAssay(atac_seurat) <- "ATAC"

# ATAC analysis add gene annotation information
annotations <- GetGRangesFromEnsDb(ensdb = EnsDb.Hsapiens.v86)
seqlevelsStyle(annotations) <- "UCSC"
genome(annotations) <- "hg38"
Annotation(atac_seurat) <- annotations
saveRDS(atac_seurat,file="E:/Projects/Brain/SEA_AD_analysis/test_results/seurat_objects/aggregated_atac_seurat.RDS")
rm(list=c("temp_counts","temp_meta","atac_data","atac_seurat"))

pbmc.rna <- readRDS("test_results/seurat_objects/aggregated_rna_seurat.RDS")
pbmc.atac <- readRDS("test_results/seurat_objects/aggregated_atac_seurat.RDS")

pbmc.rna <- NormalizeData(pbmc.rna)
pbmc.rna <- FindVariableFeatures(pbmc.rna)
pbmc.rna <- ScaleData(pbmc.rna)
pbmc.rna <- RunPCA(pbmc.rna)
pbmc.rna <- RunUMAP(pbmc.rna, dims = 1:30)

pbmc.atac <- RunTFIDF(pbmc.atac)
pbmc.atac <- FindTopFeatures(pbmc.atac, min.cutoff = "q0")
pbmc.atac <- RunSVD(pbmc.atac)
pbmc.atac <- RunUMAP(pbmc.atac, reduction = "lsi", dims = 2:30, reduction.name = "umap.atac", reduction.key = "atacUMAP_")

p1 <- DimPlot(pbmc.rna, group.by = "Subclass", label = TRUE) + NoLegend() + ggtitle("RNA")
p2 <- DimPlot(pbmc.atac, group.by = "Subclass", label = FALSE) + NoLegend() + ggtitle("ATAC")
plot <- (p1 + p2) & xlab("UMAP 1") & ylab("UMAP 2") & theme(axis.title = element_text(size = 18))
plot






