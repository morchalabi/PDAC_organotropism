# This script plots Fig. 3a (UMAP projection of cancer compartment).

library(Seurat)
options(Seurat.object.assay.version = "v3")
library(SeuratObject)
library(ggplot2)

# Reading in cancer compartment ####

load(file = 'compartments/cancer.RData')
DefaultAssay(s_objs) = 'RNA'

# Plotting UMAP ####

DimPlot(s_objs, group.by = 'cell_type', label = T, label.size = 5, pt.size = 1, repel = T) & NoLegend()
ggsave(filename = 'a.pdf', device = 'pdf', width = 7, height = 7)
