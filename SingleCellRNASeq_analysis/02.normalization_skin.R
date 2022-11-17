

all_samples <- subset(x=merged_seurat_all_samples, subset=tissue=="Skin")

split_seurat <- SplitObject(all_samples, split.by="sample")

# Normalize with SCTransform and identify the most variable genes
split_seurat$SkinF1 <- SCTransform(split_seurat$SkinF1)
varfeat_F1 <- VariableFeatures(split_seurat$SkinF1)
split_seurat$SkinF2 <- SCTransform(split_seurat$SkinF2)
varfeat_F2 <- VariableFeatures(split_seurat$SkinF2)
split_seuratSkinF3 <- SCTransform(split_seurat$SkinF3)
varfeat_F3 <- VariableFeatures(split_seurat$SkinF3)
split_seurat$SkinM1 <- SCTransform(split_seurat$SkinM1)
varfeat_M1 <- VariableFeatures(split_seurat$SkinM1)
split_seurat$SkinM2 <- SCTransform(split_seurat$SkinM2)
varfeat_M2 <- VariableFeatures(split_seurat$SkinM2)
split_seurat$SkinM3 <- SCTransform(split_seurat$SkinM3)
varfeat_M3 <- VariableFeatures(split_seurat$SkinM3)

intersect_females = intersect(intersect(varfeat_F1, varfeat_F2), varfeat_F3)
intersect_males = intersect(intersect(varfeat_M1, varfeat_M2), varfeat_M3)
intersect_all = intersect(intersect_females, intersect_males)

merged <- merge(x = split_seurat$SkinF1, 
                       y = c(split_seurat$SkinF2, split_seurat$SkinF3, split_seurat$SkinM1, split_seurat$SkinM2, split_seurat$SkinM3),
                       add.cell.id = c("F1", "F2", "F3", "M1", "M2", "M3"), merge.data=TRUE)

# Perform PCA
pca <- RunPCA(object=merged, features=intersect_all)
PCAPlot(pca, split.by="sample", raster=FALSE, pt.size=0.3)

# Determine how many PCs to use
ElbowPlot(pca)
pct <- pca[["pca"]]@stdev / sum(pca[["pca"]]@stdev) *100
cumulat_all <- cumsum(pct)
co1_all <- which(cumulat_all>90 & pct<5)[1]
col2_all <- sort(which((pct[1:length(pct)-1] - pct[2:length(pct)]) > 0.1), decreasing =T)[1] +1
col2_all
[1] 20

# UMAP visualization
umap <- RunUMAP(pca, dims=1:20, reduction="pca")
# Plot by sex
DimPlot(umap, raster=FALSE, pt.size=0.3, order=c("Lib34_filtered_feature_bc_matrix", "Lib33_filtered_feature_bc_matrix", "Lib32_filtered_feature_bc_matrix", "Lib26_filtered_feature_bc_matrix", "Lib19_filtered_feature_bc_matrix", "Lib18_filtered_feature_bc_matrix"), split.by="sex") + ggtitle('Skin')
# Plot by sample
DimPlot(umap, raster=FALSE, pt.size=0.3, order=c("Lib34_filtered_feature_bc_matrix", "Lib33_filtered_feature_bc_matrix", "Lib32_filtered_feature_bc_matrix", "Lib26_filtered_feature_bc_matrix", "Lib19_filtered_feature_bc_matrix", "Lib18_filtered_feature_bc_matrix") + ggtitle('Skin')
# Plot all together
DimPlot(umap, raster=FALSE, pt.size=0.3, order=c("Lib34_filtered_feature_bc_matrix", "Lib33_filtered_feature_bc_matrix", "Lib32_filtered_feature_bc_matrix", "Lib26_filtered_feature_bc_matrix", "Lib19_filtered_feature_bc_matrix", "Lib18_filtered_feature_bc_matrix")) + ggtitle('Skin')

# Determine which resolution to use
umap <- FindNeighbors(object=umap, dims=1:20)
umap <- FindClusters(object=umap, resolution=c(0, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9))
clustree(umap)
Idents(object=umap) <- "SCT_snn_res.0.1"
res01 <- DimPlot(umap, reduction="umap", label=TRUE, label.size=6, raster=FALSE, pt.size=0.3) + ggtitle('res0.1')
Idents(object=umap) <- "SCT_snn_res.0.2"
res02 <- DimPlot(umap, reduction="umap", label=TRUE, label.size=6, raster=FALSE, pt.size=0.3) + ggtitle('res0.2')
Idents(object=umap) <- "SCT_snn_res.0.3"
res03 <- DimPlot(umap, reduction="umap", label=TRUE, label.size=6, raster=FALSE, pt.size=0.3) + ggtitle('res0.3')
...
...
Idents(object=umap) <- "SCT_snn_res.0.9"
res09 <- DimPlot(umap, reduction="umap", label=TRUE, label.size=6, raster=FALSE, pt.size=0.3) + ggtitle('res0.9')

# Visualize clusters
Idents(object=umap) <- "SCT_snn_res.0.4"
# Plot by sex
DimPlot(umap, reduction="umap", label=TRUE, label.size=6, raster=FALSE, pt.size=0.3, split.by="sex") + ggtitle('SCTransform')
# Plot by sample
DimPlot(umap, reduction="umap", label=TRUE, label.size=6, raster=FALSE, pt.size=0.3, split.by="sample") + ggtitle('SCTransform')
# Plot all together
DimPlot(umap, reduction="umap", label=TRUE, label.size=6, raster=FALSE, pt.size=0.3) + ggtitle('SCTransform')

save(pca, file="skin_pca.RData")
save(umap, file="skin_umap.RData")
