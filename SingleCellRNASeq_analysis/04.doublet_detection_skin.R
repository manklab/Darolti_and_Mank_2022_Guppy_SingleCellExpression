
skin_f1 <- subset(x=skin_umap_rename_clusters, subset=sample=="SkinF1")
#Get expected number of doublets (no. cells = 3398, estimated doublet rate ~2.6%)
nExp_skin_f1 <- round(ncol(skin_f1)*0.026)
#pK identification
sweep.res.list_f1 <- paramSweep_v3(skin_f1, PCs=1:20, sct=T)
sweep.stats_f1 <- summarizeSweep(sweep.res.list_f1, GT=F)
bcmvn_f1 <- find.pK(sweep.stats_f1)
ggplot(bcmvn_f1, aes(pK, BCmetric, group=1)) +
	geom_point() +
	geom_line()
pK_f1 <- bcmvn_f1 %>% filter(BCmetric == max(BCmetric)) %>% select(pK)
pK_f1 <- as.numeric(as.character(pK_f1[[1]]))
skin_f1_doubl <- doubletFinder_v3(skin_f1, pN=0.25, pK=pK_f1, nExp=nExp_skin_f1, PCs=1:20, sct=TRUE)
DF.name_skin_f1 = colnames(skin_f1_doubl@meta.data)[grepl("DF.classification", colnames(skin_f1_doubl@meta.data))]
DimPlot(skin_f1_doubl, group.by = DF.name_skin_f1, order=c("Doublet", "Singlet"), cols=c("dodgerblue", "orange"))

skin_f2 <- subset(x=skin_umap_rename_clusters, subset=sample=="SkinF2")
#Get expected number of doublets (no. cells = 1690, estimated doublet rate ~1.3%)
nExp_skin_f2 <- round(ncol(skin_f2)*0.013)
#pK identification
sweep.res.list_f2 <- paramSweep_v3(skin_f2, PCs=1:20, sct=T)
sweep.stats_f2 <- summarizeSweep(sweep.res.list_f2, GT=F)
bcmvn_f2 <- find.pK(sweep.stats_f2)
ggplot(bcmvn_f2, aes(pK, BCmetric, group=1)) +
	geom_point() +
	geom_line()
pK_f2 <- bcmvn_f2 %>% filter(BCmetric == max(BCmetric)) %>% select(pK)
pK_f2 <- as.numeric(as.character(pK_f2[[1]]))
skin_f2_doubl <- doubletFinder_v3(skin_f2, pN=0.25, pK=pK_f2, nExp=nExp_skin_f2, PCs=1:20, sct=TRUE)
DF.name_skin_f2 = colnames(skin_f2_doubl@meta.data)[grepl("DF.classification", colnames(skin_f2_doubl@meta.data))]
DimPlot(skin_f2_doubl, group.by = DF.name_skin_f2, order=c("Doublet", "Singlet"), cols=c("dodgerblue", "orange"))

skin_f3 <- subset(x=skin_umap_rename_clusters, subset=sample=="SkinF3")
#Get expected number of doublets (no. cells = 3866, estimated doublet rate ~3.0%)
nExp_skin_f3 <- round(ncol(skin_f3)*0.030)
#pK identification
sweep.res.list_f3 <- paramSweep_v3(skin_f3, PCs=1:20, sct=T)
sweep.stats_f3 <- summarizeSweep(sweep.res.list_f3, GT=F)
bcmvn_f3 <- find.pK(sweep.stats_f3)
ggplot(bcmvn_f3, aes(pK, BCmetric, group=1)) +
	geom_point() +
	geom_line()
pK_f3 <- bcmvn_f3 %>% filter(BCmetric == max(BCmetric)) %>% select(pK)
pK_f3 <- as.numeric(as.character(pK_f3[[1]]))
skin_f3_doubl <- doubletFinder_v3(skin_f3, pN=0.25, pK=pK_f3, nExp=nExp_skin_f3, PCs=1:20, sct=TRUE)
DF.name_skin_f3 = colnames(skin_f3_doubl@meta.data)[grepl("DF.classification", colnames(skin_f3_doubl@meta.data))]
DimPlot(skin_f3_doubl, group.by = DF.name_skin_f3, order=c("Doublet", "Singlet"), cols=c("dodgerblue", "orange"))

skin_m1 <- subset(x=skin_umap_rename_clusters, subset=sample=="SkinM1")
#Get expected number of doublets (no. cells = 957, estimated doublet rate ~0.7%)
nExp_skin_m1 <- round(ncol(skin_m1)*0.007)
#pK identification
sweep.res.list_m1 <- paramSweep_v3(skin_m1, PCs=1:20, sct=T)
sweep.stats_m1 <- summarizeSweep(sweep.res.list_m1, GT=F)
bcmvn_m1 <- find.pK(sweep.stats_m1)
ggplot(bcmvn_m1, aes(pK, BCmetric, group=1)) +
	geom_point() +
	geom_line()
pK_m1 <- bcmvn_m1 %>% filter(BCmetric == max(BCmetric)) %>% select(pK)
pK_m1 <- as.numeric(as.character(pK_m1[[1]]))
skin_m1_doubl <- doubletFinder_v3(skin_m1, pN=0.25, pK=pK_m1, nExp=nExp_skin_m1, PCs=1:20, sct=TRUE)
DF.name_skin_m1 = colnames(skin_m1_doubl@meta.data)[grepl("DF.classification", colnames(skin_m1_doubl@meta.data))]
DimPlot(skin_m1_doubl, group.by = DF.name_skin_m1, order=c("Doublet", "Singlet"), cols=c("dodgerblue", "orange"))

skin_m2 <- subset(x=skin_umap_rename_clusters, subset=sample=="SkinM2")
#Get expected number of doublets (no. cells = 552, estimated doublet rate ~0.4%)
nExp_skin_m2 <- round(ncol(skin_m2)*0.004)
#pK identification
sweep.res.list_m2 <- paramSweep_v3(skin_m2, PCs=1:20, sct=T)
sweep.stats_m2 <- summarizeSweep(sweep.res.list_m2, GT=F)
bcmvn_m2 <- find.pK(sweep.stats_m2)
ggplot(bcmvn_m2, aes(pK, BCmetric, group=1)) +
	geom_point() +
	geom_line()
pK_m2 <- bcmvn_m2 %>% filter(BCmetric == max(BCmetric)) %>% select(pK)
pK_m2 <- as.numeric(as.character(pK_m2[[1]]))
skin_m2_doubl <- doubletFinder_v3(skin_m2, pN=0.25, pK=pK_m2, nExp=nExp_skin_m2, PCs=1:20, sct=TRUE)
DF.name_skin_m2 = colnames(skin_m2_doubl@meta.data)[grepl("DF.classification", colnames(skin_m2_doubl@meta.data))]
DimPlot(skin_m2_doubl, group.by = DF.name_skin_m2, order=c("Doublet", "Singlet"), cols=c("dodgerblue", "orange"))

skin_m3 <- subset(x=skin_umap_rename_clusters, subset=sample=="SkinM3")
#Get expected number of doublets (no. cells = 2334, estimated doublet rate ~1.8%)
nExp_skin_m3 <- round(ncol(skin_m3)*0.018)
#pK identification
sweep.res.list_m3 <- paramSweep_v3(skin_m3, PCs=1:20, sct=T)
sweep.stats_m3 <- summarizeSweep(sweep.res.list_m3, GT=F)
bcmvn_m3 <- find.pK(sweep.stats_m3)
ggplot(bcmvn_m3, aes(pK, BCmetric, group=1)) +
	geom_point() +
	geom_line()
pK_m3 <- bcmvn_m3 %>% filter(BCmetric == max(BCmetric)) %>% select(pK)
pK_m3 <- as.numeric(as.character(pK_m3[[1]]))
skin_m3_doubl <- doubletFinder_v3(skin_m3, pN=0.25, pK=pK_m3, nExp=nExp_skin_m3, PCs=1:20, sct=TRUE)
DF.name_skin_m3 = colnames(skin_m3_doubl@meta.data)[grepl("DF.classification", colnames(skin_m3_doubl@meta.data))]
DimPlot(skin_m3_doubl, group.by = DF.name_skin_m3, order=c("Doublet", "Singlet"), cols=c("dodgerblue", "orange"))


DimPlot(skin_f1_doubl, group.by = DF.name_skin_f1, order=c("Doublet", "Singlet"), cols=c("dodgerblue", "orange")) + DimPlot(skin_f2_doubl, group.by = DF.name_skin_f2, order=c("Doublet", "Singlet"), cols=c("dodgerblue", "orange")) + DimPlot(skin_f3_doubl, group.by = DF.name_skin_f3, order=c("Doublet", "Singlet"), cols=c("dodgerblue", "orange")) + DimPlot(skin_m1_doubl, group.by = DF.name_skin_m1, order=c("Doublet", "Singlet"), cols=c("dodgerblue", "orange")) + DimPlot(skin_m2_doubl, group.by = DF.name_skin_m2, order=c("Doublet", "Singlet"), cols=c("dodgerblue", "orange")) + DimPlot(skin_m3_doubl, group.by = DF.name_skin_m3, order=c("Doublet", "Singlet"), cols=c("dodgerblue", "orange"))


# Remove doublets
skin_f1_remove <- subset(x=skin_f1_doubl, subset=DF.classifications_0.25_0.01_85=="Doublet")
skin_f2_remove <- subset(x=skin_f2_doubl, subset=DF.classifications_0.25_0.02_20=="Doublet")
skin_f3_remove <- subset(x=skin_f3_doubl, subset=DF.classifications_0.25_0.1_116=="Doublet")
skin_m1_remove <- subset(x=skin_m1_doubl, subset=DF.classifications_0.25_0.13_6=="Doublet")
skin_m2_remove <- subset(x=skin_m2_doubl, subset=DF.classifications_0.25_0.04_2=="Doublet")
skin_m3_remove <- subset(x=skin_m3_doubl, subset=DF.classifications_0.25_0.28_41=="Doublet")

skin_umap_rename_clusters_dbrmv <- skin_umap_rename_clusters[,!colnames(skin_umap_rename_clusters)%in%colnames(skin_f1_remove)]
skin_umap_rename_clusters_dbrmv <- skin_umap_rename_clusters_dbrmv[,!colnames(skin_umap_rename_clusters_dbrmv)%in%colnames(skin_f2_remove)]
skin_umap_rename_clusters_dbrmv <- skin_umap_rename_clusters_dbrmv[,!colnames(skin_umap_rename_clusters_dbrmv)%in%colnames(skin_f3_remove)]
skin_umap_rename_clusters_dbrmv <- skin_umap_rename_clusters_dbrmv[,!colnames(skin_umap_rename_clusters_dbrmv)%in%colnames(skin_m1_remove)]
skin_umap_rename_clusters_dbrmv <- skin_umap_rename_clusters_dbrmv[,!colnames(skin_umap_rename_clusters_dbrmv)%in%colnames(skin_m2_remove)]
skin_umap_rename_clusters_dbrmv <- skin_umap_rename_clusters_dbrmv[,!colnames(skin_umap_rename_clusters_dbrmv)%in%colnames(skin_m3_remove)]

before <- DimPlot(skin_umap_rename_clusters, reduction="umap", label=TRUE, label.size=2.5, repel=T, label.box=T,raster=FALSE, pt.size=0.3)
after <- DimPlot(skin_umap_rename_clusters_dbrmv, reduction="umap", label=TRUE, label.size=2.5, repel=T, label.box=T,raster=FALSE, pt.size=0.3)
before + after

save(skin_umap_rename_clusters_dbrmv, file="skin_umap_rename_clusters_dbrmv.RData")

