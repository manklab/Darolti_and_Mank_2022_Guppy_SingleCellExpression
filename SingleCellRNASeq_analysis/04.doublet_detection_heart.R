
heart_f1 <- subset(x=heart_umap_rename_clusters, subset=sample=="HeartF1")
#Get expected number of doublets (no. cells = 7457, estimated doublet rate ~5.6%)
nExp_heart_f1 <- round(ncol(heart_f1)*0.056)
#pK identification
sweep.res.list_f1 <- paramSweep_v3(heart_f1, PCs=1:15, sct=T)
sweep.stats_f1 <- summarizeSweep(sweep.res.list_f1, GT=F)
bcmvn_f1 <- find.pK(sweep.stats_f1)
ggplot(bcmvn_f1, aes(pK, BCmetric, group=1)) +
	geom_point() +
	geom_line()
pK_f1 <- bcmvn_f1 %>% filter(BCmetric == max(BCmetric)) %>% select(pK)
pK_f1 <- as.numeric(as.character(pK_f1[[1]]))
heart_f1_doubl <- doubletFinder_v3(heart_f1, pN=0.25, pK=pK_f1, nExp=nExp_heart_f1, PCs=1:15, sct=TRUE)
DF.name_heart_f1 = colnames(heart_f1_doubl@meta.data)[grepl("DF.classification", colnames(heart_f1_doubl@meta.data))]
DimPlot(heart_f1_doubl, group.by = DF.name_heart_f1, order=c("Doublet", "Singlet"), cols=c("dodgerblue", "orange"))

heart_f2 <- subset(x=heart_umap_rename_clusters, subset=sample=="HeartF2")
#Get expected number of doublets (no. cells = 10007, estimated doublet rate ~7.6%)
nExp_heart_f2 <- round(ncol(heart_f2)*0.076)
#pK identification
sweep.res.list_f2 <- paramSweep_v3(heart_f2, PCs=1:15, sct=T)
sweep.stats_f2 <- summarizeSweep(sweep.res.list_f2, GT=F)
bcmvn_f2 <- find.pK(sweep.stats_f2)
ggplot(bcmvn_f2, aes(pK, BCmetric, group=1)) +
	geom_point() +
	geom_line()
pK_f2 <- bcmvn_f2 %>% filter(BCmetric == max(BCmetric)) %>% select(pK)
pK_f2 <- as.numeric(as.character(pK_f2[[1]]))
heart_f2_doubl <- doubletFinder_v3(heart_f2, pN=0.25, pK=pK_f2, nExp=nExp_heart_f2, PCs=1:15, sct=TRUE)
DF.name_heart_f2 = colnames(heart_f2_doubl@meta.data)[grepl("DF.classification", colnames(heart_f2_doubl@meta.data))]
DimPlot(heart_f2_doubl, group.by = DF.name_heart_f2, order=c("Doublet", "Singlet"), cols=c("dodgerblue", "orange"))

heart_f3 <- subset(x=heart_umap_rename_clusters, subset=sample=="HeartF3")
#Get expected number of doublets (no. cells = 6817, estimated doublet rate ~5.2%)
nExp_heart_f3 <- round(ncol(heart_f3)*0.052)
#pK identification
sweep.res.list_f3 <- paramSweep_v3(heart_f3, PCs=1:15, sct=T)
sweep.stats_f3 <- summarizeSweep(sweep.res.list_f3, GT=F)
bcmvn_f3 <- find.pK(sweep.stats_f3)
ggplot(bcmvn_f3, aes(pK, BCmetric, group=1)) +
	geom_point() +
	geom_line()
pK_f3 <- bcmvn_f3 %>% filter(BCmetric == max(BCmetric)) %>% select(pK)
pK_f3 <- as.numeric(as.character(pK_f3[[1]]))
heart_f3_doubl <- doubletFinder_v3(heart_f3, pN=0.25, pK=pK_f3, nExp=nExp_heart_f3, PCs=1:15, sct=TRUE)
DF.name_heart_f3 = colnames(heart_f3_doubl@meta.data)[grepl("DF.classification", colnames(heart_f3_doubl@meta.data))]
DimPlot(heart_f3_doubl, group.by = DF.name_heart_f3, order=c("Doublet", "Singlet"), cols=c("dodgerblue", "orange"))

heart_m1 <- subset(x=heart_umap_rename_clusters, subset=sample=="HeartM1")
#Get expected number of doublets (no. cells = 4974, estimated doublet rate ~3.8%)
nExp_heart_m1 <- round(ncol(heart_m1)*0.038)
#pK identification
sweep.res.list_m1 <- paramSweep_v3(heart_m1, PCs=1:15, sct=T)
sweep.stats_m1 <- summarizeSweep(sweep.res.list_m1, GT=F)
bcmvn_m1 <- find.pK(sweep.stats_m1)
ggplot(bcmvn_m1, aes(pK, BCmetric, group=1)) +
	geom_point() +
	geom_line()
pK_m1 <- bcmvn_m1 %>% filter(BCmetric == max(BCmetric)) %>% select(pK)
pK_m1 <- as.numeric(as.character(pK_m1[[1]]))
heart_m1_doubl <- doubletFinder_v3(heart_m1, pN=0.25, pK=pK_m1, nExp=nExp_heart_m1, PCs=1:15, sct=TRUE)
DF.name_heart_m1 = colnames(heart_m1_doubl@meta.data)[grepl("DF.classification", colnames(heart_m1_doubl@meta.data))]
DimPlot(heart_m1_doubl, group.by = DF.name_heart_m1, order=c("Doublet", "Singlet"), cols=c("dodgerblue", "orange"))

heart_m2 <- subset(x=heart_umap_rename_clusters, subset=sample=="HeartM2")
#Get expected number of doublets (no. cells = 2037, estimated doublet rate ~1.5%)
nExp_heart_m2 <- round(ncol(heart_m2)*0.015)
#pK identification
sweep.res.list_m2 <- paramSweep_v3(heart_m2, PCs=1:15, sct=T)
sweep.stats_m2 <- summarizeSweep(sweep.res.list_m2, GT=F)
bcmvn_m2 <- find.pK(sweep.stats_m2)
ggplot(bcmvn_m2, aes(pK, BCmetric, group=1)) +
	geom_point() +
	geom_line()
pK_m2 <- bcmvn_m2 %>% filter(BCmetric == max(BCmetric)) %>% select(pK)
pK_m2 <- as.numeric(as.character(pK_m2[[1]]))
heart_m2_doubl <- doubletFinder_v3(heart_m2, pN=0.25, pK=pK_m2, nExp=nExp_heart_m2, PCs=1:15, sct=TRUE)
DF.name_heart_m2 = colnames(heart_m2_doubl@meta.data)[grepl("DF.classification", colnames(heart_m2_doubl@meta.data))]
DimPlot(heart_m2_doubl, group.by = DF.name_heart_m2, order=c("Doublet", "Singlet"), cols=c("dodgerblue", "orange"))

heart_m3 <- subset(x=heart_umap_rename_clusters, subset=sample=="HeartM3")
#Get expected number of doublets (no. cells = 1426, estimated doublet rate ~1.1%)
nExp_heart_m3 <- round(ncol(heart_m3)*0.011)
#pK identification
sweep.res.list_m3 <- paramSweep_v3(heart_m3, PCs=1:15, sct=T)
sweep.stats_m3 <- summarizeSweep(sweep.res.list_m3, GT=F)
bcmvn_m3 <- find.pK(sweep.stats_m3)
ggplot(bcmvn_m3, aes(pK, BCmetric, group=1)) +
	geom_point() +
	geom_line()
pK_m3 <- bcmvn_m3 %>% filter(BCmetric == max(BCmetric)) %>% select(pK)
pK_m3 <- as.numeric(as.character(pK_m3[[1]]))
heart_m3_doubl <- doubletFinder_v3(heart_m3, pN=0.25, pK=pK_m3, nExp=nExp_heart_m3, PCs=1:15, sct=TRUE)
DF.name_heart_m3 = colnames(heart_m3_doubl@meta.data)[grepl("DF.classification", colnames(heart_m3_doubl@meta.data))]
DimPlot(heart_m3_doubl, group.by = DF.name_heart_m3, order=c("Doublet", "Singlet"), cols=c("dodgerblue", "orange"))


DimPlot(heart_f1_doubl, group.by = DF.name_heart_f1, order=c("Doublet", "Singlet"), cols=c("dodgerblue", "orange")) + DimPlot(heart_f2_doubl, group.by = DF.name_heart_f2, order=c("Doublet", "Singlet"), cols=c("dodgerblue", "orange")) + DimPlot(heart_f3_doubl, group.by = DF.name_heart_f3, order=c("Doublet", "Singlet"), cols=c("dodgerblue", "orange")) + DimPlot(heart_m1_doubl, group.by = DF.name_heart_m1, order=c("Doublet", "Singlet"), cols=c("dodgerblue", "orange")) + DimPlot(heart_m2_doubl, group.by = DF.name_heart_m2, order=c("Doublet", "Singlet"), cols=c("dodgerblue", "orange")) + DimPlot(heart_m3_doubl, group.by = DF.name_heart_m3, order=c("Doublet", "Singlet"), cols=c("dodgerblue", "orange"))


# Remove doublets
heart_f1_remove <- subset(x=heart_f1_doubl, subset=DF.classifications_0.25_0.03_400=="Doublet")
heart_f2_remove <- subset(x=heart_f2_doubl, subset=DF.classifications_0.25_0.01_737=="Doublet")
heart_f3_remove <- subset(x=heart_f3_doubl, subset=DF.classifications_0.25_0.03_337=="Doublet")
heart_m1_remove <- subset(x=heart_m1_doubl, subset=DF.classifications_0.25_0.01_182=="Doublet")
heart_m2_remove <- subset(x=heart_m2_doubl, subset=DF.classifications_0.25_0.01_18=="Doublet")
heart_m3_remove <- subset(x=heart_m3_doubl, subset=DF.classifications_0.25_0.02_14=="Doublet")

heart_umap_rename_clusters_dbrmv <- heart_umap_rename_clusters[,!colnames(heart_umap_rename_clusters)%in%colnames(heart_f1_remove)]
heart_umap_rename_clusters_dbrmv <- heart_umap_rename_clusters_dbrmv[,!colnames(heart_umap_rename_clusters_dbrmv)%in%colnames(heart_f2_remove)]
heart_umap_rename_clusters_dbrmv <- heart_umap_rename_clusters_dbrmv[,!colnames(heart_umap_rename_clusters_dbrmv)%in%colnames(heart_f3_remove)]
heart_umap_rename_clusters_dbrmv <- heart_umap_rename_clusters_dbrmv[,!colnames(heart_umap_rename_clusters_dbrmv)%in%colnames(heart_m1_remove)]
heart_umap_rename_clusters_dbrmv <- heart_umap_rename_clusters_dbrmv[,!colnames(heart_umap_rename_clusters_dbrmv)%in%colnames(heart_m2_remove)]
heart_umap_rename_clusters_dbrmv <- heart_umap_rename_clusters_dbrmv[,!colnames(heart_umap_rename_clusters_dbrmv)%in%colnames(heart_m3_remove)]

before <- DimPlot(heart_umap_rename_clusters, reduction="umap", label=TRUE, label.size=2.5, repel=T, label.box=T,raster=FALSE, pt.size=0.3)
after <- DimPlot(heart_umap_rename_clusters_dbrmv, reduction="umap", label=TRUE, label.size=2.5, repel=T, label.box=T,raster=FALSE, pt.size=0.3)
before + after

save(heart_umap_rename_clusters_dbrmv, file="heart_umap_rename_clusters_dbrmv.RData")

