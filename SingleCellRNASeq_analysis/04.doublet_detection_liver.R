
liver_f1 <- subset(x=liver_umap_rename_clusters, subset=sample=="LiverF1")
#Get expected number of doublets (no. cells = 6220, estimated doublet rate ~4.8%)
nExp_liver_f1 <- round(ncol(liver_f1)*0.048)
#pK identification
sweep.res.list_f1 <- paramSweep_v3(liver_f1, PCs=1:17, sct=T)
sweep.stats_f1 <- summarizeSweep(sweep.res.list_f1, GT=F)
bcmvn_f1 <- find.pK(sweep.stats_f1)
ggplot(bcmvn_f1, aes(pK, BCmetric, group=1)) +
	geom_point() +
	geom_line()
pK_f1 <- bcmvn_f1 %>% filter(BCmetric == max(BCmetric)) %>% select(pK)
pK_f1 <- as.numeric(as.character(pK_f1[[1]]))
liver_f1_doubl <- doubletFinder_v3(liver_f1, pN=0.25, pK=pK_f1, nExp=nExp_liver_f1, PCs=1:17, sct=TRUE)
DF.name_liver_f1 = colnames(liver_f1_doubl@meta.data)[grepl("DF.classification", colnames(liver_f1_doubl@meta.data))]
DimPlot(liver_f1_doubl, group.by = DF.name_liver_f1, order=c("Doublet", "Singlet"), cols=c("dodgerblue", "orange"))

liver_f2 <- subset(x=liver_umap_rename_clusters, subset=sample=="LiverF2")
#Get expected number of doublets (no. cells = 6629, estimated doublet rate ~5.1%)
nExp_liver_f2 <- round(ncol(liver_f2)*0.051)
#pK identification
sweep.res.list_f2 <- paramSweep_v3(liver_f2, PCs=1:17, sct=T)
sweep.stats_f2 <- summarizeSweep(sweep.res.list_f2, GT=F)
bcmvn_f2 <- find.pK(sweep.stats_f2)
ggplot(bcmvn_f2, aes(pK, BCmetric, group=1)) +
	geom_point() +
	geom_line()
pK_f2 <- bcmvn_f2 %>% filter(BCmetric == max(BCmetric)) %>% select(pK)
pK_f2 <- as.numeric(as.character(pK_f2[[1]]))
liver_f2_doubl <- doubletFinder_v3(liver_f2, pN=0.25, pK=pK_f2, nExp=nExp_liver_f2, PCs=1:17, sct=TRUE)
DF.name_liver_f2 = colnames(liver_f2_doubl@meta.data)[grepl("DF.classification", colnames(liver_f2_doubl@meta.data))]
DimPlot(liver_f2_doubl, group.by = DF.name_liver_f2, order=c("Doublet", "Singlet"), cols=c("dodgerblue", "orange"))

liver_f3 <- subset(x=liver_umap_rename_clusters, subset=sample=="LiverF3")
#Get expected number of doublets (no. cells = 7070, estimated doublet rate ~5.4%)
nExp_liver_f3 <- round(ncol(liver_f3)*0.054)
#pK identification
sweep.res.list_f3 <- paramSweep_v3(liver_f3, PCs=1:17, sct=T)
sweep.stats_f3 <- summarizeSweep(sweep.res.list_f3, GT=F)
bcmvn_f3 <- find.pK(sweep.stats_f3)
ggplot(bcmvn_f3, aes(pK, BCmetric, group=1)) +
	geom_point() +
	geom_line()
pK_f3 <- bcmvn_f3 %>% filter(BCmetric == max(BCmetric)) %>% select(pK)
pK_f3 <- as.numeric(as.character(pK_f3[[1]]))
liver_f3_doubl <- doubletFinder_v3(liver_f3, pN=0.25, pK=pK_f3, nExp=nExp_liver_f3, PCs=1:20, sct=TRUE)
DF.name_liver_f3 = colnames(liver_f3_doubl@meta.data)[grepl("DF.classification", colnames(liver_f3_doubl@meta.data))]
DimPlot(liver_f3_doubl, group.by = DF.name_liver_f3, order=c("Doublet", "Singlet"), cols=c("dodgerblue", "orange"))

liver_m1 <- subset(x=liver_umap_rename_clusters, subset=sample=="LiverM1")
#Get expected number of doublets (no. cells = 5244, estimated doublet rate ~4.1%)
nExp_liver_m1 <- round(ncol(liver_m1)*0.041)
#pK identification
sweep.res.list_m1 <- paramSweep_v3(liver_m1, PCs=1:17, sct=T)
sweep.stats_m1 <- summarizeSweep(sweep.res.list_m1, GT=F)
bcmvn_m1 <- find.pK(sweep.stats_m1)
ggplot(bcmvn_m1, aes(pK, BCmetric, group=1)) +
	geom_point() +
	geom_line()
pK_m1 <- bcmvn_m1 %>% filter(BCmetric == max(BCmetric)) %>% select(pK)
pK_m1 <- as.numeric(as.character(pK_m1[[1]]))
liver_m1_doubl <- doubletFinder_v3(liver_m1, pN=0.25, pK=pK_m1, nExp=nExp_liver_m1, PCs=1:20, sct=TRUE)
DF.name_liver_m1 = colnames(liver_m1_doubl@meta.data)[grepl("DF.classification", colnames(liver_m1_doubl@meta.data))]
DimPlot(liver_m1_doubl, group.by = DF.name_liver_m1, order=c("Doublet", "Singlet"), cols=c("dodgerblue", "orange"))

liver_m2 <- subset(x=liver_umap_rename_clusters, subset=sample=="LiverM2")
#Get expected number of doublets (no. cells = 6603, estimated doublet rate ~5.1%)
nExp_liver_m2 <- round(ncol(liver_m2)*0.051)
#pK identification
sweep.res.list_m2 <- paramSweep_v3(liver_m2, PCs=1:17, sct=T)
sweep.stats_m2 <- summarizeSweep(sweep.res.list_m2, GT=F)
bcmvn_m2 <- find.pK(sweep.stats_m2)
ggplot(bcmvn_m2, aes(pK, BCmetric, group=1)) +
	geom_point() +
	geom_line()
pK_m2 <- bcmvn_m2 %>% filter(BCmetric == max(BCmetric)) %>% select(pK)
pK_m2 <- as.numeric(as.character(pK_m2[[1]]))
liver_m2_doubl <- doubletFinder_v3(liver_m2, pN=0.25, pK=pK_m2, nExp=nExp_liver_m2, PCs=1:20, sct=TRUE)
DF.name_liver_m2 = colnames(liver_m2_doubl@meta.data)[grepl("DF.classification", colnames(liver_m2_doubl@meta.data))]
DimPlot(liver_m2_doubl, group.by = DF.name_liver_m2, order=c("Doublet", "Singlet"), cols=c("dodgerblue", "orange"))

liver_m3 <- subset(x=liver_umap_rename_clusters, subset=sample=="LiverM3")
#Get expected number of doublets (no. cells = 7483, estimated doublet rate ~5.7%)
nExp_liver_m3 <- round(ncol(liver_m3)*0.057)
#pK identification
sweep.res.list_m3 <- paramSweep_v3(liver_m3, PCs=1:17, sct=T)
sweep.stats_m3 <- summarizeSweep(sweep.res.list_m3, GT=F)
bcmvn_m3 <- find.pK(sweep.stats_m3)
ggplot(bcmvn_m3, aes(pK, BCmetric, group=1)) +
	geom_point() +
	geom_line()
pK_m3 <- bcmvn_m3 %>% filter(BCmetric == max(BCmetric)) %>% select(pK)
pK_m3 <- as.numeric(as.character(pK_m3[[1]]))
liver_m3_doubl <- doubletFinder_v3(liver_m3, pN=0.25, pK=pK_m3, nExp=nExp_liver_m3, PCs=1:20, sct=TRUE)
DF.name_liver_m3 = colnames(liver_m3_doubl@meta.data)[grepl("DF.classification", colnames(liver_m3_doubl@meta.data))]
DimPlot(liver_m3_doubl, group.by = DF.name_liver_m3, order=c("Doublet", "Singlet"), cols=c("dodgerblue", "orange"))


DimPlot(liver_f1_doubl, group.by = DF.name_liver_f1, order=c("Doublet", "Singlet"), cols=c("dodgerblue", "orange")) + DimPlot(liver_f2_doubl, group.by = DF.name_liver_f2, order=c("Doublet", "Singlet"), cols=c("dodgerblue", "orange")) + DimPlot(liver_f3_doubl, group.by = DF.name_liver_f3, order=c("Doublet", "Singlet"), cols=c("dodgerblue", "orange")) + DimPlot(liver_m1_doubl, group.by = DF.name_liver_m1, order=c("Doublet", "Singlet"), cols=c("dodgerblue", "orange")) + DimPlot(liver_m2_doubl, group.by = DF.name_liver_m2, order=c("Doublet", "Singlet"), cols=c("dodgerblue", "orange")) + DimPlot(liver_m3_doubl, group.by = DF.name_liver_m3, order=c("Doublet", "Singlet"), cols=c("dodgerblue", "orange"))


# Remove doublets
liver_f1_remove <- subset(x=liver_f1_doubl, subset=DF.classifications_0.25_0.02_293=="Doublet")
liver_f2_remove <- subset(x=liver_f2_doubl, subset=DF.classifications_0.25_0.03_322=="Doublet")
liver_f3_remove <- subset(x=liver_f3_doubl, subset=DF.classifications_0.25_0.005_347=="Doublet")
liver_m1_remove <- subset(x=liver_m1_doubl, subset=DF.classifications_0.25_0.005_213=="Doublet")
liver_m2_remove <- subset(x=liver_m2_doubl, subset=DF.classifications_0.25_0.02_326=="Doublet")
liver_m3_remove <- subset(x=liver_m3_doubl, subset=DF.classifications_0.25_0.03_413=="Doublet")

liver_umap_rename_clusters_dbrmv <- liver_umap_rename_clusters[,!colnames(liver_umap_rename_clusters)%in%colnames(liver_f1_remove)]
liver_umap_rename_clusters_dbrmv <- liver_umap_rename_clusters_dbrmv[,!colnames(liver_umap_rename_clusters_dbrmv)%in%colnames(liver_f2_remove)]
liver_umap_rename_clusters_dbrmv <- liver_umap_rename_clusters_dbrmv[,!colnames(liver_umap_rename_clusters_dbrmv)%in%colnames(liver_f3_remove)]
liver_umap_rename_clusters_dbrmv <- liver_umap_rename_clusters_dbrmv[,!colnames(liver_umap_rename_clusters_dbrmv)%in%colnames(liver_m1_remove)]
liver_umap_rename_clusters_dbrmv <- liver_umap_rename_clusters_dbrmv[,!colnames(liver_umap_rename_clusters_dbrmv)%in%colnames(liver_m2_remove)]
liver_umap_rename_clusters_dbrmv <- liver_umap_rename_clusters_dbrmv[,!colnames(liver_umap_rename_clusters_dbrmv)%in%colnames(liver_m3_remove)]

before <- DimPlot(liver_umap_rename_clusters, reduction="umap", label=TRUE, label.size=2.5, repel=T, label.box=T,raster=FALSE, pt.size=0.3)
after <- DimPlot(liver_umap_rename_clusters_dbrmv, reduction="umap", label=TRUE, label.size=2.5, repel=T, label.box=T,raster=FALSE, pt.size=0.3)
before + after


save(liver_umap_rename_clusters_dbrmv, file="liver_umap_rename_clusters_dbrmv.RData")

