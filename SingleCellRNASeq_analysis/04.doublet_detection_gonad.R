library(DoubletFinder)
libary(Seurat)

## Doublet detection ##

head(gonad_umap_3)

gonad_umap_3_renumber_clusters <- RenameIdents(object=gonad_umap_3, 
										"0" = "0",
										"1" = "1",
										"2" = "2",
										"3" = "3",
										"4" = "4",
										"5" = "5",
										"6" = "6",
										"7" = "7",
										"8" = "8",
										"9" = "9",
										"10" = "6",
										"11" = "4",
										"12" = "4",
										"13" = "13",
										"14" = "0",
										"15" = "15",
										"16" = "16",
										"17" = "15",
										"18" = "16",
										"19" = "19",
										"20" = "7",
										"21" = "21",
										"22" = "2",
										"23" = "16",
										"24" = "24",
										"25" = "6",
										"26" = "15",
										"27" = "6",
										"28" = "28",
										"29" = "6",
										"30" = "6")

gonad_f1 <- subset(x=gonad_umap_3_renumber_clusters, subset=sample=="GonadF1")
dim(gonad_f1)
[1] 19742  3259 # estimated doublet rate ~2.5%
nExp_gonad_f1 <- round(ncol(gonad_f1)*0.025)
gonad_f1_doubl <- doubletFinder_v3(gonad_f1, pN=0.25, pK=0.09, nExp=nExp_gonad_f1, PCs=1:20, sct=TRUE)
DF.name_gonad_f1 = colnames(gonad_f1_doubl@meta.data)[grepl("DF.classification", colnames(gonad_f1_doubl@meta.data))]
p1_gonad_f1 = FeaturePlot(gonad_f1, reduction="umap", features="nGene", pt.size=0.4, order=TRUE, min.cutoff='q10')
p2_gonad_f1 = FeaturePlot(gonad_f1, reduction="umap", features="nUMI", pt.size=0.4, order=TRUE, min.cutoff='q10')
p3_gonad_f1 = DimPlot(gonad_f1_doubl, group.by = DF.name_gonad_f1, order=c("Doublet", "Singlet"), cols=c("dodgerblue", "orange"))
p1_gonad_f1+p2_gonad_f1+p3_gonad_f1

gonad_f2 <- subset(x=gonad_umap_3_renumber_clusters, subset=sample=="GonadF2")
dim(gonad_f2)
[1] 19742  2100 # estimated doublet rate ~ 1.6%
nExp_gonad_f2 <- round(ncol(gonad_f2)*0.016)
gonad_f2_doubl <- doubletFinder_v3(gonad_f2, pN=0.25, pK=0.09, nExp=nExp_gonad_f2, PCs=1:20, sct=TRUE)
DF.name_gonad_f2 = colnames(gonad_f2_doubl@meta.data)[grepl("DF.classification", colnames(gonad_f2_doubl@meta.data))]
p1_gonad_f2 = FeaturePlot(gonad_f2, reduction="umap", features="nGene", pt.size=0.4, order=TRUE, min.cutoff='q10')
p2_gonad_f2 = FeaturePlot(gonad_f2, reduction="umap", features="nUMI", pt.size=0.4, order=TRUE, min.cutoff='q10')
p3_gonad_f2 = DimPlot(gonad_f2_doubl, group.by = DF.name_gonad_f2, order=c("Doublet", "Singlet"), cols=c("dodgerblue", "orange"))
p1_gonad_f2+p2_gonad_f2+p3_gonad_f2

gonad_f3 <- subset(x=gonad_umap_3_renumber_clusters, subset=sample=="GonadF3")
dim(gonad_f3)
[1] 19742  1510 # estimated doublet rate ~ 1%
nExp_gonad_f3 <- round(ncol(gonad_f3)*0.01)
gonad_f3_doubl <- doubletFinder_v3(gonad_f3, pN=0.25, pK=0.09, nExp=nExp_gonad_f3, PCs=1:20, sct=TRUE)
DF.name_gonad_f3 = colnames(gonad_f3_doubl@meta.data)[grepl("DF.classification", colnames(gonad_f3_doubl@meta.data))]
p1_gonad_f3 = FeaturePlot(gonad_f3, reduction="umap", features="nGene", pt.size=0.4, order=TRUE, min.cutoff='q10')
p2_gonad_f3 = FeaturePlot(gonad_f3, reduction="umap", features="nUMI", pt.size=0.4, order=TRUE, min.cutoff='q10')
p3_gonad_f3 = DimPlot(gonad_f3_doubl, group.by = DF.name_gonad_f3, order=c("Doublet", "Singlet"), cols=c("dodgerblue", "orange"))
p1_gonad_f3+p2_gonad_f3+p3_gonad_f3

gonad_m1 <- subset(x=gonad_umap_3_renumber_clusters, subset=sample=="GonadM4")
dim(gonad_m1)
[1] 19742  8527 # estimated doublet rate ~ 6.5%
nExp_gonad_m1 <- round(ncol(gonad_m1)*0.065)
gonad_m1_doubl <- doubletFinder_v3(gonad_m1, pN=0.25, pK=0.09, nExp=nExp_gonad_m1, PCs=1:20, sct=TRUE)
DF.name_gonad_m1 = colnames(gonad_m1_doubl@meta.data)[grepl("DF.classification", colnames(gonad_m1_doubl@meta.data))]
p1_gonad_m1 = FeaturePlot(gonad_m1, reduction="umap", features="nGene", pt.size=0.4, order=TRUE, min.cutoff='q10')
p2_gonad_m1 = FeaturePlot(gonad_m1, reduction="umap", features="nUMI", pt.size=0.4, order=TRUE, min.cutoff='q10')
p3_gonad_m1 = DimPlot(gonad_m1_doubl, group.by = DF.name_gonad_m1, order=c("Doublet", "Singlet"), cols=c("dodgerblue", "orange"))
p1_gonad_m1+p2_gonad_m1+p3_gonad_m1

gonad_m2 <- subset(x=gonad_umap_3_renumber_clusters, subset=sample=="GonadM5")
dim(gonad_m2)
[1] 17687  6396 # estimated doublet rate ~ 5%
nExp_gonad_m2 <- round(ncol(gonad_m2)*0.05)
gonad_m2_doubl <- doubletFinder_v3(gonad_m2, pN=0.25, pK=0.09, nExp=nExp_gonad_m2, PCs=1:20, sct=TRUE)
DF.name_gonad_m2 = colnames(gonad_m2_doubl@meta.data)[grepl("DF.classification", colnames(gonad_m2_doubl@meta.data))]
p1_gonad_m2 = FeaturePlot(gonad_m2, reduction="umap", features="nGene", pt.size=0.4, order=TRUE, min.cutoff='q10')
p2_gonad_m2 = FeaturePlot(gonad_m2, reduction="umap", features="nUMI", pt.size=0.4, order=TRUE, min.cutoff='q10')
p3_gonad_m2 = DimPlot(gonad_m2_doubl, group.by = DF.name_gonad_m2, order=c("Doublet", "Singlet"), cols=c("dodgerblue", "orange"))
p1_gonad_m2+p2_gonad_m2+p3_gonad_m2

gonad_m3 <- subset(x=gonad_umap_3_renumber_clusters, subset=sample=="GonadM6")
dim(gonad_m3)
[1] 19742  7186 # estimated doublet rate ~ 5.57%
nExp_gonad_m3 <- round(ncol(gonad_m3)*0.0557)
gonad_m3_doubl <- doubletFinder_v3(gonad_m3, pN=0.25, pK=0.09, nExp=nExp_gonad_m3, PCs=1:20, sct=TRUE)
DF.name_gonad_m3 = colnames(gonad_m3_doubl@meta.data)[grepl("DF.classification", colnames(gonad_m3_doubl@meta.data))]
p1_gonad_m3 = FeaturePlot(gonad_m3, reduction="umap", features="nGene", pt.size=0.4, order=TRUE, min.cutoff='q10')
p2_gonad_m3 = FeaturePlot(gonad_m3, reduction="umap", features="nUMI", pt.size=0.4, order=TRUE, min.cutoff='q10')
p3_gonad_m3 = DimPlot(gonad_m3_doubl, group.by = DF.name_gonad_m3, order=c("Doublet", "Singlet"), cols=c("dodgerblue", "orange"))
p1_gonad_m3+p2_gonad_m3+p3_gonad_m3


DimPlot(gonad_f1_doubl, group.by = DF.name_gonad_f1, order=c("Doublet", "Singlet"), cols=c("dodgerblue", "orange")) + DimPlot(gonad_f2_doubl, group.by = DF.name_gonad_f2, order=c("Doublet", "Singlet"), cols=c("dodgerblue", "orange")) + DimPlot(gonad_f3_doubl, group.by = DF.name_gonad_f3, order=c("Doublet", "Singlet"), cols=c("dodgerblue", "orange")) + DimPlot(gonad_m1_doubl, group.by = DF.name_gonad_m1, order=c("Doublet", "Singlet"), cols=c("dodgerblue", "orange")) + DimPlot(gonad_m2_doubl, group.by = DF.name_gonad_m2, order=c("Doublet", "Singlet"), cols=c("dodgerblue", "orange")) + DimPlot(gonad_m3_doubl, group.by = DF.name_gonad_m3, order=c("Doublet", "Singlet"), cols=c("dodgerblue", "orange"))


# Remove doublets
gonad_f1_remove <- subset(x=gonad_f1_doubl, subset=DF.classifications_0.25_0.09_287=="Doublet")
dim(gonad_f1_remove)
[1] 17687   287
gonad_f2_remove <- subset(x=gonad_f2_doubl, subset=DF.classifications_0.25_0.09_306=="Doublet")
dim(gonad_f2_remove)
[1] 17687   306
gonad_f3_remove <- subset(x=gonad_f3_doubl, subset=DF.classifications_0.25_0.09_317=="Doublet")
dim(gonad_f3_remove)
[1] 17687   317
gonad_m1_remove <- subset(x=gonad_m1_doubl, subset=DF.classifications_0.25_0.09_208=="Doublet")
dim(gonad_m1_remove)
[1] 17687   208
gonad_m2_remove <- subset(x=gonad_m2_doubl, subset=DF.classifications_0.25_0.09_310=="Doublet")
dim(gonad_m2_remove)
[1] 17687   310
gonad_m3_remove <- subset(x=gonad_m3_doubl, subset=DF.classifications_0.25_0.09_404=="Doublet")
dim(gonad_m3_remove)
[1] 17687   404

#cells to remove: 1832, total remaining: 35852

length(colnames(gonad_umap_3_rename_clusters))
[1] 37684
gonad_umap_3_rename_clusters_dbrmv <- gonad_umap_3_rename_clusters[,!colnames(gonad_umap_3_rename_clusters)%in%colnames(gonad_f1_remove)]
gonad_umap_3_rename_clusters_dbrmv <- gonad_umap_3_rename_clusters_dbrmv[,!colnames(gonad_umap_3_rename_clusters_dbrmv)%in%colnames(gonad_f2_remove)]
gonad_umap_3_rename_clusters_dbrmv <- gonad_umap_3_rename_clusters_dbrmv[,!colnames(gonad_umap_3_rename_clusters_dbrmv)%in%colnames(gonad_f3_remove)]
gonad_umap_3_rename_clusters_dbrmv <- gonad_umap_3_rename_clusters_dbrmv[,!colnames(gonad_umap_3_rename_clusters_dbrmv)%in%colnames(gonad_m1_remove)]
gonad_umap_3_rename_clusters_dbrmv <- gonad_umap_3_rename_clusters_dbrmv[,!colnames(gonad_umap_3_rename_clusters_dbrmv)%in%colnames(gonad_m2_remove)]
gonad_umap_3_rename_clusters_dbrmv <- gonad_umap_3_rename_clusters_dbrmv[,!colnames(gonad_umap_3_rename_clusters_dbrmv)%in%colnames(gonad_m3_remove)]
length(colnames(gonad_umap_3_rename_clusters_dbrmv))
[1] 35852

gonad_umap_3_rename_clusters_dbrmv@reductions

before <- DimPlot(gonad_umap_3_rename_clusters, reduction="umap", label=TRUE, label.size=2.5, repel=T, label.box=T,raster=FALSE, pt.size=0.3)
after <- DimPlot(gonad_umap_3_rename_clusters_dbrmv, reduction="umap", label=TRUE, label.size=2.5, repel=T, label.box=T,raster=FALSE, pt.size=0.3)
before + after


save(gonad_umap_3_rename_clusters_dbrmv, file="~/Desktop/UBC/10x/NovaSeq/cellranger_count_rawfullGTF/04.remove_doublets/gonad_filtered/gonad_umap_3_rename_clusters_dbrmv.RData")

