
load("heart_umap.RData")

heart_umap_prepmarkers <- PrepSCTFindMarkers(object=heart_umap)

markers_heart <- FindAllMarkers(object=heart_umap_prepmarkers, assay="SCT", only.pos=TRUE)

write.table(markers_heart, file="allmarkers_heart.txt",quote=F, sep=",")

heart_umap_rename_clusters <- RenameIdents(object=heart_umap, 
										"0" = "Fibroblast",
										"1" = "Macrophage",
										"2" = "Cardiomyocyte 1",
										"3" = "Erythrocyte",
										"4" = "T lymphocyte",
										"5" = "B lymphocyte",
										"6" = "Erythrocyte",
										"7" = "Endocardial",
										"8" = "Vacular endothelial",
										"9" = "Cardiomyocyte 2",
										"10" = "Granulocyte",
										"11" = "Vascular smooth muscle")

DimPlot(heart_umap_rename_clusters, reduction="umap", label=TRUE, label.size=2.5, repel=T, label.box=T, raster=F, pt.size=0.3, cols=c("#f6bc27", "#63873b", "#643b9f", "#e5243f", "#2c946f", "#a6c9a1", "#e28923", "#54c6be", "#ac94f4", "#fc7869", "#449bac")) + NoLegend() + ggtitle("Heart") + theme(plot.title = element_text(hjust = 0.5, size=15))

DoHeatmap(heart_umap_rename_clusters, size=4, features=c("ENSPREG00000022076","krt18","fabp11a.1","C1QA","ENSPREG00000003659","ENSPREG00000022933","atp5pd","desma","myl7","tpm4a","HBE1","ENSPREG00000007583","si:ch211-103n10.5","nfkbiaa","sik1","si:ch211-25d12.7","EBF1","ENSPREG00000001483","sele","ENSPREG00000016434","smoc1","ENSPREG00000014429","ryr2b","ENSPREG00000018609","ENSPREG00000018556","hsd3b7","npsn","bmp6","pdlim7","TIMP3"), group.colors=c("#f6bc27", "#63873b", "#643b9f", "#e5243f", "#2c946f", "#a6c9a1", "#e28923", "#54c6be", "#ac94f4", "#fc7869", "#449bac")) + ggtitle("Heart") + theme(plot.title = element_text(hjust = 0.5))

save(heart_umap_rename_clusters, file="heart_umap_rename_clusters.RData")