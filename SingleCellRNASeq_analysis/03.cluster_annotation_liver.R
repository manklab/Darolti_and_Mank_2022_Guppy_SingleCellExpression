
load("liver_umap.RData")

liver_umap_prepmarkers <- PrepSCTFindMarkers(object=liver_umap)

markers_liver <- FindAllMarkers(object=liver_umap_prepmarkers, assay="SCT", only.pos=TRUE)

write.table(markers_liver, file="allmarkers_liver.txt",quote=F, sep=",")

liver_umap_rename_clusters <- RenameIdents(object=liver_umap, 
										"0" = "Hepatocyte 1",
										"1" = "Hepatocyte 1",
										"2" = "T lymphocyte",
										"3" = "Hepatocyte 2",
										"4" = "Hepatocyte 2",
										"5" = "Hepatocyte 2",
										"6" = "Macrophage",
										"7" = "B lymphocyte",
										"8" = "Biliary epithelial",
										"9" = "Erythrocyte",
										"10" = "Granulocyte",
										"11" = "Hepatocyte 1",
										"12" = "Erythrocyte",
										"13" = "Endothelial",
										"14" = "Endothelial")

DimPlot(liver_umap_rename_clusters, reduction="umap", label=TRUE, label.size=2.5, repel=T, label.box=T, raster=F, pt.size=0.3, cols=c("#ac94f4", "#2c946f", "#643b9f", "#63873b", "#a6c9a1", "#f6bc27", "#e5243f", "#fc7869", "#e28923") ) + NoLegend() + ggtitle("Liver") + theme(plot.title = element_text(hjust = 0.5, size=15))

DoHeatmap(liver_umap_rename_clusters, size=4, features=c("ENSPREG00000022089","sik1","si:ch211-25d12.7","fabp10a","tfa","C1QA","HMOX1","ENSPREG00000001483","ENSPREG00000001025","anxa4","epcam","lgals2b","HBE1","ENSPREG00000007451","hsd3b7","npsn","fabp1b.1"), group.colors=c("#ac94f4", "#2c946f", "#643b9f", "#63873b", "#a6c9a1", "#f6bc27", "#e5243f", "#fc7869", "#e28923")) + ggtitle("Liver") + theme(plot.title = element_text(hjust = 0.5))

save(liver_umap_rename_clusters, file="liver_umap_rename_clusters.RData")