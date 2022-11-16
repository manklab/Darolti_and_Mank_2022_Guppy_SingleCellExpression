
load("gonad_umap.RData")

gonad_umap_renumber_clusters <- RenameIdents(object=gonad_umap, 
										"0" = "0",
										"1" = "0",
										"2" = "0",
										"3" = "0",
										"4" = "4",
										"5" = "5",
										"6" = "6",
										"7" = "7",
										"8" = "8",
										"9" = "9",
										"10" = "10",
										"11" = "5",
										"12" = "7",
										"13" = "5",
										"14" = "14",
										"15" = "15",
										"16" = "0",
										"17" = "17",
										"18" = "0",
										"19" = "19",
										"20" = "20",
										"21" = "21",
										"22" = "22",
										"23" = "0",
										"24" = "7",
										"25" = "0",
										"26" = "0",
										"27" = "0",
										"28" = "0",
										"29" = "0",
										"30" = "0",
										"31" = "22")

gonad_umap_prepmarkers <- PrepSCTFindMarkers(object=gonad_umap_renumber_clusters)

markers_gonad <- FindAllMarkers(object=gonad_umap_prepmarkers, assay="SCT", only.pos=TRUE)

write.table(markers_gonad, file="allmarkers_gonad.txt",quote=F, sep=",")

gonad_umap_rename_clusters <- RenameIdents(object=gonad_umap, 
										"0" = "Spermatocyte",										
										"5" = "Spermatocyte",
										"9" = "Spermatocyte",
										"15" = "Spermatocyte",
										"20" = "Spermatocyte",
										"22" = "Spermatocyte",
										"7" = "Spermatozoa",
										"4" = "Granulosa",
										"6" = "Macrophage",
										"8" = "Endocrine",
										"10" = "Germ",
										"14" = "Mitochondrial-rich",
										"21" = "Embryonic Yolk Sac",
										"19" = "Endocrine",
										"17" = "Spermatocyte")

DimPlot(gonad_umap_rename_clusters, reduction="umap", label=TRUE, label.size=2.5, repel=T, label.box=T, raster=F, pt.size=0.3,cols=c("#fc7869", "#ac94f4","#e5243f", "#63873b", "#54c6be", "#449bac", "#e28923", "#f6bc27")) + NoLegend() + ggtitle("Gonad") + theme(plot.title = element_text(hjust = 0.5, size=0.5, size=15))

DoHeatmap(gonad_umap_rename_clusters, size=4, features=c("march11","tmem254","myo10","ENSPREG00000014611","ATP6","COX1","CYTB","ND2","fosab","ENSPREG00000005206","C1QA","ENSPREG00000009958","ENSPREG00000022933","col1a1b","tagln","si:ch211-156l18.7","rps26","rps27a","HBE1","ENSPREG00000007583"), group.colors=c("#fc7869", "#ac94f4","#e5243f", "#63873b", "#54c6be", "#449bac", "#e28923", "#f6bc27"))

save(gonad_umap_rename_clusters, file="gonad_umap_rename_clusters.RData")