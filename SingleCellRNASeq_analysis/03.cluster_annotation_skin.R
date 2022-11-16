
load("skin_umap.RData")

skin_umap_prepmarkers <- PrepSCTFindMarkers(object=skin_umap)

markers_skin <- FindAllMarkers(object=skin_umap_prepmarkers, assay="SCT", only.pos=TRUE)

write.table(markers_skin, file="allmarkers_skin.txt",quote=F, sep=",")

skin_umap_rename_clusters <- RenameIdents(object=skin_umap, 
										"0" = "T lymphocyte",
										"1" = "Fibroblast",
										"2" = "Stromal",
										"3" = "Fibroblast",
										"4" = "Mitochondrial rich",
										"5" = "Epidermal",
										"6" = "Keratinoycte",
										"7" = "Mesenchymal stromal",
										"8" = "Macrophage",
										"9" = "Melanocyte",
										"10" = "T lymphocyte",
										"11" = "Granulocyte",
										"12" = "B lymphocyte",
										"13" = "Epidermal (myelin-rich)",
										"14" = "Skeletal muscle",
										"15" = "Mesenchymal stromal")
										
DimPlot(skin_umap_rename_clusters, reduction="umap", label=TRUE, label.size=2.5, repel=T, label.box=T,raster=FALSE, pt.size=0.3, cols=c("#2c946f", "#f6bc27", "#e28923", "#54c6be", "#ac94f4", "tan3", "#e5243f", "#63873b", "azure4", "#fc7869", "#a6c9a1", "#643b9f", "plum1")) + NoLegend() + ggtitle("Skin") + theme(plot.title = element_text(hjust = 0.5, size=15))

DoHeatmap(skin_umap_rename_clusters, size=4, features=c("si:ch211-25d12.7","apoeb","zgc:86896","cxcl12a","soul5","CYTB","krt8","krt15","krt4","si:ch211-157c3.4","sparc","col1a2","C1QA","ENSPREG00000022933","gstm.1","pmela","tyrp1a","hsd3b7","npsn","timp2b","ENSPREG00000001025","ENSPREG00000001483","ENSPREG00000001539","mbpb","myo1cb","ENSPREG00000009304","CKM","ckma","mylz3","pvalb4"), group.colors=c("#2c946f", "#f6bc27", "#e28923", "#54c6be", "#ac94f4", "tan3", "#e5243f", "#63873b", "azure4", "#fc7869", "#a6c9a1", "#643b9f", "plum1")) + ggtitle("Skin") + theme(plot.title = element_text(hjust = 0.5))

save(skin_umap_rename_clusters, file="skin_umap_rename_clusters.RData")