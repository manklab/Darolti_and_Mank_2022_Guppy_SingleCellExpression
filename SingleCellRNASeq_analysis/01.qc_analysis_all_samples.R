Load data
for (file in c("Lib2_filtered_feature_bc_matrix", "Lib5_filtered_feature_bc_matrix", "Lib6_filtered_feature_bc_matrix", "Lib7_filtered_feature_bc_matrix", "Lib8_filtered_feature_bc_matrix", "Lib9_filtered_feature_bc_matrix", "Lib10_filtered_feature_bc_matrix", "Lib11_filtered_feature_bc_matrix", "Lib12_filtered_feature_bc_matrix", "Lib14_filtered_feature_bc_matrix", "Lib16_filtered_feature_bc_matrix", "Lib18_filtered_feature_bc_matrix", "Lib19_filtered_feature_bc_matrix", "Lib21_filtered_feature_bc_matrix", "Lib23_filtered_feature_bc_matrix", "Lib25_filtered_feature_bc_matrix", "Lib26_filtered_feature_bc_matrix", "Lib27_filtered_feature_bc_matrix", "Lib32_filtered_feature_bc_matrix", "Lib33_filtered_feature_bc_matrix", "Lib34_filtered_feature_bc_matrix", "Lib35_filtered_feature_bc_matrix", "Lib36_filtered_feature_bc_matrix", "Lib37_filtered_feature_bc_matrix")){
        seurat_data <- Read10X(data.dir = paste0("~/filtered_barcodes/", file))
        seurat_obj <- CreateSeuratObject(counts = seurat_data, min.cells = 3, min.features = 100, project = file)
        assign(file, seurat_obj)
}

merged_seurat_all_samples <- merge(x = Lib2_filtered_feature_bc_matrix, 
                       y = c(Lib5_filtered_feature_bc_matrix, Lib6_filtered_feature_bc_matrix, Lib7_filtered_feature_bc_matrix, Lib8_filtered_feature_bc_matrix, Lib9_filtered_feature_bc_matrix, Lib10_filtered_feature_bc_matrix, Lib11_filtered_feature_bc_matrix, Lib12_filtered_feature_bc_matrix, Lib14_filtered_feature_bc_matrix, Lib16_filtered_feature_bc_matrix, Lib18_filtered_feature_bc_matrix, Lib19_filtered_feature_bc_matrix, Lib21_filtered_feature_bc_matrix, Lib23_filtered_feature_bc_matrix, Lib25_filtered_feature_bc_matrix, Lib26_filtered_feature_bc_matrix, Lib27_filtered_feature_bc_matrix, Lib32_filtered_feature_bc_matrix, Lib33_filtered_feature_bc_matrix, Lib34_filtered_feature_bc_matrix, Lib35_filtered_feature_bc_matrix, Lib36_filtered_feature_bc_matrix, Lib37_filtered_feature_bc_matrix),
                       add.cell.id = c("Lib2", "Lib5", "Lib6", "Lib7", "Lib8", "Lib9", "Lib10", "Lib11", "Lib12", "Lib14", "Lib16", "Lib18", "Lib19", "Lib21", "Lib23", "Lib25", "Lib26", "Lib27", "Lib32", "Lib33", "Lib34", "Lib35", "Lib36", "Lib37"))

metadata_all_samples <- merged_seurat_all_samples@meta.data

# Rename for ease of understanding
metadata_all_samples <- metadata_all_samples %>% dplyr::rename(seq_folder = orig.ident, nUMI = nCount_RNA, nGene = nFeature_RNA)

# Add cells info
metadata_all_samples$cells <- rownames(metadata_all_samples)

# Add sample info
metadata_all_samples$sample <- NA
metadata_all_samples$sample[which(str_detect(metadata_all_samples$cells, "^Lib2_"))] <- "LiverM1"  
metadata_all_samples$sample[which(str_detect(metadata_all_samples$cells, "^Lib5_"))] <- "HeartM1" 
metadata_all_samples$sample[which(str_detect(metadata_all_samples$cells, "^Lib6_"))] <- "HeartF1" 
metadata_all_samples$sample[which(str_detect(metadata_all_samples$cells, "^Lib7_"))] <- "HeartF2" 
metadata_all_samples$sample[which(str_detect(metadata_all_samples$cells, "^Lib8_"))] <- "LiverF1" 
metadata_all_samples$sample[which(str_detect(metadata_all_samples$cells, "^Lib9_"))] <- "LiverF2" 
metadata_all_samples$sample[which(str_detect(metadata_all_samples$cells, "^Lib10_"))] <- "LiverF3" 
metadata_all_samples$sample[which(str_detect(metadata_all_samples$cells, "^Lib11_"))] <- "GonadF1" 
metadata_all_samples$sample[which(str_detect(metadata_all_samples$cells, "^Lib12_"))] <- "LiverM2" 
metadata_all_samples$sample[which(str_detect(metadata_all_samples$cells, "^Lib14_"))] <- "LiverM3" 
metadata_all_samples$sample[which(str_detect(metadata_all_samples$cells, "^Lib16_"))] <- "HeartM2" 
metadata_all_samples$sample[which(str_detect(metadata_all_samples$cells, "^Lib18_"))] <- "SkinM1" 
metadata_all_samples$sample[which(str_detect(metadata_all_samples$cells, "^Lib19_"))] <- "SkinF1" 
metadata_all_samples$sample[which(str_detect(metadata_all_samples$cells, "^Lib21_"))] <- "HeartF3" 
metadata_all_samples$sample[which(str_detect(metadata_all_samples$cells, "^Lib23_"))] <- "GonadF2" 
metadata_all_samples$sample[which(str_detect(metadata_all_samples$cells, "^Lib25_"))] <- "GonadF3" 
metadata_all_samples$sample[which(str_detect(metadata_all_samples$cells, "^Lib26_"))] <- "SkinF2" 
metadata_all_samples$sample[which(str_detect(metadata_all_samples$cells, "^Lib27_"))] <- "HeartM3" 
metadata_all_samples$sample[which(str_detect(metadata_all_samples$cells, "^Lib32_"))] <- "SkinM2" 
metadata_all_samples$sample[which(str_detect(metadata_all_samples$cells, "^Lib33_"))] <- "SkinM3" 
metadata_all_samples$sample[which(str_detect(metadata_all_samples$cells, "^Lib34_"))] <- "SkinF3" 
metadata_all_samples$sample[which(str_detect(metadata_all_samples$cells, "^Lib35_"))] <- "GonadM4" 
metadata_all_samples$sample[which(str_detect(metadata_all_samples$cells, "^Lib36_"))] <- "GonadM5" 
metadata_all_samples$sample[which(str_detect(metadata_all_samples$cells, "^Lib37_"))] <- "GonadM6" 

# Add sex info
metadata_all_samples$sex <- NA
metadata_all_samples$sex[which(str_detect(metadata_all_samples$sample, "M"))] <- "Male"
metadata_all_samples$sex[which(str_detect(metadata_all_samples$sample, "F"))] <- "Female" 

# Add tissue info
metadata_all_samples$tissue <- NA
metadata_all_samples$tissue[which(str_detect(metadata_all_samples$sample, "Heart"))] <- "Heart"
metadata_all_samples$tissue[which(str_detect(metadata_all_samples$sample, "Skin"))] <- "Skin"
metadata_all_samples$tissue[which(str_detect(metadata_all_samples$sample, "Liver"))] <- "Liver"
metadata_all_samples$tissue[which(str_detect(metadata_all_samples$sample, "Gonad"))] <- "Gonad"

# Save merged seurat
merged_seurat_all_samples@meta.data <- metadata_all_samples
save(merged_seurat_all_samples, file="merged_seurat_all_samples.RData")
