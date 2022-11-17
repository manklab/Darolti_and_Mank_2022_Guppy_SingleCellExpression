
counts_heart <- heart_umap_rename_clusters_dbrmv@assays$RNA@counts
metadata_heart <- heart_umap_rename_clusters_dbrmv@meta.data
metadata_heart$CellType <- factor(heart_umap_rename_clusters_dbrmv@active.ident)
sce_heart <- SingleCellExperiment(assays=list(counts=counts_heart), colData=metadata_heart)

# Determine the number of clusters (nk) and the cluster names (kids), and the number of samples (ns) and sample names (sids)
colData(sce_heart)$sample_id <- as.factor(colData(sce_heart)$sample)
colData(sce_heart)$sex_id <- as.factor(colData(sce_heart)$sex)
colData(sce_heart)$cluster_id <- as.factor(colData(sce_heart)$CellType)
kids_heart <- purrr::set_names(levels(sce_heart$cluster_id))
nk_heart <- length(kids_heart)
sids_heart <- purrr::set_names(levels(sce_heart$sample_id))
ns_heart <- length(sids_heart)

# Turn named vector into a numeric vector
n_cells_heart <- as.numeric(table(sce_heart$sample_id))
 
# Reorder samples (rows) of the metadata to match the order of the sample names
m_heart <- match(sids_heart, sce_heart$sample_id)

# Create sample level metadata by combining the reordered metadata with the number of cells
ei_heart <- data.frame(colData(sce_heart)[m_heart, ], n_cells_heart, row.names=NULL) %>% select(-"cluster_id")

# Remove lowly expressed genes
sce_heart <- sce_heart[rowSums(counts(sce_heart)>1)>=10,]	

# Aggregate counts per sample_id and cluster_id
groups_heart <- colData(sce_heart)[, c("cluster_id", "sample_id")]
pb_heart <- aggregate.Matrix(t(counts(sce_heart)), groupings=groups_heart, fun="sum")

# Split data by cell type and Transform matrix so that genes are row names and samples are column names
splitf_heart <- sapply(stringr::str_split(rownames(pb_heart), pattern="_", n=2), '[', 1) # Check that each cluster is present in all samples
pb_heart <- split.data.frame(pb_heart, factor(splitf_heart)) %>% lapply(function(u) set_colnames(t(u), stringr::str_extract(rownames(u), "(?<=_)[:alnum:]+")))

# Get sample names for each cluster
get_sample_ids_heart <- function(x){
	pb_heart[[x]] %>% colnames()
}
de_samples_heart <- map(1:length(kids_heart), get_sample_ids_heart) %>% unlist()

# Get cluster IDs
samples_list_heart <- map(1:length(kids_heart), get_sample_ids_heart)
get_cluster_ids_heart <- function(x){
	rep(names(pb_heart)[x], each=length(samples_list_heart[[x]]))
}
de_cluster_ids_heart <- map(1:length(kids_heart), get_cluster_ids_heart) %>% unlist()

# Create dataframe
gg_df_heart <- data.frame(cluster_id=de_cluster_ids_heart, sample_id=de_samples_heart)
gg_df_heart <- left_join(gg_df_heart, ei_heart[, c("sample_id", "sex_id", "n_cells_heart")])
metadata_heart <- gg_df_heart %>% dplyr::select(cluster_id, sample_id, sex_id, n_cells_heart)

metadata_heart$cluster_id <- as.factor(metadata_heart$cluster_id)

clusters_heart <- levels(metadata_heart$cluster_id)

# Run DESeq2 on all cell type clusters
get_dds_resultsAvsB <- function(x, A, B){
        cluster_metadata_heart <- metadata_heart[which(metadata_heart$cluster_id == clusters_heart[x]), ]
        rownames(cluster_metadata_heart) <- cluster_metadata_heart$sample_id
        counts_heart <- pb_heart[[clusters_heart[x]]]
        cluster_counts_heart <- data.frame(counts_heart[, which(colnames(counts_heart) %in% rownames(cluster_metadata_heart))])
        
        all(rownames(cluster_metadata_heart) == colnames(cluster_counts_heart))        
        
        str("Create DESeq2 object")
        dds_heart <- DESeqDataSetFromMatrix(cluster_counts_heart, 
                                      colData = cluster_metadata_heart, 
                                      design = ~ sex_id)
        raw_counts <- counts(dds_heart)
        write.csv(raw_counts, paste0("~/separate_clusters/",clusters_heart[x],"_dds_raw_counts_heart.csv"), quote = FALSE, row.names = T)
        cpm_heart <- calculateCPM(dds_heart)
        write.csv(cpm_heart, paste0("~/separate_clusters/", clusters_heart[x], "_dds_cpm_heart.csv"), quote = FALSE, row.names = T)
        
        str("Transform counts for data visualization")
        rld_heart <- rlog(dds_heart, blind=TRUE)
        
        str("Run DESeq2 differential expression analysis")
        dds_deseq_heart <- DESeq(dds_heart)
        
        str("Output results of Wald test for contrast A vs B")
        contrast_heart <- c("sex_id", levels(cluster_metadata_heart$sex_id)[A], levels(cluster_metadata_heart$sex_id)[B])
        
        str("resultsNames(dds)")
        res_heart <- results(dds_deseq_heart, 
                       contrast = contrast_heart,
                       alpha = 0.05)
        
        str("Set thresholds")
        padj_cutoff <- 0.05
        log2fc_cutoff <- 1
        
        str("Turn the results object into a tibble for use with tidyverse functions")
        res_tbl_heart <- res_heart %>%
                data.frame() %>%
                rownames_to_column(var="gene") %>%
                as_tibble()

        str("Subset the significant results")
        sig_res_heart <- dplyr::filter(res_tbl_heart, padj < padj_cutoff) %>% dplyr::arrange(padj) 
        write.csv(sig_res_heart, paste0("~/separate_clusters/", clusters_heart[x], "_deseq_sig_genes_heart.csv"), quote = FALSE, row.names = FALSE)
        
        sig_up_fc <- dplyr::filter(sig_res_heart, log2FoldChange >= log2fc_cutoff)
		sig_down_fc <- dplyr::filter(sig_res_heart, log2FoldChange <= -log2fc_cutoff)
		write.csv(sig_up_fc, paste0("~/separate_clusters/", clusters_heart[x], "_deseq_sig_genes_fcpass_malebiased_heart.csv"), quote = FALSE, row.names = FALSE)
		write.csv(sig_down_fc, paste0("~/separate_clusters/", clusters_heart[x], "_deseq_sig_genes_fcpass_femalebiased_heart.csv"), quote = FALSE, row.names = FALSE)
	
}      

map(1:length(clusters_heart), get_dds_resultsAvsB, A = 2, B = 1)






