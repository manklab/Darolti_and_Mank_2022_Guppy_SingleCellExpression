
counts_skin <- skin_umap_rename_clusters_dbrmv@assays$RNA@counts
metadata_skin <- skin_umap_rename_clusters_dbrmv@meta.data
metadata_skin$CellType <- factor(skin_umap_rename_clusters_dbrmv@active.ident)
sce_skin <- SingleCellExperiment(assays=list(counts=counts_skin), colData=metadata_skin)

# Determine the number of clusters (nk) and the cluster names (kids), and the number of samples (ns) and sample names (sids)
colData(sce_skin)$sample_id <- as.factor(colData(sce_skin)$sample)
colData(sce_skin)$sex_id <- as.factor(colData(sce_skin)$sex)
colData(sce_skin)$cluster_id <- as.factor(colData(sce_skin)$CellType)
kids_skin <- purrr::set_names(levels(sce_skin$cluster_id))
nk_skin <- length(kids_skin)
sids_skin <- purrr::set_names(levels(sce_skin$sample_id))
ns_skin <- length(sids_skin)

# Turn named vector into a numeric vector
n_cells_skin <- as.numeric(table(sce_skin$sample_id))
 
# Reorder samples (rows) of the metadata to match the order of the sample names
m_skin <- match(sids_skin, sce_skin$sample_id)

# Create sample level metadata by combining the reordered metadata with the number of cells
ei_skin <- data.frame(colData(sce_skin)[m_skin, ], n_cells_skin, row.names=NULL) %>% select(-"cluster_id")

# Remove lowly expressed genes
sce_skin <- sce_skin[rowSums(counts(sce_skin)>1)>=10,]	

# Aggregate counts per sample_id and cluster_id
groups_skin <- colData(sce_skin)[, c("cluster_id", "sample_id")]
pb_skin <- aggregate.Matrix(t(counts(sce_skin)), groupings=groups_skin, fun="sum")

# Split data by cell type and Transform matrix so that genes are row names and samples are column names
splitf_skin <- sapply(stringr::str_split(rownames(pb_skin), pattern="_", n=2), '[', 1) # Check that each cluster is present in all samples
pb_skin <- split.data.frame(pb_skin, factor(splitf_skin)) %>% lapply(function(u) set_colnames(t(u), stringr::str_extract(rownames(u), "(?<=_)[:alnum:]+")))

# Get sample names for each cluster
get_sample_ids_skin <- function(x){
	pb_skin[[x]] %>% colnames()
}
de_samples_skin <- map(1:length(kids_skin), get_sample_ids_skin) %>% unlist()

# Get cluster IDs
samples_list_skin <- map(1:length(kids_skin), get_sample_ids_skin)
get_cluster_ids_skin <- function(x){
	rep(names(pb_skin)[x], each=length(samples_list_skin[[x]]))
}
de_cluster_ids_skin <- map(1:length(kids_skin), get_cluster_ids_skin) %>% unlist()

# Create dataframe
gg_df_skin <- data.frame(cluster_id=de_cluster_ids_skin, sample_id=de_samples_skin)
gg_df_skin <- left_join(gg_df_skin, ei_skin[, c("sample_id", "sex_id", "n_cells_skin")])
metadata_skin <- gg_df_skin %>% dplyr::select(cluster_id, sample_id, sex_id, n_cells_skin)

metadata_skin$cluster_id <- as.factor(metadata_skin$cluster_id)

clusters_skin <- levels(metadata_skin$cluster_id)

# Run DESeq2 on all cell type clusters
get_dds_resultsAvsB <- function(x, A, B){
        cluster_metadata_skin <- metadata_skin[which(metadata_skin$cluster_id == clusters_skin[x]), ]
        rownames(cluster_metadata_skin) <- cluster_metadata_skin$sample_id
        counts_skin <- pb_skin[[clusters_skin[x]]]
        cluster_counts_skin <- data.frame(counts_skin[, which(colnames(counts_skin) %in% rownames(cluster_metadata_skin))])
        
        all(rownames(cluster_metadata_skin) == colnames(cluster_counts_skin))        
        
        str("Create DESeq2 object")
        dds_skin <- DESeqDataSetFromMatrix(cluster_counts_skin, 
                                      colData = cluster_metadata_skin, 
                                      design = ~ sex_id)
        raw_counts <- counts(dds_skin)
        write.csv(raw_counts, paste0("~/separate_clusters/",clusters_skin[x],"_dds_raw_counts_skin.csv"), quote = FALSE, row.names = T)
        cpm_skin <- calculateCPM(dds_skin)
        write.csv(cpm_skin, paste0("~/separate_clusters/", clusters_skin[x], "_dds_cpm_skin.csv"), quote = FALSE, row.names = T)
        
        str("Transform counts for data visualization")
        rld_skin <- rlog(dds_skin, blind=TRUE)
        
        str("Run DESeq2 differential expression analysis")
        dds_deseq_skin <- DESeq(dds_skin)
        
        str("Output results of Wald test for contrast A vs B")
        contrast_skin <- c("sex_id", levels(cluster_metadata_skin$sex_id)[A], levels(cluster_metadata_skin$sex_id)[B])
        
        str("resultsNames(dds)")
        res_skin <- results(dds_deseq_skin, 
                       contrast = contrast_skin,
                       alpha = 0.05)
        
        str("Set thresholds")
        padj_cutoff <- 0.05
        log2fc_cutoff <- 1
        
        str("Turn the results object into a tibble for use with tidyverse functions")
        res_tbl_skin <- res_skin %>%
                data.frame() %>%
                rownames_to_column(var="gene") %>%
                as_tibble()

        str("Subset the significant results")
        sig_res_skin <- dplyr::filter(res_tbl_skin, padj < padj_cutoff) %>% dplyr::arrange(padj) 
        write.csv(sig_res_skin, paste0("~/separate_clusters/", clusters_skin[x], "_deseq_sig_genes_skin.csv"), quote = FALSE, row.names = FALSE)
        
        sig_up_fc <- dplyr::filter(sig_res_skin, log2FoldChange >= log2fc_cutoff)
		sig_down_fc <- dplyr::filter(sig_res_skin, log2FoldChange <= -log2fc_cutoff)
		write.csv(sig_up_fc, paste0("~/separate_clusters/", clusters_skin[x], "_deseq_sig_genes_fcpass_malebiased_skin.csv"), quote = FALSE, row.names = FALSE)
		write.csv(sig_down_fc, paste0("~/separate_clusters/", clusters_skin[x], "_deseq_sig_genes_fcpass_femalebiased_skin.csv"), quote = FALSE, row.names = FALSE)
	
}      

map(1:length(clusters_skin), get_dds_resultsAvsB, A = 2, B = 1)






