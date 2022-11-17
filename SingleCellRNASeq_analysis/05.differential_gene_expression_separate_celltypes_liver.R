
counts_liver <- liver_umap_rename_clusters_dbrmv@assays$RNA@counts
metadata_liver <- liver_umap_rename_clusters_dbrmv@meta.data
metadata_liver$CellType <- factor(liver_umap_rename_clusters_dbrmv@active.ident)
sce_liver <- SingleCellExperiment(assays=list(counts=counts_liver), colData=metadata_liver)

# Determine the number of clusters (nk) and the cluster names (kids), and the number of samples (ns) and sample names (sids)
colData(sce_liver)$sample_id <- as.factor(colData(sce_liver)$sample)
colData(sce_liver)$sex_id <- as.factor(colData(sce_liver)$sex)
colData(sce_liver)$cluster_id <- as.factor(colData(sce_liver)$CellType)
kids_liver <- purrr::set_names(levels(sce_liver$cluster_id))
nk_liver <- length(kids_liver)
sids_liver <- purrr::set_names(levels(sce_liver$sample_id))
ns_liver <- length(sids_liver)

# Turn named vector into a numeric vector
n_cells_liver <- as.numeric(table(sce_liver$sample_id))
 
# Reorder samples (rows) of the metadata to match the order of the sample names
m_liver <- match(sids_liver, sce_liver$sample_id)

# Create sample level metadata by combining the reordered metadata with the number of cells
ei_liver <- data.frame(colData(sce_liver)[m_liver, ], n_cells_liver, row.names=NULL) %>% select(-"cluster_id")

# Remove lowly expressed genes
sce_liver <- sce_liver[rowSums(counts(sce_liver)>1)>=10,]	

# Aggregate counts per sample_id and cluster_id
groups_liver <- colData(sce_liver)[, c("cluster_id", "sample_id")]
pb_liver <- aggregate.Matrix(t(counts(sce_liver)), groupings=groups_liver, fun="sum")

# Split data by cell type and Transform matrix so that genes are row names and samples are column names
splitf_liver <- sapply(stringr::str_split(rownames(pb_liver), pattern="_", n=2), '[', 1) # Check that each cluster is present in all samples
pb_liver <- split.data.frame(pb_liver, factor(splitf_liver)) %>% lapply(function(u) set_colnames(t(u), stringr::str_extract(rownames(u), "(?<=_)[:alnum:]+")))

# Get sample names for each cluster
get_sample_ids_liver <- function(x){
	pb_liver[[x]] %>% colnames()
}
de_samples_liver <- map(1:length(kids_liver), get_sample_ids_liver) %>% unlist()

# Get cluster IDs
samples_list_liver <- map(1:length(kids_liver), get_sample_ids_liver)
get_cluster_ids_liver <- function(x){
	rep(names(pb_liver)[x], each=length(samples_list_liver[[x]]))
}
de_cluster_ids_liver <- map(1:length(kids_liver), get_cluster_ids_liver) %>% unlist()

# Create dataframe
gg_df_liver <- data.frame(cluster_id=de_cluster_ids_liver, sample_id=de_samples_liver)
gg_df_liver <- left_join(gg_df_liver, ei_liver[, c("sample_id", "sex_id", "n_cells_liver")])
metadata_liver <- gg_df_liver %>% dplyr::select(cluster_id, sample_id, sex_id, n_cells_liver)

metadata_liver$cluster_id <- as.factor(metadata_liver$cluster_id)

clusters_liver <- levels(metadata_liver$cluster_id)

# Run DESeq2 on all cell type clusters
get_dds_resultsAvsB <- function(x, A, B){
        cluster_metadata_liver <- metadata_liver[which(metadata_liver$cluster_id == clusters_liver[x]), ]
        rownames(cluster_metadata_liver) <- cluster_metadata_liver$sample_id
        counts_liver <- pb_liver[[clusters_liver[x]]]
        cluster_counts_liver <- data.frame(counts_liver[, which(colnames(counts_liver) %in% rownames(cluster_metadata_liver))])
        
        all(rownames(cluster_metadata_liver) == colnames(cluster_counts_liver))        
        
        str("Create DESeq2 object")
        dds_liver <- DESeqDataSetFromMatrix(cluster_counts_liver, 
                                      colData = cluster_metadata_liver, 
                                      design = ~ sex_id)
        raw_counts <- counts(dds_liver)
        write.csv(raw_counts, paste0("~/separate_clusters/",clusters_liver[x],"_dds_raw_counts_liver.csv"), quote = FALSE, row.names = T)
        cpm_liver <- calculateCPM(dds_liver)
        write.csv(cpm_liver, paste0("~/separate_clusters/", clusters_liver[x], "_dds_cpm_liver.csv"), quote = FALSE, row.names = T)
        
        str("Transform counts for data visualization")
        rld_liver <- rlog(dds_liver, blind=TRUE)
        
        str("Run DESeq2 differential expression analysis")
        dds_deseq_liver <- DESeq(dds_liver)
        
        str("Output results of Wald test for contrast A vs B")
        contrast_liver <- c("sex_id", levels(cluster_metadata_liver$sex_id)[A], levels(cluster_metadata_liver$sex_id)[B])
        
        str("resultsNames(dds)")
        res_liver <- results(dds_deseq_liver, 
                       contrast = contrast_liver,
                       alpha = 0.05)
        
        str("Set thresholds")
        padj_cutoff <- 0.05
        log2fc_cutoff <- 1
        
        str("Turn the results object into a tibble for use with tidyverse functions")
        res_tbl_liver <- res_liver %>%
                data.frame() %>%
                rownames_to_column(var="gene") %>%
                as_tibble()

        str("Subset the significant results")
        sig_res_liver <- dplyr::filter(res_tbl_liver, padj < padj_cutoff) %>% dplyr::arrange(padj) 
        write.csv(sig_res_liver, paste0("~/separate_clusters/", clusters_liver[x], "_deseq_sig_genes_liver.csv"), quote = FALSE, row.names = FALSE)
        
        sig_up_fc <- dplyr::filter(sig_res_liver, log2FoldChange >= log2fc_cutoff)
		sig_down_fc <- dplyr::filter(sig_res_liver, log2FoldChange <= -log2fc_cutoff)
		write.csv(sig_up_fc, paste0("~/separate_clusters/", clusters_liver[x], "_deseq_sig_genes_fcpass_malebiased_liver.csv"), quote = FALSE, row.names = FALSE)
		write.csv(sig_down_fc, paste0("~/separate_clusters/", clusters_liver[x], "_deseq_sig_genes_fcpass_femalebiased_liver.csv"), quote = FALSE, row.names = FALSE)
	
}      

map(1:length(clusters_liver), get_dds_resultsAvsB, A = 2, B = 1)






