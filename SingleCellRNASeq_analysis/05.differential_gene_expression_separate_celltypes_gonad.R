
gonad_subset <- subset(gonad_umap_rename_clusters_dbrmv, idents="Germ", invert=T)
counts_gonad <- gonad_subset@assays$RNA@counts
metadata_gonad <- gonad_subset@meta.data
metadata_gonad$CellType <- factor(gonad_subset@active.ident)
sce_gonad <- SingleCellExperiment(assays=list(counts=counts_gonad), colData=metadata_gonad)

# Determine the number of clusters (nk) and the cluster names (kids), and the number of samples (ns) and sample names (sids)
colData(sce_gonad)$sample_id <- as.factor(colData(sce_gonad)$sample)
colData(sce_gonad)$sex_id <- as.factor(colData(sce_gonad)$sex)
colData(sce_gonad)$cluster_id <- as.factor(colData(sce_gonad)$CellType)
kids_gonad <- purrr::set_names(levels(sce_gonad$cluster_id))
nk_gonad <- length(kids_gonad)
sids_gonad <- purrr::set_names(levels(sce_gonad$sample_id))
ns_gonad <- length(sids_gonad)

# Turn named vector into a numeric vector
n_cells_gonad <- as.numeric(table(sce_gonad$sample_id))
 
# Reorder samples (rows) of the metadata to match the order of the sample names
m_gonad <- match(sids_gonad, sce_gonad$sample_id)

# Create sample level metadata by combining the reordered metadata with the number of cells
ei_gonad <- data.frame(colData(sce_gonad)[m_gonad, ], n_cells_gonad, row.names=NULL) %>% select(-"cluster_id")

# Remove lowly expressed genes
sce_gonad <- sce_gonad[rowSums(counts(sce_gonad)>1)>=10,]	

# Aggregate counts per sample_id and cluster_id
groups_gonad <- colData(sce_gonad)[, c("cluster_id", "sample_id")]
pb_gonad <- aggregate.Matrix(t(counts(sce_gonad)), groupings=groups_gonad, fun="sum")

# Split data by cell type and Transform matrix so that genes are row names and samples are column names
splitf_gonad <- sapply(stringr::str_split(rownames(pb_gonad), pattern="_", n=2), '[', 1) # Check that each cluster is present in all samples
pb_gonad <- split.data.frame(pb_gonad, factor(splitf_gonad)) %>% lapply(function(u) set_colnames(t(u), stringr::str_extract(rownames(u), "(?<=_)[:alnum:]+")))

# Get sample names for each cluster
get_sample_ids_gonad <- function(x){
	pb_gonad[[x]] %>% colnames()
}
de_samples_gonad <- map(1:length(kids_gonad), get_sample_ids_gonad) %>% unlist()

# Get cluster IDs
samples_list_gonad <- map(1:length(kids_gonad), get_sample_ids_gonad)
get_cluster_ids_gonad <- function(x){
	rep(names(pb_gonad)[x], each=length(samples_list_gonad[[x]]))
}
de_cluster_ids_gonad <- map(1:length(kids_gonad), get_cluster_ids_gonad) %>% unlist()

# Create dataframe
gg_df_gonad <- data.frame(cluster_id=de_cluster_ids_gonad, sample_id=de_samples_gonad)
gg_df_gonad <- left_join(gg_df_gonad, ei_gonad[, c("sample_id", "sex_id", "n_cells_gonad")])
metadata_gonad <- gg_df_gonad %>% dplyr::select(cluster_id, sample_id, sex_id, n_cells_gonad)

metadata_gonad$cluster_id <- as.factor(metadata_gonad$cluster_id)

clusters_gonad <- levels(metadata_gonad$cluster_id)

# Run DESeq2 on all cell type clusters
get_dds_resultsAvsB <- function(x, A, B){
        cluster_metadata_gonad <- metadata_gonad[which(metadata_gonad$cluster_id == clusters_gonad[x]), ]
        rownames(cluster_metadata_gonad) <- cluster_metadata_gonad$sample_id
        counts_gonad <- pb_gonad[[clusters_gonad[x]]]
        cluster_counts_gonad <- data.frame(counts_gonad[, which(colnames(counts_gonad) %in% rownames(cluster_metadata_gonad))])
        
        all(rownames(cluster_metadata_gonad) == colnames(cluster_counts_gonad))        
        
        str("Create DESeq2 object")
        dds_gonad <- DESeqDataSetFromMatrix(cluster_counts_gonad, 
                                      colData = cluster_metadata_gonad, 
                                      design = ~ sex_id)
        raw_counts <- counts(dds_gonad)
        write.csv(raw_counts, paste0("~/separate_clusters/",clusters_gonad[x],"_dds_raw_counts_gonad.csv"), quote = FALSE, row.names = T)
        cpm_gonad <- calculateCPM(dds_gonad)
        write.csv(cpm_gonad, paste0("~/separate_clusters/", clusters_gonad[x], "_dds_cpm_gonad.csv"), quote = FALSE, row.names = T)
        
        str("Transform counts for data visualization")
        rld_gonad <- rlog(dds_gonad, blind=TRUE)
        
        str("Run DESeq2 differential expression analysis")
        dds_deseq_gonad <- DESeq(dds_gonad)
        
        str("Output results of Wald test for contrast A vs B")
        contrast_gonad <- c("sex_id", levels(cluster_metadata_gonad$sex_id)[A], levels(cluster_metadata_gonad$sex_id)[B])
        
        str("resultsNames(dds)")
        res_gonad <- results(dds_deseq_gonad, 
                       contrast = contrast_gonad,
                       alpha = 0.05)
        
        str("Set thresholds")
        padj_cutoff <- 0.05
        log2fc_cutoff <- 1
        
        str("Turn the results object into a tibble for use with tidyverse functions")
        res_tbl_gonad <- res_gonad %>%
                data.frame() %>%
                rownames_to_column(var="gene") %>%
                as_tibble()

        str("Subset the significant results")
        sig_res_gonad <- dplyr::filter(res_tbl_gonad, padj < padj_cutoff) %>% dplyr::arrange(padj) 
        write.csv(sig_res_gonad, paste0("~/separate_clusters/", clusters_gonad[x], "_deseq_sig_genes_gonad.csv"), quote = FALSE, row.names = FALSE)
        
        sig_up_fc <- dplyr::filter(sig_res_gonad, log2FoldChange >= log2fc_cutoff)
		sig_down_fc <- dplyr::filter(sig_res_gonad, log2FoldChange <= -log2fc_cutoff)
		write.csv(sig_up_fc, paste0("~/separate_clusters/", clusters_gonad[x], "_deseq_sig_genes_fcpass_malebiased_gonad.csv"), quote = FALSE, row.names = FALSE)
		write.csv(sig_down_fc, paste0("~/separate_clusters/", clusters_gonad[x], "_deseq_sig_genes_fcpass_femalebiased_gonad.csv"), quote = FALSE, row.names = FALSE)
	
}      

map(1:length(clusters_gonad), get_dds_resultsAvsB, A = 2, B = 1)






