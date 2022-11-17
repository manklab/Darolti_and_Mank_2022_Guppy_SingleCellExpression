
gonad_umap_rename_clusters_dbrmv_PB <- gonad_umap_rename_clusters_dbrmv
gonad_umap_rename_clusters_dbrmv_PB$CellType <- "All"

counts_gonad_PB <- gonad_umap_rename_clusters_dbrmv_PB@assays$RNA@counts
metadata_gonad_PB <- gonad_umap_rename_clusters_dbrmv_PB@meta.data
metadata_gonad_PB$CellType <- "All"
sce_gonad_PB <- SingleCellExperiment(assays=list(counts=counts_gonad_PB), colData=metadata_gonad_PB)

# Determine the number of clusters (nk) and the cluster names (kids), and the number of samples (ns) and sample names (sids)
colData(sce_gonad_PB)$sample_id <- as.factor(colData(sce_gonad_PB)$sample)
colData(sce_gonad_PB)$sex_id <- as.factor(colData(sce_gonad_PB)$sex)
colData(sce_gonad_PB)$cluster_id <- as.factor(colData(sce_gonad_PB)$CellType)
kids_gonad_PB <- purrr::set_names(levels(sce_gonad_PB$cluster_id))
nk_gonad_PB <- length(kids_gonad_PB)
sids_gonad_PB <- purrr::set_names(levels(sce_gonad_PB$sample_id))
ns_gonad_PB <- length(sids_gonad_PB)

# Turn named vector into a numeric vector
n_cells_gonad_PB <- as.numeric(table(sce_gonad_PB$sample_id))
 
# Reorder samples (rows) of the metadata to match the order of the sample names
m_gonad_PB <- match(sids_gonad_PB, sce_gonad_PB$sample_id)

# Create sample level metadata by combining the reordered metadata with the number of cells
ei_gonad_PB <- data.frame(colData(sce_gonad_PB)[m_gonad_PB, ], n_cells_gonad_PB, row.names=NULL) %>% select(-"cluster_id")

# Remove lowly expressed genes
sce_gonad_PB <- sce_gonad_PB[rowSums(counts(sce_gonad_PB)>1)>=10,]

# Aggregate counts per sample_id and cluster_id
groups_gonad_PB <- colData(sce_gonad_PB)[, c("cluster_id", "sample_id")]
pb_gonad_PB <- aggregate.Matrix(t(counts(sce_gonad_PB)), groupings=groups_gonad_PB, fun="sum")

# Split data by cell type and Transform matrix so that genes are row names and samples are column names
splitf_gonad_PB <- sapply(stringr::str_split(rownames(pb_gonad_PB), pattern="_", n=2), '[', 1) # Check that each cluster is present in all samples
pb_gonad_PB <- split.data.frame(pb_gonad_PB, factor(splitf_gonad_PB)) %>% lapply(function(u) set_colnames(t(u), stringr::str_extract(rownames(u), "(?<=_)[:alnum:]+")))

# Get sample names
get_sample_ids_gonad_PB <- function(x){
	pb_gonad_PB[[x]] %>% colnames()
}
de_samples_gonad_PB <- map(1:length(kids_gonad_PB), get_sample_ids_gonad_PB) %>% unlist()

# Get cluster IDs (in this case there is just one as we are merging counts across all cells)
samples_list_gonad_PB <- map(1:length(kids_gonad_PB), get_sample_ids_gonad_PB)
get_cluster_ids_gonad_PB <- function(x){
	rep(names(pb_gonad_PB)[x], each=length(samples_list_gonad_PB[[x]]))
}
de_cluster_ids_gonad_PB <- map(1:length(kids_gonad_PB), get_cluster_ids_gonad_PB) %>% unlist()

# Create dataframe
gg_df_gonad_PB <- data.frame(cluster_id=de_cluster_ids_gonad_PB, sample_id=de_samples_gonad_PB)
gg_df_gonad_PB <- left_join(gg_df_gonad_PB, ei_gonad_PB[, c("sample_id", "sex_id", "n_cells_gonad_PB")])
metadata_gonad_PB <- gg_df_gonad_PB %>% dplyr::select(cluster_id, sample_id, sex_id, n_cells_gonad_PB)

metadata_gonad_PB$cluster_id <- as.factor(metadata_gonad_PB$cluster_id)

clusters_gonad_PB <- levels(metadata_gonad_PB$cluster_id)

cluster_metadata_gonad_PB <- metadata_gonad_PB[which(metadata_gonad_PB$cluster_id == clusters_gonad_PB[1]), ]
rownames(cluster_metadata_gonad_PB) <- cluster_metadata_gonad_PB$sample_id
counts_gonad_PB <- pb_gonad_PB[[clusters_gonad_PB[1]]]
cluster_counts_gonad_PB <- data.frame(counts_gonad_PB[, which(colnames(counts_gonad_PB) %in% rownames(cluster_metadata_gonad_PB))])
all(rownames(cluster_metadata_gonad_PB) == colnames(cluster_counts_gonad_PB))    
		[1] TRUE

# Create DESeq2 object
dds_gonad_PB <- DESeqDataSetFromMatrix(cluster_counts_gonad_PB, 
				colData = cluster_metadata_gonad_PB, 
                design = ~ sex_id)
raw_counts_PB <- counts(dds_gonad_PB)
write.csv(raw_counts_PB, paste0("~/bulk/dds_raw_counts_PB_gonad.csv"), quote = FALSE, row.names = T)
cpm_PB <- calculateCPM(dds_gonad_PB)
write.csv(cpm_PB, paste0("~/bulk/dds_cpm_PB_gonad.csv"), quote = FALSE, row.names = T)

# Transform counts for data visualization
rld_gonad_PB <- rlog(dds_gonad_PB, blind=TRUE)

# Run DESeq2 differential expression analysis
dds_deseq_gonad_PB <- DESeq(dds_gonad_PB)

# Output results of Wald test for contrast A vs B
contrast_gonad_PB <- c("sex_id", levels(cluster_metadata_gonad_PB$sex_id)[2], levels(cluster_metadata_gonad_PB$sex_id)[1])
        [1] "sex_id" "Male"   "Female"
# resultsNames(dds)
res_gonad_PB <- results(dds_deseq_gonad_PB, 
                contrast = contrast_gonad_PB,
                alpha = 0.05)
                      
# Set thresholds
padj_cutoff <- 0.05
log2fc_cutoff <- 1
        
# Turn the results object into a tibble for use with tidyverse functions
res_tbl_gonad_PB <- res_gonad_PB %>%
                 data.frame() %>%
                 rownames_to_column(var="gene") %>%
                 as_tibble()
        
# Subset the significant results
sig_res_gonad_PB <- dplyr::filter(res_tbl_gonad_PB, padj < padj_cutoff) %>% dplyr::arrange(padj)    
write.csv(sig_res_gonad_PB, paste0("~/bulk/deseq_sig_genes_gonad.csv"), quote = FALSE, row.names = FALSE)

sig_up_fc_PB_gonad <- dplyr::filter(sig_res_gonad_PB, log2FoldChange >= log2fc_cutoff)
sig_down_fc_PB_gonad <- dplyr::filter(sig_res_gonad_PB, log2FoldChange <= -log2fc_cutoff)
write.csv(sig_up_fc_PB_gonad, paste0("~/bulk/deseq_sig_genes_fcpass_malebiased_gonad.csv"), quote = FALSE, row.names = FALSE)
write.csv(sig_down_fc_PB_gonad, paste0("~/bulk/deseq_sig_genes_fcpass_femalebiased_gonad.csv"), quote = FALSE, row.names = FALSE)

