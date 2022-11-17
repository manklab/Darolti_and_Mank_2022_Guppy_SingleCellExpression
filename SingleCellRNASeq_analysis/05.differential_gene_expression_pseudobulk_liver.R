
liver_umap_rename_clusters_dbrmv_PB <- liver_umap_rename_clusters_dbrmv
liver_umap_rename_clusters_dbrmv_PB$CellType <- "All"

counts_liver_PB <- liver_umap_rename_clusters_dbrmv_PB@assays$RNA@counts
metadata_liver_PB <- liver_umap_rename_clusters_dbrmv_PB@meta.data
metadata_liver_PB$CellType <- "All"
sce_liver_PB <- SingleCellExperiment(assays=list(counts=counts_liver_PB), colData=metadata_liver_PB)

# Determine the number of clusters (nk) and the cluster names (kids), and the number of samples (ns) and sample names (sids)
colData(sce_liver_PB)$sample_id <- as.factor(colData(sce_liver_PB)$sample)
colData(sce_liver_PB)$sex_id <- as.factor(colData(sce_liver_PB)$sex)
colData(sce_liver_PB)$cluster_id <- as.factor(colData(sce_liver_PB)$CellType)
kids_liver_PB <- purrr::set_names(levels(sce_liver_PB$cluster_id))
nk_liver_PB <- length(kids_liver_PB)
sids_liver_PB <- purrr::set_names(levels(sce_liver_PB$sample_id))
ns_liver_PB <- length(sids_liver_PB)

# Turn named vector into a numeric vector
n_cells_liver_PB <- as.numeric(table(sce_liver_PB$sample_id))
 
# Reorder samples (rows) of the metadata to match the order of the sample names
m_liver_PB <- match(sids_liver_PB, sce_liver_PB$sample_id)

# Create sample level metadata by combining the reordered metadata with the number of cells
ei_liver_PB <- data.frame(colData(sce_liver_PB)[m_liver_PB, ], n_cells_liver_PB, row.names=NULL) %>% select(-"cluster_id")

# Remove lowly expressed genes
sce_liver_PB <- sce_liver_PB[rowSums(counts(sce_liver_PB)>1)>=10,]

# Aggregate counts per sample_id and cluster_id
groups_liver_PB <- colData(sce_liver_PB)[, c("cluster_id", "sample_id")]
pb_liver_PB <- aggregate.Matrix(t(counts(sce_liver_PB)), groupings=groups_liver_PB, fun="sum")

# Split data by cell type and Transform matrix so that genes are row names and samples are column names
splitf_liver_PB <- sapply(stringr::str_split(rownames(pb_liver_PB), pattern="_", n=2), '[', 1) # Check that each cluster is present in all samples
pb_liver_PB <- split.data.frame(pb_liver_PB, factor(splitf_liver_PB)) %>% lapply(function(u) set_colnames(t(u), stringr::str_extract(rownames(u), "(?<=_)[:alnum:]+")))

# Get sample names
get_sample_ids_liver_PB <- function(x){
	pb_liver_PB[[x]] %>% colnames()
}
de_samples_liver_PB <- map(1:length(kids_liver_PB), get_sample_ids_liver_PB) %>% unlist()

# Get cluster IDs (in this case there is just one as we are merging counts across all cells)
samples_list_liver_PB <- map(1:length(kids_liver_PB), get_sample_ids_liver_PB)
get_cluster_ids_liver_PB <- function(x){
	rep(names(pb_liver_PB)[x], each=length(samples_list_liver_PB[[x]]))
}
de_cluster_ids_liver_PB <- map(1:length(kids_liver_PB), get_cluster_ids_liver_PB) %>% unlist()

# Create dataframe
gg_df_liver_PB <- data.frame(cluster_id=de_cluster_ids_liver_PB, sample_id=de_samples_liver_PB)
gg_df_liver_PB <- left_join(gg_df_liver_PB, ei_liver_PB[, c("sample_id", "sex_id", "n_cells_liver_PB")])
metadata_liver_PB <- gg_df_liver_PB %>% dplyr::select(cluster_id, sample_id, sex_id, n_cells_liver_PB)

metadata_liver_PB$cluster_id <- as.factor(metadata_liver_PB$cluster_id)

clusters_liver_PB <- levels(metadata_liver_PB$cluster_id)

cluster_metadata_liver_PB <- metadata_liver_PB[which(metadata_liver_PB$cluster_id == clusters_liver_PB[1]), ]
rownames(cluster_metadata_liver_PB) <- cluster_metadata_liver_PB$sample_id
counts_liver_PB <- pb_liver_PB[[clusters_liver_PB[1]]]
cluster_counts_liver_PB <- data.frame(counts_liver_PB[, which(colnames(counts_liver_PB) %in% rownames(cluster_metadata_liver_PB))])
all(rownames(cluster_metadata_liver_PB) == colnames(cluster_counts_liver_PB))    
		[1] TRUE

# Create DESeq2 object
dds_liver_PB <- DESeqDataSetFromMatrix(cluster_counts_liver_PB, 
				colData = cluster_metadata_liver_PB, 
                design = ~ sex_id)
raw_counts_PB <- counts(dds_liver_PB)
write.csv(raw_counts_PB, paste0("~/bulk/dds_raw_counts_PB_liver.csv"), quote = FALSE, row.names = T)
cpm_PB <- calculateCPM(dds_liver_PB)
write.csv(cpm_PB, paste0("~/bulk/dds_cpm_PB_liver.csv"), quote = FALSE, row.names = T)

# Transform counts for data visualization
rld_liver_PB <- rlog(dds_liver_PB, blind=TRUE)

# Run DESeq2 differential expression analysis
dds_deseq_liver_PB <- DESeq(dds_liver_PB)

# Output results of Wald test for contrast A vs B
contrast_liver_PB <- c("sex_id", levels(cluster_metadata_liver_PB$sex_id)[2], levels(cluster_metadata_liver_PB$sex_id)[1])
        [1] "sex_id" "Male"   "Female"
# resultsNames(dds)
res_liver_PB <- results(dds_deseq_liver_PB, 
                contrast = contrast_liver_PB,
                alpha = 0.05)
                      
# Set thresholds
padj_cutoff <- 0.05
log2fc_cutoff <- 1
        
# Turn the results object into a tibble for use with tidyverse functions
res_tbl_liver_PB <- res_liver_PB %>%
                 data.frame() %>%
                 rownames_to_column(var="gene") %>%
                 as_tibble()
        
# Subset the significant results
sig_res_liver_PB <- dplyr::filter(res_tbl_liver_PB, padj < padj_cutoff) %>% dplyr::arrange(padj)    
write.csv(sig_res_liver_PB, paste0("~/bulk/deseq_sig_genes_liver.csv"), quote = FALSE, row.names = FALSE)

sig_up_fc_PB_liver <- dplyr::filter(sig_res_liver_PB, log2FoldChange >= log2fc_cutoff)
sig_down_fc_PB_liver <- dplyr::filter(sig_res_liver_PB, log2FoldChange <= -log2fc_cutoff)
write.csv(sig_up_fc_PB_liver, paste0("~/bulk/deseq_sig_genes_fcpass_malebiased_liver.csv"), quote = FALSE, row.names = FALSE)
write.csv(sig_down_fc_PB_liver, paste0("~/bulk/deseq_sig_genes_fcpass_femalebiased_liver.csv"), quote = FALSE, row.names = FALSE)

