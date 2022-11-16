library(SingleCellExperiment)
library(Seurat)
library(Matrix.utils)
library(Matrix)
library(dplyr)
library(magrittr)
library(purrr)
library(stringr)
library(DESeq2)
library(pheatmap)
library(apeglm)
library(tibble)
library(tidyr)
library(ggplot2)
library(RColorBrewer)
library(scuttle)
library(ashr)

load("~/Desktop/UBC/10x/NovaSeq/cellranger_count_rawfullGTF/04.remove_doublets/liver_filtered/liver_umap_3_rename_clusters_dbrmv.RData")

liver_umap_3_rename_clusters_dbrmv_PB <- liver_umap_3_rename_clusters_dbrmv
liver_umap_3_rename_clusters_dbrmv_PB$CellType <- "All"

head(liver_umap_3_rename_clusters_dbrmv_PB)

counts_liver_PB <- liver_umap_3_rename_clusters_dbrmv_PB@assays$RNA@counts
metadata_liver_PB <- liver_umap_3_rename_clusters_dbrmv_PB@meta.data
metadata_liver_PB$CellType <- "All"

sce_liver_PB <- SingleCellExperiment(assays=list(counts=counts_liver_PB), colData=metadata_liver_PB)
sce_liver_PB
class: SingleCellExperiment 
dim: 21975 35770 
metadata(0):
assays(1): counts
rownames(21975): ENSPREG00000006268 arl13b ... guca1d ENSPREG00000006408
rowData names(0):
colnames(35770): F1_Lib8_AAACCCAAGAGTCTTC-1 F1_Lib8_AAACCCAAGGATTTCC-1 ...
  M3_Lib14_TTTGTTGTCGAGTCCG-1 M3_Lib14_TTTGTTGTCTTACACT-1
colData names(27): seq_folder nUMI ... seurat_clusters CellType
reducedDimNames(0):
altExpNames(0):

assays(sce_liver_PB)
		List of length 1
		names(1): counts
dim(counts(sce_liver_PB))
[1] 21975 35770
counts(sce_liver_PB)[1:2, 1:2]
		2 x 2 sparse Matrix of class "dgCMatrix"
		                   F1_Lib8_AAACCCAAGAGTCTTC-1 F1_Lib8_AAACCCAAGGATTTCC-1
		ENSPREG00000006268                          .                          .
		arl13b                                      .                          .
dim(colData(sce_liver_PB))
[1] 35770    27

# Determine the number of clusters and the cluster names
colData(sce_liver_PB)$sample_id <- as.factor(colData(sce_liver_PB)$sample)
colData(sce_liver_PB)$sex_id <- as.factor(colData(sce_liver_PB)$sex)
colData(sce_liver_PB)$cluster_id <- as.factor(colData(sce_liver_PB)$CellType)
kids_liver_PB <- purrr::set_names(levels(sce_liver_PB$cluster_id))
kids_liver_PB
		  All 
		"All" 
nk_liver_PB <- length(kids_liver_PB)
nk_liver_PB
		[1] 1

sids_liver_PB <- purrr::set_names(levels(sce_liver_PB$sample_id))
sids_liver_PB
		  LiverF1   LiverF2   LiverF3   LiverM1   LiverM2   LiverM3 
		"LiverF1" "LiverF2" "LiverF3" "LiverM1" "LiverM2" "LiverM3" 
ns_liver_PB <- length(sids_liver_PB)
ns_liver_PB
		[1] 6

# Determine no. cells per sample
table(sce_liver_PB$sample_id)
		LiverF1 LiverF2 LiverF3 LiverM1 LiverM2 LiverM3 
   5814    5991    6075    4985    6070    6835

# Turn named vector into a numeric vector
n_cells_liver_PB <- as.numeric(table(sce_liver_PB$sample_id))
 
# Reorder samples (rows) of the metadata to match the order of the sample names
m_liver_PB <- match(sids_liver_PB, sce_liver_PB$sample_id)
m_liver_PB
[1]     1  5815 11806 17881 22866 28936

# Create sample level metadata by combining the reordered metadata with the number of cells
ei_liver_PB <- data.frame(colData(sce_liver_PB)[m_liver_PB, ], n_cells_liver_PB, row.names=NULL) %>% select(-"cluster_id")

# Remove lowly expressed genes
dim(sce_liver_PB)
[1] 21975 35770
sce_liver_PB <- sce_liver_PB[rowSums(counts(sce_liver_PB)>1)>=10,]
dim(sce_liver_PB)
[1] 11478 35770

# Aggregate counts per sample_id and cluster_id
groups_liver_PB <- colData(sce_liver_PB)[, c("cluster_id", "sample_id")]
pb_liver_PB <- aggregate.Matrix(t(counts(sce_liver_PB)), groupings=groups_liver_PB, fun="sum")
class(pb_liver_PB)
		[1] "dgCMatrix"
		attr(,"package")
		[1] "Matrix"
dim(pb_liver_PB)
[1]     6 11478

pb_liver_PB[1:6, 1:6]
6 x 6 sparse Matrix of class "dgCMatrix"
            ENSPREG00000006268 gart n6amt1 timmdc1 VANGL1 ENSPREG00000006674
All_LiverF1                  5   11     13      77     29               6254
All_LiverF2                 26   20     11      42     33               5237
All_LiverF3                 23   19     16      56     71               4645
All_LiverM1                 48   46     30      82     24               1680
All_LiverM2                 22   42     14      61     68              20003
All_LiverM3                 60   40     55     133    204               2395


## Split data by cell type and Transform matrix so that genes are row names and samples are column names
splitf_liver_PB <- sapply(stringr::str_split(rownames(pb_liver_PB), pattern="_", n=2), '[', 1) # Check that each cluster is present in all samples
pb_liver_PB <- split.data.frame(pb_liver_PB, factor(splitf_liver_PB)) %>% lapply(function(u) set_colnames(t(u), stringr::str_extract(rownames(u), "(?<=_)[:alnum:]+")))

class(pb_liver_PB)
		[1] "list"
str(pb_liver_PB)

options(width=100)
table(sce_liver_PB$cluster_id, sce_liver_PB$sample_id)
		      LiverF1 LiverF2 LiverF3 LiverM1 LiverM2 LiverM3
  All    5814    5991    6075    4985    6070    6835
  
get_sample_ids_liver_PB <- function(x){
	pb_liver_PB[[x]] %>% colnames()
}
de_samples_liver_PB <- map(1:length(kids_liver_PB), get_sample_ids_liver_PB) %>% unlist()

samples_list_liver_PB <- map(1:length(kids_liver_PB), get_sample_ids_liver_PB)

get_cluster_ids_liver_PB <- function(x){
	rep(names(pb_liver_PB)[x], each=length(samples_list_liver_PB[[x]]))
}

de_cluster_ids_liver_PB <- map(1:length(kids_liver_PB), get_cluster_ids_liver_PB) %>% unlist()

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
write.csv(raw_counts_PB, paste0("~/Desktop/UBC/10x/NovaSeq/cellranger_count_rawfullGTF/05.differential_gene_expression/liver_filtered/DESeq2_filterlowexpr/pseudobulk/dds_raw_counts_PB.csv"), quote = FALSE, row.names = T)

#-----------
fpm_PB <- fpm(dds_liver_PB)
head(fpm_PB)
                       LiverF1     LiverF2    LiverF3    LiverM1     LiverM2    LiverM3
ENSPREG00000006268   0.3754481   2.2361327   1.848998   4.623911    1.814487   3.324094
gart                 0.8259858   1.7201021   1.527433   4.431248    3.464020   2.216063
n6amt1               0.9761651   0.9460561   1.286260   2.889944    1.154673   3.047086
timmdc1              5.7819008   3.6122144   4.501908   7.899181    5.031077   7.368408
VANGL1               2.1775990   2.8381684   5.707777   2.311955    5.608414  11.301919
ENSPREG00000006674 469.6104879 450.4087301 373.417216 161.836873 1649.780897 132.686750

write.csv(fpm_PB, paste0("~/Desktop/UBC/10x/NovaSeq/cellranger_count_rawfullGTF/05.differential_gene_expression/liver_filtered/DESeq2_filterlowexpr/pseudobulk/dds_fpm_PB.csv"), quote = FALSE, row.names = T)

cpm_PB <- calculateCPM(dds_liver_PB)
head(cpm_PB)
              LiverF1    LiverF2    LiverF3    LiverM1    LiverM2    LiverM3
ENSPREG00000006268   0.3494988   2.1066220   1.667732   6.040856    1.4062216   4.310594
gart                 0.7688973   1.6204784   1.377692   5.789153    2.6846049   2.873730
n6amt1               0.9086968   0.8912631   1.160162   3.775535    0.8948683   3.951378
timmdc1              5.3822811   3.4030047   4.060566  10.319795    3.8990691   9.555151
VANGL1               2.0270929   2.6737894   5.148218   3.020428    4.3465032  14.656021
ENSPREG00000006674 437.1530613 424.3222795 336.809444 211.429954 1278.5750628 172.064561

write.csv(cpm_PB, paste0("~/Desktop/UBC/10x/NovaSeq/cellranger_count_rawfullGTF/05.differential_gene_expression/liver_filtered/DESeq2_filterlowexpr/pseudobulk/dds_cpm_PB.csv"), quote = FALSE, row.names = T)
#--------------


# Transform counts for data visualization
rld_liver_PB <- rlog(dds_liver_PB, blind=TRUE)
colData(rld_liver_PB)
DataFrame with 6 rows and 5 columns
        cluster_id   sample_id   sex_id n_cells_liver_PB sizeFactor
          <factor> <character> <factor>        <numeric>  <numeric>
LiverF1        All     LiverF1   Female             5814   1.047885
LiverF2        All     LiverF2   Female             5991   0.914891
LiverF3        All     LiverF3   Female             6075   0.978780
LiverM1        All     LiverM1   Male               4985   0.816818
LiverM2        All     LiverM2   Male               6070   0.954031
LiverM3        All     LiverM3   Male               6835   1.420272
        
# Plot PCA
DESeq2::plotPCA(rld_liver_PB, ntop=100,intgroup = "sex_id")
ggsave(paste0("~/Desktop/UBC/10x/NovaSeq/cellranger_count_rawfullGTF/05.differential_gene_expression/liver_filtered/DESeq2_filterlowexpr/pseudobulk/plots/pseudobulk_PCAplot.png"))

        
        DESeq2::plotPCA(rld_liver_PB, ntop=100, intgroup = "n_cells_liver_PB")
        ggsave(paste0("~/Desktop/UBC/10x/NovaSeq/cellranger_count_rawfullGTF/05.differential_gene_expression/liver_filtered/DESeq2_filterlowexpr/pseudobulk/plots/pseudobulk_PCAplot_cellcount.png"))

# Extract the rlog matrix from the object and compute pairwise correlation values
rld_mat_liver_PB <- assay(rld_liver_PB)
head(rld_mat_liver_PB)
                     LiverF1   LiverF2   LiverF3   LiverM1   LiverM2   LiverM3
ENSPREG00000006268  4.223859  4.794013  4.701428  5.221848  4.692856  5.018354
gart                4.419940  4.687395  4.635803  5.206945  5.054253  4.806787
n6amt1              4.175666  4.168260  4.270932  4.643388  4.232427  4.681310
timmdc1             6.170414  5.858734  5.998043  6.401727  6.072986  6.351448
VANGL1              5.385566  5.532731  5.989602  5.423270  5.976644  6.537391
ENSPREG00000006674 12.458272 12.416332 12.230139 11.446658 13.786159 11.272982

write.csv(rld_mat_liver_PB, paste0("~/Desktop/UBC/10x/NovaSeq/cellranger_count_rawfullGTF/05.differential_gene_expression/liver_filtered/DESeq2_filterlowexpr/pseudobulk/rlogtransformed_counts_PB.csv"), quote = FALSE, row.names = T)
rld_cor_liver_PB <- cor(rld_mat_liver_PB)
write.csv(rld_cor_liver_PB, paste0("~/Desktop/UBC/10x/NovaSeq/cellranger_count_rawfullGTF/05.differential_gene_expression/liver_filtered/DESeq2_filterlowexpr/pseudobulk/deseq2_pairwise_correlations_PB.csv"), quote = FALSE, row.names = T)
        
# Plot heatmap
png(paste0("~/Desktop/UBC/10x/NovaSeq/cellranger_count_rawfullGTF/05.differential_gene_expression/liver_filtered/DESeq2_filterlowexpr/pseudobulk/plots/pseudobulk_heatmap.png"), height=6, width=7.5, units="in", res=300)
pheatmap(rld_cor_liver_PB, annotation = cluster_metadata_liver_PB[, c("sex_id"), drop=F])
dev.off()
        
# Run DESeq2 differential expression analysis
dds_deseq_liver_PB <- DESeq(dds_liver_PB)
#estimating size factors
#estimating dispersions
#gene-wise dispersion estimates
#mean-dispersion relationship
#final dispersion estimates
#fitting model and testing
        
# Plot dispersion estimates
png(paste0("~/Desktop/UBC/10x/NovaSeq/cellranger_count_rawfullGTF/05.differential_gene_expression/liver_filtered/DESeq2_filterlowexpr/pseudobulk/plots/pseudobulk_dispersion_plot.png"))
plotDispEsts(dds_deseq_liver_PB)
dev.off()

# Output results of Wald test for contrast for A vs B
contrast_liver_PB <- c("sex_id", levels(cluster_metadata_liver_PB$sex_id)[2], levels(cluster_metadata_liver_PB$sex_id)[1])
        [1] "sex_id" "Male"   "Female"
# resultsNames(dds)
res_liver_PB <- results(dds_deseq_liver_PB, 
                contrast = contrast_liver_PB,
                alpha = 0.05)

res_liver_PB_ashr <- lfcShrink(dds_deseq_liver_PB, 
                contrast=contrast_liver_PB,
                res=res_liver_PB, type="ashr")
                     
 # Try other adjsuting methods
 res_liver_PB_bon <- results(dds_deseq_liver_PB, 
                contrast = contrast_liver_PB,
                alpha = 0.05,
                pAdjustMethod="bonferroni")                   
                     
                      
# Set thresholds
padj_cutoff <- 0.05
log2fc_cutoff <- 1
        
# Turn the results object into a tibble for use with tidyverse functions

res_tbl_liver_PB <- res_liver_PB %>%
                 data.frame() %>%
                 rownames_to_column(var="gene") %>%
                 as_tibble()
write.csv(res_tbl_liver_PB, paste0("~/Desktop/UBC/10x/NovaSeq/cellranger_count_rawfullGTF/05.differential_gene_expression/liver_filtered/DESeq2_filterlowexpr/pseudobulk/deseq_beforeshrink_all_genes.csv"), quote = FALSE, row.names = FALSE)


res_tbl_liver_PB_ashr <- res_liver_PB_ashr %>%
                 data.frame() %>%
                 rownames_to_column(var="gene") %>%
                 as_tibble()
write.csv(res_tbl_liver_PB_ashr, paste0("~/Desktop/UBC/10x/NovaSeq/cellranger_count_rawfullGTF/05.differential_gene_expression/liver_filtered/DESeq2_filterlowexpr/pseudobulk/deseq_shrinkashr_all_genes.csv"), quote = FALSE, row.names = FALSE)
        
# Subset the significant results
sig_res_liver_PB <- dplyr::filter(res_tbl_liver_PB, padj < padj_cutoff) %>% dplyr::arrange(padj)    
write.csv(sig_res_liver_PB, paste0("~/Desktop/UBC/10x/NovaSeq/cellranger_count_rawfullGTF/05.differential_gene_expression/liver_filtered/DESeq2_filterlowexpr/pseudobulk/deseq_beforeshrink_sig_genes.csv"), quote = FALSE, row.names = FALSE)

sig_res_liver_PB_ashr <- dplyr::filter(res_tbl_liver_PB_ashr, padj < padj_cutoff) %>% dplyr::arrange(padj)    
write.csv(sig_res_liver_PB_ashr, paste0("~/Desktop/UBC/10x/NovaSeq/cellranger_count_rawfullGTF/05.differential_gene_expression/liver_filtered/DESeq2_filterlowexpr/pseudobulk/deseq_shrinkashr_sig_genes.csv"), quote = FALSE, row.names = FALSE)

sig_up_fc_PB_liver <- dplyr::filter(sig_res_liver_PB, log2FoldChange >= log2fc_cutoff)
sig_down_fc_PB_liver <- dplyr::filter(sig_res_liver_PB, log2FoldChange <= -log2fc_cutoff)
write.csv(sig_up_fc_PB_liver, paste0("~/Desktop/UBC/10x/NovaSeq/cellranger_count_rawfullGTF/05.differential_gene_expression/liver_filtered/DESeq2_filterlowexpr/pseudobulk/deseq_beforeshrink_sig_genes_fcpass_malebiased.csv"), quote = FALSE, row.names = FALSE)
write.csv(sig_down_fc_PB_liver, paste0("~/Desktop/UBC/10x/NovaSeq/cellranger_count_rawfullGTF/05.differential_gene_expression/liver_filtered/DESeq2_filterlowexpr/pseudobulk/deseq_beforeshrink_sig_genes_fcpass_femalebiased.csv"), quote = FALSE, row.names = FALSE)

sig_up_fc_PB_ashr_liver <- dplyr::filter(sig_res_liver_PB_ashr, log2FoldChange >= log2fc_cutoff)
sig_down_fc_PB_ashr_liver <- dplyr::filter(sig_res_liver_PB_ashr, log2FoldChange <= -log2fc_cutoff)
write.csv(sig_up_fc_PB_ashr_liver, paste0("~/Desktop/UBC/10x/NovaSeq/cellranger_count_rawfullGTF/05.differential_gene_expression/liver_filtered/DESeq2_filterlowexpr/pseudobulk/deseq_shrinkashr_sig_genes_fcpass_malebiased.csv"), quote = FALSE, row.names = FALSE)
write.csv(sig_down_fc_PB_ashr_liver, paste0("~/Desktop/UBC/10x/NovaSeq/cellranger_count_rawfullGTF/05.differential_gene_expression/liver_filtered/DESeq2_filterlowexpr/pseudobulk/deseq_shrinkashr_sig_genes_fcpass_femalebiased.csv"), quote = FALSE, row.names = FALSE)

## ggplot of top genes
normalized_counts_liver_PB <- counts(dds_deseq_liver_PB, normalized = TRUE)
write.csv(normalized_counts_liver_PB, paste0("~/Desktop/UBC/10x/NovaSeq/cellranger_count_rawfullGTF/05.differential_gene_expression/liver_filtered/DESeq2_filterlowexpr/pseudobulk/dds_normalized_unlogged_counts_PB.csv"), quote = FALSE, row.names = T)       
    
## Order results by padj values
top20_sig_genes_liver_PB_ashr <- sig_res_liver_PB_ashr %>%
                dplyr::arrange(padj) %>%
                dplyr::pull(gene) %>%
                head(n=20)

top20_sig_norm_liver_PB_ashr <- data.frame(normalized_counts_liver_PB) %>%
                rownames_to_column(var = "gene") %>%
                dplyr::filter(gene %in% top20_sig_genes_liver_PB_ashr)
        
gathered_top20_sig_liver_PB_ashr <- top20_sig_norm_liver_PB_ashr %>%
                gather(colnames(top20_sig_norm_liver_PB_ashr)[2:length(colnames(top20_sig_norm_liver_PB_ashr))], key = "samplename", value = "normalized_counts")
        
gathered_top20_sig_liver_PB_ashr <- inner_join(ei_liver_PB[, c("sample_id", "sex_id" )], gathered_top20_sig_liver_PB_ashr, by = c("sample_id" = "samplename"))
        
## plot using ggplot2    
ggplot(gathered_top20_sig_liver_PB_ashr) +
                geom_point(aes(x = gene, 
                     y = normalized_counts, 
                     color = sex_id), 
                     position=position_jitter(w=0.1,h=0)) +
                scale_y_log10() +
                xlab("Genes") +
                ylab("log10 Normalized Counts") +
                ggtitle("Top 20 Significant DE Genes") +
                theme_bw() +
                theme(axis.text.x = element_text(angle = 45, hjust = 1)) + 
                theme(plot.title = element_text(hjust = 0.5))
ggsave(paste0("~/Desktop/UBC/10x/NovaSeq/cellranger_count_rawfullGTF/05.differential_gene_expression/liver_filtered/DESeq2_filterlowexpr/pseudobulk/plots/pseudobulk_shrinkashr_top20_DE_genes.png"))  

        str("volcano plot")
		res_tbl_thresh <- res_tbl_liver_PB_ashr[!is.na(res_tbl_liver_PB_ashr$padj), ] %>% mutate(threshold=padj<0.05 & abs(log2FoldChange) >=1)
		min(log10(res_tbl_thresh$padj))
		
		ggplot(res_tbl_thresh) +
			geom_point(aes(x=log2FoldChange, y=-log10(padj), colour=threshold)) +
			ggtitle("Volcano plot") +
			xlab("log2 fold change") +
			ylab("-log10 adjusted p-value")+
			scale_color_manual(values=c("grey60", "red3")) +
			theme(legend.position="none", plot.title=element_text(size=rel(1.3), hjust=0.5), axis.title=element_text(size=rel(1.15)))
		ggsave(paste0("~/Desktop/UBC/10x/NovaSeq/cellranger_count_rawfullGTF/05.differential_gene_expression/liver_filtered/DESeq2_filterlowexpr/pseudobulk/plots/pseudobulk_volcano_plot.png"))



