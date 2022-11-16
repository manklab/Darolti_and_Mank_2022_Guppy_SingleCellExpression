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

load("~/Desktop/UBC/10x/NovaSeq/cellranger_count_rawfullGTF/03.marker_identification/gonad_filtered_dbrmv_reevaluated/gonad_umap_3_rename_clusters.RData")

gonad_umap_3_rename_clusters_dbrmv_PB <- gonad_umap_3_rename_clusters
gonad_umap_3_rename_clusters_dbrmv_PB$CellType <- "All"

head(gonad_umap_3_rename_clusters_dbrmv_PB)

counts_gonad_PB <- gonad_umap_3_rename_clusters_dbrmv_PB@assays$RNA@counts
metadata_gonad_PB <- gonad_umap_3_rename_clusters_dbrmv_PB@meta.data
metadata_gonad_PB$CellType <- "All"

sce_gonad_PB <- SingleCellExperiment(assays=list(counts=counts_gonad_PB), colData=metadata_gonad_PB)
sce_gonad_PB
class: SingleCellExperiment 
dim: 21793 28014 
metadata(0):
assays(1): counts
rownames(21793): ENSPREG00000006268 ENSPREG00000006273 ... ENSPREG00000011737 BDKRB2.2
rowData names(0):
colnames(28014): F1_Lib11_AAACCCAAGGTACAGC-1 F1_Lib11_AAACCCAAGTCTGGAG-1 ...
  M6_Lib37_TTTGTTGTCGCACGAC-1 M6_Lib37_TTTGTTGTCTATACTC-1
colData names(17): seq_folder nUMI ... seurat_clusters CellType
reducedDimNames(0):
altExpNames(0):

assays(sce_gonad_PB)
		List of length 1
		names(1): counts
dim(counts(sce_gonad_PB))
[1] 21793 28014
counts(sce_gonad_PB)[1:2, 1:2]
		2 x 2 sparse Matrix of class "dgCMatrix"
		                   F1_Lib8_AAACCCAAGAGTCTTC-1 F1_Lib8_AAACCCAAGGATTTCC-1
		ENSPREG00000006268                          .                          .
		arl13b                                      .                          .
dim(colData(sce_gonad_PB))
[1] 28014    17

# Determine the number of clusters and the cluster names
colData(sce_gonad_PB)$sample_id <- as.factor(colData(sce_gonad_PB)$sample)
colData(sce_gonad_PB)$sex_id <- as.factor(colData(sce_gonad_PB)$sex)
colData(sce_gonad_PB)$cluster_id <- as.factor(colData(sce_gonad_PB)$CellType)
kids_gonad_PB <- purrr::set_names(levels(sce_gonad_PB$cluster_id))
kids_gonad_PB
		  All 
		"All" 
nk_gonad_PB <- length(kids_gonad_PB)
nk_gonad_PB
		[1] 1

sids_gonad_PB <- purrr::set_names(levels(sce_gonad_PB$sample_id))
sids_gonad_PB
		  GonadF1   GonadF2   GonadF3   GonadM1   GonadM2   GonadM3 
		"GonadF1" "GonadF2" "GonadF3" "GonadM1" "GonadM2" "GonadM3" 
ns_gonad_PB <- length(sids_gonad_PB)
ns_gonad_PB
		[1] 6

# Determine no. cells per sample
table(sce_gonad_PB$sample_id)
		GonadF1 GonadF2 GonadF3 GonadM1 GonadM2 GonadM3 
   3162    2063    1492    7973    6541    6783

# Turn named vector into a numeric vector
n_cells_gonad_PB <- as.numeric(table(sce_gonad_PB$sample_id))
 
# Reorder samples (rows) of the metadata to match the order of the sample names
m_gonad_PB <- match(sids_gonad_PB, sce_gonad_PB$sample_id)
m_gonad_PB
[1]     1  3163  5226  6718 14691 21232

# Create sample level metadata by combining the reordered metadata with the number of cells
ei_gonad_PB <- data.frame(colData(sce_gonad_PB)[m_gonad_PB, ], n_cells_gonad_PB, row.names=NULL) %>% select(-"cluster_id")

# Remove lowly expressed genes
dim(sce_gonad_PB)
[1] 21793 28014
sce_gonad_PB <- sce_gonad_PB[rowSums(counts(sce_gonad_PB)>1)>=10,]
dim(sce_gonad_PB)
[1] 13617 28014

# Aggregate counts per sample_id and cluster_id
groups_gonad_PB <- colData(sce_gonad_PB)[, c("cluster_id", "sample_id")]
pb_gonad_PB <- aggregate.Matrix(t(counts(sce_gonad_PB)), groupings=groups_gonad_PB, fun="sum")
class(pb_gonad_PB)
		[1] "dgCMatrix"
		attr(,"package")
		[1] "Matrix"
dim(pb_gonad_PB)
[1]     6 13617

pb_gonad_PB[1:6, 1:6]
6 x 6 sparse Matrix of class "dgCMatrix"
            arl13b zgc:152904 gart n6amt1 timmdc1 VANGL1
All_GonadF1     22         86   98     19      93     43
All_GonadF2     46         85  146     43      61     59
All_GonadF3     74         44  291     25      56     38
All_GonadM4     93        464   47     43      58     31
All_GonadM5    127        377   65     53      60     62
All_GonadM6    135        614   52     63      67     24


## Split data by cell type and Transform matrix so that genes are row names and samples are column names
splitf_gonad_PB <- sapply(stringr::str_split(rownames(pb_gonad_PB), pattern="_", n=2), '[', 1) # Check that each cluster is present in all samples
pb_gonad_PB <- split.data.frame(pb_gonad_PB, factor(splitf_gonad_PB)) %>% lapply(function(u) set_colnames(t(u), stringr::str_extract(rownames(u), "(?<=_)[:alnum:]+")))

class(pb_gonad_PB)
		[1] "list"
str(pb_gonad_PB)

options(width=100)
table(sce_gonad_PB$cluster_id, sce_gonad_PB$sample_id)
		      GonadF1 GonadF2 GonadF3 GonadM1 GonadM2 GonadM3
  All    3162    2063    1492    7973    6541    6783
  
get_sample_ids_gonad_PB <- function(x){
	pb_gonad_PB[[x]] %>% colnames()
}
de_samples_gonad_PB <- map(1:length(kids_gonad_PB), get_sample_ids_gonad_PB) %>% unlist()

samples_list_gonad_PB <- map(1:length(kids_gonad_PB), get_sample_ids_gonad_PB)

get_cluster_ids_gonad_PB <- function(x){
	rep(names(pb_gonad_PB)[x], each=length(samples_list_gonad_PB[[x]]))
}

de_cluster_ids_gonad_PB <- map(1:length(kids_gonad_PB), get_cluster_ids_gonad_PB) %>% unlist()

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
write.csv(raw_counts_PB, paste0("~/Desktop/UBC/10x/NovaSeq/cellranger_count_rawfullGTF/05.differential_gene_expression/gonad_filtered/DESeq2_filterlowexpr/pseudobulk/dds_raw_counts_PB.csv"), quote = FALSE, row.names = T)

#-----------
fpm_PB <- fpm(dds_gonad_PB)
head(fpm_PB)
             GonadF1   GonadF2   GonadF3   GonadM4   GonadM5   GonadM6
arl13b      2.568094  6.274860 11.250085 12.150384 15.495140 16.088398
zgc:152904 10.038912 11.594849  6.689240 60.621269 45.997383 73.172417
gart       11.439690 19.915858 44.240200  6.140516  7.930583  6.197013
n6amt1      2.217899  5.865630  3.800704  5.617919  6.466476  7.507919
timmdc1    10.856032  8.321009  8.513578  7.577659  7.320538  7.984612
VANGL1      5.019456  8.048189  5.777071  4.050128  7.564556  2.860160

write.csv(fpm_PB, paste0("~/Desktop/UBC/10x/NovaSeq/cellranger_count_rawfullGTF/05.differential_gene_expression/gonad_filtered/DESeq2_filterlowexpr/pseudobulk/dds_fpm_PB.csv"), quote = FALSE, row.names = T)

cpm_PB <- calculateCPM(dds_gonad_PB)
head(cpm_PB)
            GonadF1   GonadF2   GonadF3   GonadM4   GonadM5   GonadM6
arl13b     1.931988  6.252001 12.233407 13.426540 18.675056 14.450295
zgc:152904 7.552318 11.552611  7.273918 66.988326 55.436978 65.722081
gart       8.606129 19.843309 48.107046  6.785455  9.558100  5.566039
n6amt1     1.668535  5.844262  4.132908  6.207970  7.793527  6.743471
timmdc1    8.167041  8.290697  9.257713  8.373541  8.822861  7.171628
VANGL1     3.776159  8.018871  6.282020  4.475513  9.116957  2.568941

write.csv(cpm_PB, paste0("~/Desktop/UBC/10x/NovaSeq/cellranger_count_rawfullGTF/05.differential_gene_expression/gonad_filtered/DESeq2_filterlowexpr/pseudobulk/dds_cpm_PB.csv"), quote = FALSE, row.names = T)
#--------------


# Transform counts for data visualization
rld_gonad_PB <- rlog(dds_gonad_PB, blind=TRUE)
colData(rld_gonad_PB)
DataFrame with 6 rows and 5 columns
        cluster_id   sample_id   sex_id n_cells_gonad_PB sizeFactor
          <factor> <character> <factor>        <numeric>  <numeric>
GonadF1        All     GonadF1   Female             3162   1.100066
GonadF2        All     GonadF2   Female             2063   0.941371
GonadF3        All     GonadF3   Female             1492   0.844662
GonadM4        All     GonadM4   Male               7973   0.982878
GonadM5        All     GonadM5   Male               6541   1.052483
GonadM6        All     GonadM6   Male               6783   1.077526
        
# Plot PCA
DESeq2::plotPCA(rld_gonad_PB, ntop=100,intgroup = "sex_id")
ggsave(paste0("~/Desktop/UBC/10x/NovaSeq/cellranger_count_rawfullGTF/05.differential_gene_expression/gonad_filtered/DESeq2_filterlowexpr/pseudobulk/plots/pseudobulk_PCAplot.png"))

        
        DESeq2::plotPCA(rld_gonad_PB, ntop=100, intgroup = "n_cells_gonad_PB")
        ggsave(paste0("~/Desktop/UBC/10x/NovaSeq/cellranger_count_rawfullGTF/05.differential_gene_expression/gonad_filtered/DESeq2_filterlowexpr/pseudobulk/plots/pseudobulk_PCAplot_cellcount.png"))

# Extract the rlog matrix from the object and compute pairwise correlation values
rld_mat_gonad_PB <- assay(rld_gonad_PB)
head(rld_mat_gonad_PB)
            GonadF1  GonadF2  GonadF3  GonadM4  GonadM5  GonadM6
arl13b     4.891621 5.764233 6.377669 6.460649 6.725301 6.766562
zgc:152904 6.635636 6.781848 6.237558 8.572965 8.261750 8.787434
gart       6.502935 7.103255 8.000809 5.859975 6.119359 5.868895
n6amt1     4.471423 5.442231 4.996969 5.397194 5.544770 5.703147
timmdc1    6.300746 6.019775 6.043833 5.921900 5.886065 5.976474
VANGL1     5.308096 5.804175 5.453666 5.089602 5.738009 4.747944

write.csv(rld_mat_gonad_PB, paste0("~/Desktop/UBC/10x/NovaSeq/cellranger_count_rawfullGTF/05.differential_gene_expression/gonad_filtered/DESeq2_filterlowexpr/pseudobulk/rlogtransformed_counts_PB.csv"), quote = FALSE, row.names = T)
rld_cor_gonad_PB <- cor(rld_mat_gonad_PB)
write.csv(rld_cor_gonad_PB, paste0("~/Desktop/UBC/10x/NovaSeq/cellranger_count_rawfullGTF/05.differential_gene_expression/gonad_filtered/DESeq2_filterlowexpr/pseudobulk/deseq2_pairwise_correlations_PB.csv"), quote = FALSE, row.names = T)
        
# Plot heatmap
png(paste0("~/Desktop/UBC/10x/NovaSeq/cellranger_count_rawfullGTF/05.differential_gene_expression/gonad_filtered/DESeq2_filterlowexpr/pseudobulk/plots/pseudobulk_heatmap.png"), height=6, width=7.5, units="in", res=300)
pheatmap(rld_cor_gonad_PB, annotation = cluster_metadata_gonad_PB[, c("sex_id"), drop=F])
dev.off()
        
# Run DESeq2 differential expression analysis
dds_deseq_gonad_PB <- DESeq(dds_gonad_PB)
#estimating size factors
#estimating dispersions
#gene-wise dispersion estimates
#mean-dispersion relationship
#final dispersion estimates
#fitting model and testing
        
# Plot dispersion estimates
png(paste0("~/Desktop/UBC/10x/NovaSeq/cellranger_count_rawfullGTF/05.differential_gene_expression/gonad_filtered/DESeq2_filterlowexpr/pseudobulk/plots/pseudobulk_dispersion_plot.png"))
plotDispEsts(dds_deseq_gonad_PB)
dev.off()

# Output results of Wald test for contrast for A vs B
contrast_gonad_PB <- c("sex_id", levels(cluster_metadata_gonad_PB$sex_id)[2], levels(cluster_metadata_gonad_PB$sex_id)[1])
        [1] "sex_id" "Male"   "Female"
# resultsNames(dds)
res_gonad_PB <- results(dds_deseq_gonad_PB, 
                contrast = contrast_gonad_PB,
                alpha = 0.05)

res_gonad_PB_ashr <- lfcShrink(dds_deseq_gonad_PB, 
                contrast=contrast_gonad_PB,
                res=res_gonad_PB, type="ashr")
                     
 # Try other adjsuting methods
 res_gonad_PB_bon <- results(dds_deseq_gonad_PB, 
                contrast = contrast_gonad_PB,
                alpha = 0.05,
                pAdjustMethod="bonferroni")                   
                     
                      
# Set thresholds
padj_cutoff <- 0.05
log2fc_cutoff <- 1
        
# Turn the results object into a tibble for use with tidyverse functions

res_tbl_gonad_PB <- res_gonad_PB %>%
                 data.frame() %>%
                 rownames_to_column(var="gene") %>%
                 as_tibble()
write.csv(res_tbl_gonad_PB, paste0("~/Desktop/UBC/10x/NovaSeq/cellranger_count_rawfullGTF/05.differential_gene_expression/gonad_filtered/DESeq2_filterlowexpr/pseudobulk/deseq_beforeshrink_all_genes.csv"), quote = FALSE, row.names = FALSE)


res_tbl_gonad_PB_ashr <- res_gonad_PB_ashr %>%
                 data.frame() %>%
                 rownames_to_column(var="gene") %>%
                 as_tibble()
write.csv(res_tbl_gonad_PB_ashr, paste0("~/Desktop/UBC/10x/NovaSeq/cellranger_count_rawfullGTF/05.differential_gene_expression/gonad_filtered/DESeq2_filterlowexpr/pseudobulk/deseq_shrinkashr_all_genes.csv"), quote = FALSE, row.names = FALSE)
        
# Subset the significant results
sig_res_gonad_PB <- dplyr::filter(res_tbl_gonad_PB, padj < padj_cutoff) %>% dplyr::arrange(padj)    
write.csv(sig_res_gonad_PB, paste0("~/Desktop/UBC/10x/NovaSeq/cellranger_count_rawfullGTF/05.differential_gene_expression/gonad_filtered/DESeq2_filterlowexpr/pseudobulk/deseq_beforeshrink_sig_genes.csv"), quote = FALSE, row.names = FALSE)

sig_res_gonad_PB_ashr <- dplyr::filter(res_tbl_gonad_PB_ashr, padj < padj_cutoff) %>% dplyr::arrange(padj)    
write.csv(sig_res_gonad_PB_ashr, paste0("~/Desktop/UBC/10x/NovaSeq/cellranger_count_rawfullGTF/05.differential_gene_expression/gonad_filtered/DESeq2_filterlowexpr/pseudobulk/deseq_shrinkashr_sig_genes.csv"), quote = FALSE, row.names = FALSE)

sig_up_fc_PB_gonad <- dplyr::filter(sig_res_gonad_PB, log2FoldChange >= log2fc_cutoff)
sig_down_fc_PB_gonad <- dplyr::filter(sig_res_gonad_PB, log2FoldChange <= -log2fc_cutoff)
write.csv(sig_up_fc_PB_gonad, paste0("~/Desktop/UBC/10x/NovaSeq/cellranger_count_rawfullGTF/05.differential_gene_expression/gonad_filtered/DESeq2_filterlowexpr/pseudobulk/deseq_beforeshrink_sig_genes_fcpass_malebiased.csv"), quote = FALSE, row.names = FALSE)
write.csv(sig_down_fc_PB_gonad, paste0("~/Desktop/UBC/10x/NovaSeq/cellranger_count_rawfullGTF/05.differential_gene_expression/gonad_filtered/DESeq2_filterlowexpr/pseudobulk/deseq_beforeshrink_sig_genes_fcpass_femalebiased.csv"), quote = FALSE, row.names = FALSE)

sig_up_fc_PB_ashr_gonad <- dplyr::filter(sig_res_gonad_PB_ashr, log2FoldChange >= log2fc_cutoff)
sig_down_fc_PB_ashr_gonad <- dplyr::filter(sig_res_gonad_PB_ashr, log2FoldChange <= -log2fc_cutoff)
write.csv(sig_up_fc_PB_ashr_gonad, paste0("~/Desktop/UBC/10x/NovaSeq/cellranger_count_rawfullGTF/05.differential_gene_expression/gonad_filtered/DESeq2_filterlowexpr/pseudobulk/deseq_shrinkashr_sig_genes_fcpass_malebiased.csv"), quote = FALSE, row.names = FALSE)
write.csv(sig_down_fc_PB_ashr_gonad, paste0("~/Desktop/UBC/10x/NovaSeq/cellranger_count_rawfullGTF/05.differential_gene_expression/gonad_filtered/DESeq2_filterlowexpr/pseudobulk/deseq_shrinkashr_sig_genes_fcpass_femalebiased.csv"), quote = FALSE, row.names = FALSE)

## ggplot of top genes
normalized_counts_gonad_PB <- counts(dds_deseq_gonad_PB, normalized = TRUE)
write.csv(normalized_counts_gonad_PB, paste0("~/Desktop/UBC/10x/NovaSeq/cellranger_count_rawfullGTF/05.differential_gene_expression/gonad_filtered/DESeq2_filterlowexpr/pseudobulk/dds_normalized_unlogged_counts_PB.csv"), quote = FALSE, row.names = T)       
    
## Order results by padj values
top20_sig_genes_gonad_PB_ashr <- sig_res_gonad_PB_ashr %>%
                dplyr::arrange(padj) %>%
                dplyr::pull(gene) %>%
                head(n=20)

top20_sig_norm_gonad_PB_ashr <- data.frame(normalized_counts_gonad_PB) %>%
                rownames_to_column(var = "gene") %>%
                dplyr::filter(gene %in% top20_sig_genes_gonad_PB_ashr)
        
gathered_top20_sig_gonad_PB_ashr <- top20_sig_norm_gonad_PB_ashr %>%
                gather(colnames(top20_sig_norm_gonad_PB_ashr)[2:length(colnames(top20_sig_norm_gonad_PB_ashr))], key = "samplename", value = "normalized_counts")
        
gathered_top20_sig_gonad_PB_ashr <- inner_join(ei_gonad_PB[, c("sample_id", "sex_id" )], gathered_top20_sig_gonad_PB_ashr, by = c("sample_id" = "samplename"))
        
## plot using ggplot2    
ggplot(gathered_top20_sig_gonad_PB_ashr) +
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
ggsave(paste0("~/Desktop/UBC/10x/NovaSeq/cellranger_count_rawfullGTF/05.differential_gene_expression/gonad_filtered/DESeq2_filterlowexpr/pseudobulk/plots/pseudobulk_shrinkashr_top20_DE_genes.png"))  

        str("volcano plot")
		res_tbl_thresh <- res_tbl_gonad_PB_ashr[!is.na(res_tbl_gonad_PB_ashr$padj), ] %>% mutate(threshold=padj<0.05 & abs(log2FoldChange) >=1)
		min(log10(res_tbl_thresh$padj))
		
		ggplot(res_tbl_thresh) +
			geom_point(aes(x=log2FoldChange, y=-log10(padj), colour=threshold)) +
			ggtitle("Volcano plot") +
			xlab("log2 fold change") +
			ylab("-log10 adjusted p-value")+
			scale_color_manual(values=c("grey60", "red3")) +
			theme(legend.position="none", plot.title=element_text(size=rel(1.3), hjust=0.5), axis.title=element_text(size=rel(1.15)))
		ggsave(paste0("~/Desktop/UBC/10x/NovaSeq/cellranger_count_rawfullGTF/05.differential_gene_expression/gonad_filtered/DESeq2_filterlowexpr/pseudobulk/plots/pseudobulk_volcano_plot.png"))



