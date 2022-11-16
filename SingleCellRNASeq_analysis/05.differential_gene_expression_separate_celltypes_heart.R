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

load("~/Desktop/UBC/10x/NovaSeq/cellranger_count_rawfullGTF/04.remove_doublets/heart_filtered/heart_umap_3_rename_clusters_dbrmv.RData")

head(heart_umap_3_rename_clusters_dbrmv)

counts_heart <- heart_umap_3_rename_clusters_dbrmv@assays$RNA@counts
metadata_heart <- heart_umap_3_rename_clusters_dbrmv@meta.data
metadata_heart$CellType <- factor(heart_umap_3_rename_clusters_dbrmv@active.ident)

sce_heart <- SingleCellExperiment(assays=list(counts=counts_heart), colData=metadata_heart)
sce_heart 
class: SingleCellExperiment 
dim: 21975 28923 
metadata(0):
assays(1): counts
rownames(21975): ENSPREG00000006268 arl13b ... guca1d ENSPREG00000006408
rowData names(0):
colnames(28923): F1_Lib6_AAACCCACAACCTATG-1 F1_Lib6_AAACCCACAACGGCTC-1 ...
  M3_Lib27_TTTGGTTGTATATGGA-1 M3_Lib27_TTTGGTTGTGTTACTG-1
colData names(27): seq_folder nUMI ... seurat_clusters CellType
reducedDimNames(0):
altExpNames(0):

assays(sce_heart)
		List of length 1
		names(1): counts
dim(counts(sce_heart))
[1] 21975 28923
counts(sce_heart)[1:2, 1:2]
		2 x 2 sparse Matrix of class "dgCMatrix"
		                   F1_Lib8_AAACCCAAGAGTCTTC-1 F1_Lib8_AAACCCAAGGATTTCC-1
		ENSPREG00000006268                          .                          .
		arl13b                                      .                          .
dim(colData(sce_heart))
[1] 28923    27


# Determine the number of clusters (nk) and the cluster names (kids), and the number of samples (ns) and sample names (sids)
colData(sce_heart)$sample_id <- as.factor(colData(sce_heart)$sample)
colData(sce_heart)$sex_id <- as.factor(colData(sce_heart)$sex)
colData(sce_heart)$cluster_id <- as.factor(colData(sce_heart)$CellType)

kids_heart <- purrr::set_names(levels(sce_heart$cluster_id))
kids_heart
              Fibroblast               Macrophage          Cardiomyocyte 1              Erythrocyte 
            "Fibroblast"             "Macrophage"        "Cardiomyocyte 1"            "Erythrocyte" 
            T lymphocyte             B lymphocyte              Endocardial      Vacular endothelial 
          "T lymphocyte"           "B lymphocyte"            "Endocardial"    "Vacular endothelial" 
         Cardiomyocyte 2              Granulocyte   Vascular smooth muscle 
       "Cardiomyocyte 2"            "Granulocyte" "Vascular smooth muscle" 

nk_heart <- length(kids_heart)
nk_heart
		[1] 11

sids_heart <- purrr::set_names(levels(sce_heart$sample_id))
sids_heart
		   HeartF1   HeartF2   HeartF3   HeartM1   HeartM2   HeartM3 
		 "HeartF1" "HeartF2" "HeartF3" "HeartM1" "HeartM2" "HeartM3" 
ns_heart <- length(sids_heart)
ns_heart
		[1] 6

# Number of cells per sample
table(sce_heart$sample_id)
HeartF1 HeartF2 HeartF3 HeartM1 HeartM2 HeartM3 
   6735    8966    6150    4615    1163    1294
  
# Turn named vector into a numeric vector
n_cells_heart <- as.numeric(table(sce_heart$sample_id))
 
# Reorder samples (rows) of the metadata to match the order of the sample names
m_heart <- match(sids_heart, sce_heart$sample_id)
m_heart
[1]     1  6736 15702 21852 26467 27630

# Create sample level metadata by combining the reordered metadata with the number of cells
ei_heart <- data.frame(colData(sce_heart)[m_heart, ], n_cells_heart, row.names=NULL) %>% select(-"cluster_id")

# Remove lowly expressed genes
dim(sce_heart)
[1] 21975 28923
sce_heart <- sce_heart[rowSums(counts(sce_heart)>0)>0,]
dim(sce_heart)	
[1] 18276 28923
sce_heart <- sce_heart[rowSums(counts(sce_heart)>1)>=10,]
dim(sce_heart)
[1]  9860 28923		

# Aggregate counts per sample_id and cluster_id
groups_heart <- colData(sce_heart)[, c("cluster_id", "sample_id")]
pb_heart <- aggregate.Matrix(t(counts(sce_heart)), groupings=groups_heart, fun="sum")
class(pb_heart)
		[1] "dgCMatrix"
		attr(,"package")
		[1] "Matrix"
dim(pb_heart)
[1]   66 9860

pb_heart[1:6, 1:6]
		6 x 6 sparse Matrix of class "dgCMatrix"
                        zgc:152904 gart VANGL1 sap18 rtraf gpr137c
Fibroblast_HeartF1              21   14      5    84    11      26
Macrophage_HeartF1               2   23      4   392    67      54
Cardiomyocyte 1_HeartF1          3    .      2    44     5       4
Erythrocyte_HeartF1              2    .      1    45     7      37
T lymphocyte_HeartF1             1    .      1   139    17      15
B lymphocyte_HeartF1             .    3     36   144    19       7


## Split data by cell type and Transform matrix so that genes are row names and samples are column names
splitf_heart <- sapply(stringr::str_split(rownames(pb_heart), pattern="_", n=2), '[', 1) # Check that each cluster is present in all samples
pb_heart <- split.data.frame(pb_heart, factor(splitf_heart)) %>% lapply(function(u) set_colnames(t(u), stringr::str_extract(rownames(u), "(?<=_)[:alnum:]+")))
class(pb_heart)
		[1] "list"
str(pb_heart)

# Number of cells for each cell type and sample
options(width=100)
table(sce_heart$cluster_id, sce_heart$sample_id)
                         HeartF1 HeartF2 HeartF3 HeartM1 HeartM2 HeartM3
  Fibroblast                1791    3016     512     736     372     156
  Macrophage                1121     968    1053     713     141     111
  Cardiomyocyte 1            677    1596    1268     649     109     232
  Erythrocyte               1164    1380    1366    1127     193     410
  T lymphocyte               892    1252     562     496      91     114
  B lymphocyte               690     193     420     171      13      10
  Endocardial                135     217     364     253      58     126
  Vacular endothelial         22      46     376     247      16      14
  Cardiomyocyte 2            128     114     110     143     127      94
  Granulocyte                 32     147      25      59       6      15
  Vascular smooth muscle      83      37      94      21      37      12


get_sample_ids_heart <- function(x){
	pb_heart[[x]] %>% colnames()
}
de_samples_heart <- map(1:length(kids_heart), get_sample_ids_heart) %>% unlist()

samples_list_heart <- map(1:length(kids_heart), get_sample_ids_heart)

get_cluster_ids_heart <- function(x){
	rep(names(pb_heart)[x], each=length(samples_list_heart[[x]]))
}

de_cluster_ids_heart <- map(1:length(kids_heart), get_cluster_ids_heart) %>% unlist()

gg_df_heart <- data.frame(cluster_id=de_cluster_ids_heart, sample_id=de_samples_heart)
gg_df_heart <- left_join(gg_df_heart, ei_heart[, c("sample_id", "sex_id", "n_cells_heart")])
metadata_heart <- gg_df_heart %>% dplyr::select(cluster_id, sample_id, sex_id, n_cells_heart)

metadata_heart$cluster_id <- as.factor(metadata_heart$cluster_id)

clusters_heart <- levels(metadata_heart$cluster_id)

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
        write.csv(raw_counts, paste0("~/Desktop/UBC/10x/NovaSeq/cellranger_count_rawfullGTF/05.differential_gene_expression/heart_filtered/DESeq2_filterlowexpr/separate_clusters/",clusters_heart[x],"_dds_raw_counts.csv"), quote = FALSE, row.names = T)
        
        fpm_heart <- fpm(dds_heart)
        write.csv(fpm_heart, paste0("~/Desktop/UBC/10x/NovaSeq/cellranger_count_rawfullGTF/05.differential_gene_expression/heart_filtered/DESeq2_filterlowexpr/separate_clusters/", clusters_heart[x], "_dds_fpm.csv"), quote = FALSE, row.names = T)
        
        cpm_heart <- calculateCPM(dds_heart)
        write.csv(cpm_heart, paste0("~/Desktop/UBC/10x/NovaSeq/cellranger_count_rawfullGTF/05.differential_gene_expression/heart_filtered/DESeq2_filterlowexpr/separate_clusters/", clusters_heart[x], "_dds_cpm.csv"), quote = FALSE, row.names = T)
        
        
        str("Transform counts for data visualization")
        rld_heart <- rlog(dds_heart, blind=TRUE)
        
        str("Plot PCA")
        DESeq2::plotPCA(rld_heart, ntop=100, intgroup = "sex_id")
        ggsave(paste0("~/Desktop/UBC/10x/NovaSeq/cellranger_count_rawfullGTF/05.differential_gene_expression/heart_filtered/DESeq2_filterlowexpr/separate_clusters/plots/", clusters_heart[x], "_specific_PCAplot.png"))
        
        DESeq2::plotPCA(rld_heart, ntop=100, intgroup = "n_cells_heart")
        ggsave(paste0("~/Desktop/UBC/10x/NovaSeq/cellranger_count_rawfullGTF/05.differential_gene_expression/heart_filtered/DESeq2_filterlowexpr/separate_clusters/plots/", clusters_heart[x], "_specific_PCAplot_cellcount.png"))
        
        str("Extract the rlog matrix from the object and compute pairwise correlation values")
        rld_mat_heart <- assay(rld_heart)
        write.csv(rld_mat_heart, paste0("~/Desktop/UBC/10x/NovaSeq/cellranger_count_rawfullGTF/05.differential_gene_expression/heart_filtered/DESeq2_filterlowexpr/separate_clusters/", clusters_heart[x], "_rlogtransformed_counts.csv"), quote = FALSE, row.names = T)
        rld_cor_heart <- cor(rld_mat_heart)
        write.csv(rld_cor_heart, paste0("~/Desktop/UBC/10x/NovaSeq/cellranger_count_rawfullGTF/05.differential_gene_expression/heart_filtered/DESeq2_filterlowexpr/separate_clusters/", clusters_heart[x], "_pairwise_correlations.csv"), quote = FALSE, row.names = T)
        
        
        str("Plot heatmap")
        png(paste0("~/Desktop/UBC/10x/NovaSeq/cellranger_count_rawfullGTF/05.differential_gene_expression/heart_filtered/DESeq2_filterlowexpr/separate_clusters/plots/", clusters_heart[x], "_specific_heatmap_correl.png"), height=6, width=7.5, units="in", res=300)
        pheatmap(rld_cor_heart, annotation = cluster_metadata_heart[, c("sex_id"), drop=F])
        dev.off()

        str("Run DESeq2 differential expression analysis")
        dds_deseq_heart <- DESeq(dds_heart)
        
        str("Plot dispersion estimates")
        png(paste0("~/Desktop/UBC/10x/NovaSeq/cellranger_count_rawfullGTF/05.differential_gene_expression/heart_filtered/DESeq2_filterlowexpr/separate_clusters/plots/", clusters_heart[x], "_dispersion_plot.png"), height=4, width=6, units="in", res=300)
        plotDispEsts(dds_deseq_heart)
        dev.off()

        str("Output results of Wald test for contrast for A vs B")
        contrast_heart <- c("sex_id", levels(cluster_metadata_heart$sex_id)[A], levels(cluster_metadata_heart$sex_id)[B])
        str(contrast_heart)
        
        str("resultsNames(dds)")
        res_heart <- results(dds_deseq_heart, 
                       contrast = contrast_heart,
                       alpha = 0.05)
        
        res_heart_ashr <- lfcShrink(dds_deseq_heart, 
                         contrast=contrast_heart,
                         res=res_heart, type="ashr")
        str("Set thresholds")
        padj_cutoff <- 0.05
        log2fc_cutoff <- 1
        
        str("Turn the results object into a tibble for use with tidyverse functions")
        res_tbl_heart <- res_heart %>%
                data.frame() %>%
                rownames_to_column(var="gene") %>%
                as_tibble()
         write.csv(res_tbl_heart, paste0("~/Desktop/UBC/10x/NovaSeq/cellranger_count_rawfullGTF/05.differential_gene_expression/heart_filtered/DESeq2_filterlowexpr/separate_clusters/", clusters_heart[x], "_deseq_beforeshrink_all_genes.csv"), quote = FALSE, row.names = FALSE)
         
        res_tbl_heart_ashr <- res_heart_ashr %>%
                data.frame() %>%
                rownames_to_column(var="gene") %>%
                as_tibble()
        write.csv(res_tbl_heart_ashr, paste0("~/Desktop/UBC/10x/NovaSeq/cellranger_count_rawfullGTF/05.differential_gene_expression/heart_filtered/DESeq2_filterlowexpr/separate_clusters/", clusters_heart[x], "_", levels(cluster_metadata_heart$sex_id)[2], "_vs_", levels(cluster_metadata_heart$sex_id)[1], "_shrinkashr_all_genes.csv"), quote = FALSE, row.names = FALSE)
        
        str("Subset the significant results")
        sig_res_heart <- dplyr::filter(res_tbl_heart, padj < padj_cutoff) %>% dplyr::arrange(padj) 
        write.csv(sig_res_heart, paste0("~/Desktop/UBC/10x/NovaSeq/cellranger_count_rawfullGTF/05.differential_gene_expression/heart_filtered/DESeq2_filterlowexpr/separate_clusters/", clusters_heart[x], "_deseq_beforeshrink_sig_genes.csv"), quote = FALSE, row.names = FALSE)
        
        sig_res_heart_ashr <- dplyr::filter(res_tbl_heart_ashr, padj < padj_cutoff) %>% dplyr::arrange(padj)
        write.csv(sig_res_heart_ashr, paste0("~/Desktop/UBC/10x/NovaSeq/cellranger_count_rawfullGTF/05.differential_gene_expression/heart_filtered/DESeq2_filterlowexpr/separate_clusters/", clusters_heart[x], "_", levels(cluster_metadata_heart$sex_id)[A], "_vs_", levels(cluster_metadata_heart$sex_id)[B], "_shrinkashr_sig_genes.csv"), quote = FALSE, row.names = FALSE)
        
        sig_up_fc <- dplyr::filter(sig_res_heart, log2FoldChange >= log2fc_cutoff)
sig_down_fc <- dplyr::filter(sig_res_heart, log2FoldChange <= -log2fc_cutoff)
write.csv(sig_up_fc, paste0("~/Desktop/UBC/10x/NovaSeq/cellranger_count_rawfullGTF/05.differential_gene_expression/heart_filtered/DESeq2_filterlowexpr/separate_clusters/", clusters_heart[x], "_deseq_beforeshrink_sig_genes_fcpass_malebiased.csv"), quote = FALSE, row.names = FALSE)
write.csv(sig_down_fc, paste0("~/Desktop/UBC/10x/NovaSeq/cellranger_count_rawfullGTF/05.differential_gene_expression/heart_filtered/DESeq2_filterlowexpr/separate_clusters/", clusters_heart[x], "_deseq_beforeshrink_sig_genes_fcpass_femalebiased.csv"), quote = FALSE, row.names = FALSE)

sig_up_fc_ashr <- dplyr::filter(sig_res_heart_ashr, log2FoldChange >= log2fc_cutoff)
sig_down_fc_ashr <- dplyr::filter(sig_res_heart_ashr, log2FoldChange <= -log2fc_cutoff)
write.csv(sig_up_fc_ashr, paste0("~/Desktop/UBC/10x/NovaSeq/cellranger_count_rawfullGTF/05.differential_gene_expression/heart_filtered/DESeq2_filterlowexpr/separate_clusters/", clusters_heart[x], "_deseq_shrinkashr_sig_genes_fcpass_malebiased.csv"), quote = FALSE, row.names = FALSE)
write.csv(sig_down_fc_ashr, paste0("~/Desktop/UBC/10x/NovaSeq/cellranger_count_rawfullGTF/05.differential_gene_expression/heart_filtered/DESeq2_filterlowexpr/separate_clusters/", clusters_heart[x], "_deseq_shrinkashr_sig_genes_fcpass_femalebiased.csv"), quote = FALSE, row.names = FALSE)
        
        str("ggplot of top genes")
        normalized_counts_heart <- counts(dds_deseq_heart, normalized = TRUE)
        write.csv(normalized_counts_heart, paste0("~/Desktop/UBC/10x/NovaSeq/cellranger_count_rawfullGTF/05.differential_gene_expression/heart_filtered/DESeq2_filterlowexpr/separate_clusters/", clusters_heart[x], "_rnormalized_counts.csv"), quote = FALSE, row.names = T)
        
        
        str("Order results by padj values")
        top20_sig_genes_heart_ashr <- sig_res_heart_ashr %>%
                dplyr::arrange(padj) %>%
                dplyr::pull(gene) %>%
                head(n=20)
        
        top20_sig_norm_heart_ashr <- data.frame(normalized_counts_heart) %>%
                rownames_to_column(var = "gene") %>%
                dplyr::filter(gene %in% top20_sig_genes_heart_ashr)
        
        gathered_top20_sig_heart_ashr <- top20_sig_norm_heart_ashr %>%
                gather(colnames(top20_sig_norm_heart_ashr)[2:length(colnames(top20_sig_norm_heart_ashr))], key = "samplename", value = "normalized_counts")
        
        gathered_top20_sig_heart_ashr <- inner_join(ei_heart[, c("sample_id", "sex_id" )], gathered_top20_sig_heart_ashr, by = c("sample_id" = "samplename"))
        
        str("plot using ggplot2")    
         ggplot(gathered_top20_sig_heart_ashr) +
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
        ggsave(paste0("~/Desktop/UBC/10x/NovaSeq/cellranger_count_rawfullGTF/05.differential_gene_expression/heart_filtered/DESeq2_filterlowexpr/separate_clusters/plots/", clusters_heart[x], "_", levels(cluster_metadata_heart$sex_id)[A], "_vs_", levels(cluster_metadata_heart$sex_id)[B], "_shrinkashr_top20_DE_genes.png"))
		   
        str("volcano plot")
		res_tbl_thresh <- res_tbl_heart_ashr[!is.na(res_tbl_heart_ashr$padj), ] %>% mutate(threshold=padj<0.05 & abs(log2FoldChange) >= 1)
		min(log10(res_tbl_thresh$padj))
		
		ggplot(res_tbl_thresh) +
			geom_point(aes(x=log2FoldChange, y=-log10(padj), colour=threshold)) +
			ggtitle("Volcano plot") +
			xlab("log2 fold change") +
			ylab("-log10 adjusted p-value")+
			scale_color_manual(values=c("grey60", "red3")) +
			theme(legend.position="none", plot.title=element_text(size=rel(1.3), hjust=0.5), axis.title=element_text(size=rel(1.15)))
		ggsave(paste0("~/Desktop/UBC/10x/NovaSeq/cellranger_count_rawfullGTF/05.differential_gene_expression/heart_filtered/DESeq2_filterlowexpr/separate_clusters/plots/", clusters_heart[x], "_volcano_plot.png"))
		   

		
}      

map(1:length(clusters_heart), get_dds_resultsAvsB, A = 2, B = 1)






