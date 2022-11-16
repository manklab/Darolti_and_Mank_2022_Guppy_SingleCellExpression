
gonad_umap_rename_clusters_dbrmv_celltypes <- gonad_umap_3_rename_clusters
gonad_umap_rename_clusters_dbrmv_celltypes$CellType <- NA
gonad_umap_rename_clusters_dbrmv_celltypes$CellType[which(str_detect(gonad_umap_rename_clusters_dbrmv_celltypes$SCT_snn_res.0.5, "0"))] <- "Spermatocyte"
gonad_umap_rename_clusters_dbrmv_celltypes$CellType[which(str_detect(gonad_umap_rename_clusters_dbrmv_celltypes$SCT_snn_res.0.5, "1"))] <- "Spermatocyte"  
gonad_umap_rename_clusters_dbrmv_celltypes$CellType[which(str_detect(gonad_umap_rename_clusters_dbrmv_celltypes$SCT_snn_res.0.5, "2"))] <- "Spermatocyte"  
gonad_umap_rename_clusters_dbrmv_celltypes$CellType[which(str_detect(gonad_umap_rename_clusters_dbrmv_celltypes$SCT_snn_res.0.5, "3"))] <- "Spermatocyte"  
gonad_umap_rename_clusters_dbrmv_celltypes$CellType[which(str_detect(gonad_umap_rename_clusters_dbrmv_celltypes$SCT_snn_res.0.5, "4"))] <- "Granulosa"  
gonad_umap_rename_clusters_dbrmv_celltypes$CellType[which(str_detect(gonad_umap_rename_clusters_dbrmv_celltypes$SCT_snn_res.0.5, "5"))] <- "Spermatocyte"  
gonad_umap_rename_clusters_dbrmv_celltypes$CellType[which(str_detect(gonad_umap_rename_clusters_dbrmv_celltypes$SCT_snn_res.0.5, "6"))] <- "Macrophage"  
gonad_umap_rename_clusters_dbrmv_celltypes$CellType[which(str_detect(gonad_umap_rename_clusters_dbrmv_celltypes$SCT_snn_res.0.5, "7"))] <- "Spermatozoa"  
gonad_umap_rename_clusters_dbrmv_celltypes$CellType[which(str_detect(gonad_umap_rename_clusters_dbrmv_celltypes$SCT_snn_res.0.5, "8"))] <- "Endocrine"  
gonad_umap_rename_clusters_dbrmv_celltypes$CellType[which(str_detect(gonad_umap_rename_clusters_dbrmv_celltypes$SCT_snn_res.0.5, "9"))] <- "Spermatocyte"  
gonad_umap_rename_clusters_dbrmv_celltypes$CellType[which(str_detect(gonad_umap_rename_clusters_dbrmv_celltypes$SCT_snn_res.0.5, "10"))] <- "Germ"  
gonad_umap_rename_clusters_dbrmv_celltypes$CellType[which(str_detect(gonad_umap_rename_clusters_dbrmv_celltypes$SCT_snn_res.0.5, "11"))] <- "Spermatocyte"  
gonad_umap_rename_clusters_dbrmv_celltypes$CellType[which(str_detect(gonad_umap_rename_clusters_dbrmv_celltypes$SCT_snn_res.0.5, "12"))] <- "Spermatozoa"  
gonad_umap_rename_clusters_dbrmv_celltypes$CellType[which(str_detect(gonad_umap_rename_clusters_dbrmv_celltypes$SCT_snn_res.0.5, "13"))] <- "Spermatocyte"  
gonad_umap_rename_clusters_dbrmv_celltypes$CellType[which(str_detect(gonad_umap_rename_clusters_dbrmv_celltypes$SCT_snn_res.0.5, "14"))] <- "Mitochondrial-rich"  
gonad_umap_rename_clusters_dbrmv_celltypes$CellType[which(str_detect(gonad_umap_rename_clusters_dbrmv_celltypes$SCT_snn_res.0.5, "15"))] <- "Spermatocyte"  
gonad_umap_rename_clusters_dbrmv_celltypes$CellType[which(str_detect(gonad_umap_rename_clusters_dbrmv_celltypes$SCT_snn_res.0.5, "16"))] <- "Spermatocyte"  
gonad_umap_rename_clusters_dbrmv_celltypes$CellType[which(str_detect(gonad_umap_rename_clusters_dbrmv_celltypes$SCT_snn_res.0.5, "17"))] <- "Spermatocyte"  
gonad_umap_rename_clusters_dbrmv_celltypes$CellType[which(str_detect(gonad_umap_rename_clusters_dbrmv_celltypes$SCT_snn_res.0.5, "18"))] <- "Spermatocyte"  
gonad_umap_rename_clusters_dbrmv_celltypes$CellType[which(str_detect(gonad_umap_rename_clusters_dbrmv_celltypes$SCT_snn_res.0.5, "19"))] <- "Endocrine"  
gonad_umap_rename_clusters_dbrmv_celltypes$CellType[which(str_detect(gonad_umap_rename_clusters_dbrmv_celltypes$SCT_snn_res.0.5, "20"))] <- "Spermatocyte"  
gonad_umap_rename_clusters_dbrmv_celltypes$CellType[which(str_detect(gonad_umap_rename_clusters_dbrmv_celltypes$SCT_snn_res.0.5, "21"))] <- "Embryonic Yolk Sac"  
gonad_umap_rename_clusters_dbrmv_celltypes$CellType[which(str_detect(gonad_umap_rename_clusters_dbrmv_celltypes$SCT_snn_res.0.5, "22"))] <- "Spermatocyte"  
gonad_umap_rename_clusters_dbrmv_celltypes$CellType[which(str_detect(gonad_umap_rename_clusters_dbrmv_celltypes$SCT_snn_res.0.5, "23"))] <- "Spermatocyte"  
gonad_umap_rename_clusters_dbrmv_celltypes$CellType[which(str_detect(gonad_umap_rename_clusters_dbrmv_celltypes$SCT_snn_res.0.5, "24"))] <- "Spermatozoa"  
gonad_umap_rename_clusters_dbrmv_celltypes$CellType[which(str_detect(gonad_umap_rename_clusters_dbrmv_celltypes$SCT_snn_res.0.5, "25"))] <- "Spermatocyte"  
gonad_umap_rename_clusters_dbrmv_celltypes$CellType[which(str_detect(gonad_umap_rename_clusters_dbrmv_celltypes$SCT_snn_res.0.5, "26"))] <- "Spermatocyte"  
gonad_umap_rename_clusters_dbrmv_celltypes$CellType[which(str_detect(gonad_umap_rename_clusters_dbrmv_celltypes$SCT_snn_res.0.5, "27"))] <- "Spermatocyte"  
gonad_umap_rename_clusters_dbrmv_celltypes$CellType[which(str_detect(gonad_umap_rename_clusters_dbrmv_celltypes$SCT_snn_res.0.5, "28"))] <- "Spermatocyte"  
gonad_umap_rename_clusters_dbrmv_celltypes$CellType[which(str_detect(gonad_umap_rename_clusters_dbrmv_celltypes$SCT_snn_res.0.5, "29"))] <- "Spermatocyte"
gonad_umap_rename_clusters_dbrmv_celltypes$CellType[which(str_detect(gonad_umap_rename_clusters_dbrmv_celltypes$SCT_snn_res.0.5, "30"))] <- "Spermatocyte"  
gonad_umap_rename_clusters_dbrmv_celltypes$CellType[which(str_detect(gonad_umap_rename_clusters_dbrmv_celltypes$SCT_snn_res.0.5, "31"))] <- "Spermatocyte"    
  
table(gonad_umap_rename_clusters@active.ident)
      Spermatocyte        Spermatozoa          Granulosa         Macrophage          Endocrine               Germ Mitochondrial-rich Embryonic Yolk Sac 
             17906               2754               1583               1478               1852               1319                762                360 

table(gonad_umap_rename_clusters_dbrmv_celltypes$sample)
GonadF1 GonadF2 GonadF3 GonadM4 GonadM5 GonadM6 
   3162    2063    1492    7973    6541    6783


## Proportions test (prop.test) or z-test method (similar to chisq.test) - Females:Males

tab_gonad_endocr <- matrix(c(1389, 463, 5328, 20834), ncol=2)
prop.test(tab_gonad_endocr)
data:  tab_gonad_endocr
X-squared = 2829.2, df = 1, p-value < 2.2e-16
alternative hypothesis: two.sided
95 percent confidence interval:
 0.1750692 0.1950280
sample estimates:
    prop 1     prop 2 
0.20678874 0.02174015 

tab_gonad_germ <- matrix(c(1319, 0, 5398, 21297), ncol=2)
prop.test(tab_gonad_germ)
data:  tab_gonad_germ
X-squared = 4384.3, df = 1, p-value < 2.2e-16
alternative hypothesis: two.sided
95 percent confidence interval:
 0.1867695 0.2059653
sample estimates:
   prop 1    prop 2 
0.1963674 0.0000000 

tab_gonad_granul <- matrix(c(1469, 114, 5248, 21183), ncol=2)
prop.test(tab_gonad_granul)
data:  tab_gonad_granul
X-squared = 4355.6, df = 1, p-value < 2.2e-16
alternative hypothesis: two.sided
95 percent confidence interval:
 0.2033142 0.2233777
sample estimates:
     prop 1      prop 2 
0.218698824 0.005352867 

tab_gonad_embryo <- matrix(c(271, 89, 6446, 21208), ncol=2)
prop.test(tab_gonad_embryo)
data:  tab_gonad_embryo
X-squared = 523.68, df = 1, p-value < 2.2e-16
alternative hypothesis: two.sided
95 percent confidence interval:
 0.03128379 0.04104901
sample estimates:
     prop 1      prop 2 
0.040345392 0.004178992 

tab_gonad_macro <- matrix(c(1409, 69, 5308, 21228), ncol=2)
prop.test(tab_gonad_macro)
data:  tab_gonad_macro
X-squared = 4354.1, df = 1, p-value < 2.2e-16
alternative hypothesis: two.sided
95 percent confidence interval:
 0.1966620 0.2163907
sample estimates:
     prop 1      prop 2 
0.209766265 0.003239893 

tab_gonad_mito <- matrix(c(745, 17, 5972, 21280), ncol=2)
prop.test(tab_gonad_mito)
data:  tab_gonad_mito
X-squared = 2335.8, df = 1, p-value < 2.2e-16
alternative hypothesis: two.sided
95 percent confidence interval:
 0.1024972 0.1177316
sample estimates:
      prop 1       prop 2 
0.1109126098 0.0007982345 

tab_gonad_spermatocyte <- matrix(c(111, 17795, 6606, 3502), ncol=2)
prop.test(tab_gonad_spermatocyte)
data:  tab_gonad_spermatocyte
X-squared = 14849, df = 1, p-value < 2.2e-16
alternative hypothesis: two.sided
95 percent confidence interval:
 -0.8249740 -0.8131029
sample estimates:
    prop 1     prop 2 
0.01652523 0.83556369 
		
tab_gonad_spermatozoa <- matrix(c(4, 2750, 6713, 18547), ncol=2)
prop.test(tab_gonad_spermatozoa)	
data:  tab_gonad_spermatozoa
X-squared = 950.22, df = 1, p-value < 2.2e-16
alternative hypothesis: two.sided
95 percent confidence interval:
 -0.1331700 -0.1238914
sample estimates:
      prop 1       prop 2 
0.0005955039 0.1291261680 

