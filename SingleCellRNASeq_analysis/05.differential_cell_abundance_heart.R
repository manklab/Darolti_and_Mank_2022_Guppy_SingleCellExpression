
heart_umap_rename_clusters_dbrmv_celltypes <- heart_umap_rename_clusters_dbrmv
heart_umap_rename_clusters_dbrmv_celltypes$CellType <- NA
heart_umap_rename_clusters_dbrmv_celltypes$CellType[which(str_detect(heart_umap_rename_clusters_dbrmv_celltypes$SCT_snn_res.0.2, "0"))] <- "Fibroblast"  
heart_umap_rename_clusters_dbrmv_celltypes$CellType[which(str_detect(heart_umap_rename_clusters_dbrmv_celltypes$SCT_snn_res.0.2, "1"))] <- "Macrophage" 
heart_umap_rename_clusters_dbrmv_celltypes$CellType[which(str_detect(heart_umap_rename_clusters_dbrmv_celltypes$SCT_snn_res.0.2, "2"))] <- "Cardiomyocyte 1" 
heart_umap_rename_clusters_dbrmv_celltypes$CellType[which(str_detect(heart_umap_rename_clusters_dbrmv_celltypes$SCT_snn_res.0.2, "3"))] <- "Erythrocyte" 
heart_umap_rename_clusters_dbrmv_celltypes$CellType[which(str_detect(heart_umap_rename_clusters_dbrmv_celltypes$SCT_snn_res.0.2, "4"))] <- "T lymphocyte" 
heart_umap_rename_clusters_dbrmv_celltypes$CellType[which(str_detect(heart_umap_rename_clusters_dbrmv_celltypes$SCT_snn_res.0.2, "5"))] <- "B lymphocyte" 
heart_umap_rename_clusters_dbrmv_celltypes$CellType[which(str_detect(heart_umap_rename_clusters_dbrmv_celltypes$SCT_snn_res.0.2, "6"))] <- "Erythrocyte" 
heart_umap_rename_clusters_dbrmv_celltypes$CellType[which(str_detect(heart_umap_rename_clusters_dbrmv_celltypes$SCT_snn_res.0.2, "7"))] <- "Endocardial" 
heart_umap_rename_clusters_dbrmv_celltypes$CellType[which(str_detect(heart_umap_rename_clusters_dbrmv_celltypes$SCT_snn_res.0.2, "8"))] <- "Vascular endothelial" 
heart_umap_rename_clusters_dbrmv_celltypes$CellType[which(str_detect(heart_umap_rename_clusters_dbrmv_celltypes$SCT_snn_res.0.2, "9"))] <- "Cardiomyocyte 2" 
heart_umap_rename_clusters_dbrmv_celltypes$CellType[which(str_detect(heart_umap_rename_clusters_dbrmv_celltypes$SCT_snn_res.0.2, "10"))] <- "Granulocyte" 
heart_umap_rename_clusters_dbrmv_celltypes$CellType[which(str_detect(heart_umap_rename_clusters_dbrmv_celltypes$SCT_snn_res.0.2, "11"))] <- "Vascular smooth muscle" 

table(heart_umap_rename_clusters_dbrmv_celltypes$CellType)
          B lymphocyte        Cardiomyocyte 1        Cardiomyocyte 2            Endocardial            Erythrocyte             Fibroblast            Granulocyte             Macrophage           T lymphocyte 
                  1497                   4531                    716                   1153                   5640                   6583                    284                   4107                   3407 
  Vascular endothelial Vascular smooth muscle 
                   721                    284 

table(heart_umap_rename_clusters_dbrmv_celltypes$sample)
HeartF1 HeartF2 HeartF3 HeartM1 HeartM2 HeartM3 
   6735    8966    6150    4615    1163    1294 


## Proportions test (prop.test) or z-test method (similar to chisq.test) - Females:Males

tab_heart_cardio1 <- matrix(c(3541, 990, 18310, 6082), ncol=2)
prop.test(tab_heart_cardio1)
data:  tab_heart_cardio1
X-squared = 19.519, df = 1, p-value = 9.958e-06
alternative hypothesis: two.sided
95 percent confidence interval:
 0.01252161 0.03160517
sample estimates:
   prop 1    prop 2 
0.1620521 0.1399887 

tab_heart_cardio2 <- matrix(c(352, 364, 21499, 6708), ncol=2)
prop.test(tab_heart_cardio2)
data:  tab_heart_cardio2
X-squared = 275.26, df = 1, p-value < 2.2e-16
alternative hypothesis: two.sided
95 percent confidence interval:
 -0.04086855 -0.02985442
sample estimates:
    prop 1     prop 2 
0.01610910 0.05147059 

tab_heart_endocard <- matrix(c(716, 437, 21135, 6635), ncol=2)
prop.test(tab_heart_endocard)
data:  tab_heart_endocard
X-squared = 116.84, df = 1, p-value < 2.2e-16
alternative hypothesis: two.sided
95 percent confidence interval:
 -0.03520715 -0.02284407
sample estimates:
    prop 1     prop 2 
0.03276738 0.06179299 

tab_heart_erythro <- matrix(c(3910, 1730, 17941, 5342), ncol=2)
prop.test(tab_heart_erythro)
data:  tab_heart_erythro
X-squared = 146.44, df = 1, p-value < 2.2e-16
alternative hypothesis: two.sided
95 percent confidence interval:
 -0.07701509 -0.05435995
sample estimates:
   prop 1    prop 2 
0.1789392 0.2446267 

tab_heart_fibro <- matrix(c(5319, 1264, 16532, 5808), ncol=2)
prop.test(tab_heart_fibro)
data:  tab_heart_fibro
X-squared = 126.81, df = 1, p-value < 2.2e-16
alternative hypothesis: two.sided
95 percent confidence interval:
 0.05400649 0.07537015
sample estimates:
   prop 1    prop 2 
0.2434214 0.1787330 

tab_heart_granulo <- matrix(c(204, 80, 21647, 6992), ncol=2)
prop.test(tab_heart_granulo)
data:  tab_heart_granulo
X-squared = 1.9477, df = 1, p-value = 0.1628
alternative hypothesis: two.sided
95 percent confidence interval:
 -0.0048449408  0.0008924208
sample estimates:
     prop 1      prop 2 
0.009335957 0.011312217 

tab_heart_macro <- matrix(c(3142, 965, 18709, 6107), ncol=2)
prop.test(tab_heart_macro)
data:  tab_heart_macro
X-squared = 2.3018, df = 1, p-value = 0.1292
alternative hypothesis: two.sided
95 percent confidence interval:
 -0.002009921  0.016686774
sample estimates:
   prop 1    prop 2 
0.1437920 0.1364536 

tab_heart_blymph <- matrix(c(1303, 194, 20548, 6878), ncol=2)
prop.test(tab_heart_blymph)
data:  tab_heart_blymph
X-squared = 112.21, df = 1, p-value < 2.2e-16
alternative hypothesis: two.sided
95 percent confidence interval:
 0.02717082 0.03722720
sample estimates:
    prop 1     prop 2 
0.05963114 0.02743213 

tab_heart_vascendo <- matrix(c(444, 277, 21407, 6795), ncol=2)
prop.test(tab_heart_vascendo)
data:  tab_heart_vascendo
X-squared = 77.322, df = 1, p-value < 2.2e-16
alternative hypothesis: two.sided
95 percent confidence interval:
 -0.02383579 -0.01386244
sample estimates:
    prop 1     prop 2 
0.02031944 0.03916855 

tab_heart_vascmusc <- matrix(c(214, 70, 21637, 7002), ncol=2)
prop.test(tab_heart_vascmusc)
data:  tab_heart_vascmusc
X-squared = 6.6544e-05, df = 1, p-value = 0.9935
alternative hypothesis: two.sided
95 percent confidence interval:
 -0.002849262  0.002640087
sample estimates:
     prop 1      prop 2 
0.009793602 0.009898190 
		
tab_heart_tlym <- matrix(c(2706, 701, 19145, 6371), ncol=2)
prop.test(tab_heart_tlym)	
data:  tab_heart_tlym
X-squared = 31.168, df = 1, p-value = 2.366e-08
alternative hypothesis: two.sided
95 percent confidence interval:
 0.01640107 0.03302977
sample estimates:
   prop 1    prop 2 
0.1238387 0.0991233
