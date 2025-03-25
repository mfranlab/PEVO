# ---
#   title: "NATCANCER_code"
# author: "Miguel Francisco"
# date: "24/03/2025"
# ---

#library prep
library(tidyverse)
library(ggplot2)
library(ComplexHeatmap)
require(RColorBrewer)

#data prep
test <- read.csv("pathto/oncoprint_baseline_mat_for_NatCancer.csv")
test_2 <- read.csv("pathto/oncoprint_pathways_baseline_mat_for_NatCancer.csv")
p_data_ordered <- read.csv("pathto/p_data_for_NatCancer.csv") 

## Figure 4
#colors
col_anno <- c("Anus" = "lemonchiffon1", 
              "Cervix" = "lightblue1",
              "missense" = "lightgoldenrod1",
              "Head/Neck" = "lightpink1",
              "Lung" = "lightsalmon1",
              "Penis" = "lightseagreen",
              "Vulva/Vagina" = "lightsalmon4",
              "NEGATIVE" = "gray0",
              "POSITIVE" = "tan1",
              "MSS" = "dodgerblue",
              "MSI" = "firebrick1",
              "low" = "forestgreen",
              "high" = "gold2",
              "-" = "gray50",
              "+" = "gray0",
              "Yes" = "red",
              "No" = "blue",
              "Objective response" = "lightblue",
              "No objective response" = "lightcoral",
              "Aflatoxin" = "lightgreen",
              "APOBEC" = "khaki",
              "Aristol ac" = "lightblue2",
              "Chewing tobacco" = "tan4",
              "Deamin 5MC" = "peru",
              "HRD" = "khaki4",
              "MMR" = "indianred2",
              "Platin" = "mediumpurple2",
              "POL" = "lightsalmon2",
              "ROS" = "palevioletred2",
              "Tobacco" = "tomato",
              "UV" = "royalblue1")

col_anno_2 <- c("No" = "red",
              "Yes" = "blue")

#sample split
mat_for_onco_baseline_ORR_res <- test[,colnames(test) %in% (subset(p_data_ordered, `CR_PR` == "Yes", select = sample_id_clinsight) %>% pull())]
mat_for_onco_baseline_ORR_nonres <- test[,colnames(test) %in% (subset(p_data_ordered, `CR_PR` == "No", select = sample_id_clinsight) %>% pull())]

#preparation for annotations
p_data_ordered_ORR_res <- p_data_ordered[p_data_ordered$sample_id_clinsight %in% (subset(p_data_ordered, `CR_PR` == "Yes", select = sample_id_clinsight) %>% pull()),]
p_data_ordered_ORR_nonres <- p_data_ordered[p_data_ordered$sample_id_clinsight %in% (subset(p_data_ordered, `CR_PR` == "No", select = sample_id_clinsight) %>% pull()),]

#column order
onco.id <- match(p_data_ordered_ORR_res$sample_id_clinsight, colnames(mat_for_onco_baseline_ORR_res))
onco.id
onco_colnames_ORR_res <- colnames(mat_for_onco_baseline_ORR_res)[onco.id]
mat_for_onco_baseline_ORR_res_test <- mat_for_onco_baseline_ORR_res[onco.id]
onco_colnames_ORR_res

onco.id <- match(p_data_ordered_ORR_nonres$sample_id_clinsight, colnames(mat_for_onco_baseline_ORR_nonres))
onco.id
onco_colnames_ORR_nonres <- colnames(mat_for_onco_baseline_ORR_nonres)[onco.id]
mat_for_onco_baseline_ORR_nonres_test <- mat_for_onco_baseline_ORR_nonres[onco.id]
onco_colnames_ORR_nonres

#annotations
p_data_test_anno_ORR_res <- HeatmapAnnotation(`Primary cancer site` = p_data_ordered_ORR_res$primary_cancer_site,
                                              TMB = p_data_ordered_ORR_res$`TMB (high and low)`,
                                              MSI = p_data_ordered_ORR_res$`MSI/MSS`,
                                              AS = p_data_ordered_ORR_res$AS_bi,
                                              HPV = p_data_ordered_ORR_res$`HPV (positive/negative)`,
                                              `CR/PR` = p_data_ordered_ORR_res$`CR_PR`,
                                              `PD-L1` = p_data_ordered_ORR_res$`PD-L1 IHC: positive/negative`,
                                              `Major mutational signature` = p_data_ordered_ORR_res$Short_names,
                                              col = list(`Primary cancer site` = col_anno,
                                                         HPV = col_anno,
                                                         MSI = col_anno,
                                                         AS = col_anno,
                                                         TMB = col_anno,
                                                         `PD-L1` = col_anno,
                                                         `CR/PR` = col_anno,
                                                         `Major mutational signature` = col_anno),
                                              show_annotation_name = TRUE)

p_data_test_anno_ORR_nonres <- HeatmapAnnotation(`Primary cancer site` = p_data_ordered_ORR_nonres$primary_cancer_site,
                                                 TMB = p_data_ordered_ORR_nonres$`TMB (high and low)`,
                                                 MSI = p_data_ordered_ORR_nonres$`MSI/MSS`,
                                                 AS = p_data_ordered_ORR_nonres$AS_bi,
                                                 HPV = p_data_ordered_ORR_nonres$`HPV (positive/negative)`,
                                                 `CR/PR` = p_data_ordered_ORR_nonres$CR_PR,
                                                 `PD-L1` = p_data_ordered_ORR_nonres$`PD-L1 IHC: positive/negative`,
                                                 `Major mutational signature` = p_data_ordered_ORR_nonres$Short_names,
                                                 col = list(`Primary cancer site` = col_anno,
                                                            HPV = col_anno,
                                                            MSI = col_anno,
                                                            AS = col_anno,
                                                            TMB = col_anno,
                                                            `PD-L1` = col_anno,
                                                            `CR/PR` = col_anno,
                                                            `Major mutational signature` = col_anno),
                                                 show_annotation_name = FALSE)

#oncoprint figure 4
oncoprint <- oncoPrint(mat_for_onco_baseline_ORR_nonres_test,
                       alter_fun = alter_fun, 
                       col = col.oncoprint,
                       show_column_names = FALSE,
                       remove_empty_columns = FALSE,
                       show_row_names = FALSE,
                       column_order = onco_colnames_ORR_nonres,
                       bottom_annotation = p_data_test_anno_ORR_nonres,
                       column_title = "Patients not experiencing an objective response",
                       row_split = oncogene[match(rownames(mat_for_onco_baseline_ORR_nonres_test), oncogene$Hugo.Symbol),"class"],
                       show_heatmap_legend = FALSE) +
  oncoPrint(mat_for_onco_baseline_ORR_res_test,
            alter_fun = alter_fun, 
            col = col.oncoprint,
            show_column_names = FALSE,
            show_row_names = TRUE,
            remove_empty_columns = FALSE,
            column_order = onco_colnames_ORR_res,
            bottom_annotation = p_data_test_anno_ORR_res,
            column_title = "Patients experiencing an objective response",
            row_split = oncogene[match(rownames(mat_for_onco_baseline_ORR_res_test), oncogene$Hugo.Symbol),"class"],
            show_heatmap_legend = TRUE)

draw(oncoprint)

## Figure 5
mat_for_onco_baseline_pathways_ORR_res <- mat_for_onco_baseline_pathways[,colnames(mat_for_onco_baseline_pathways) %in% (subset(p_data_ordered, `CR_PR` == "Yes", select = sample_id_clinsight) %>% pull())]
mat_for_onco_baseline_pathways_ORR_nonres <- mat_for_onco_baseline_pathways[,colnames(mat_for_onco_baseline_pathways) %in% (subset(p_data_ordered, `CR_PR` == "No", select = sample_id_clinsight) %>% pull())]

p_data_ordered_ORR_res <- p_data_ordered[p_data_ordered$sample_id_clinsight %in% (subset(p_data_ordered, `CR_PR` == "Yes", select = sample_id_clinsight) %>% pull()),]
p_data_ordered_ORR_nonres <- p_data_ordered[p_data_ordered$sample_id_clinsight %in% (subset(p_data_ordered, `CR_PR` == "No", select = sample_id_clinsight) %>% pull()),]

#for filtration
#to arrange columns by patient ID
colnames(mat_for_onco_baseline_pathways)
rownames(p_data_ordered)

p_data_ordered$`ORR_simplified`

onco.id <- match(p_data_ordered_ORR_res$sample_id_clinsight, colnames(mat_for_onco_baseline_pathways_ORR_res))
onco.id
onco_colnames_ORR_res <- colnames(mat_for_onco_baseline_pathways_ORR_res)[onco.id]
mat_for_onco_baseline_pathways_ORR_res_test <- mat_for_onco_baseline_pathways_ORR_res[onco.id]
onco_colnames_ORR_res

onco.id <- match(p_data_ordered_ORR_nonres$sample_id_clinsight, colnames(mat_for_onco_baseline_pathways_ORR_nonres))
onco.id
onco_colnames_ORR_nonres <- colnames(mat_for_onco_baseline_pathways_ORR_nonres)[onco.id]
mat_for_onco_baseline_pathways_ORR_nonres_test <- mat_for_onco_baseline_pathways_ORR_nonres[onco.id]
onco_colnames_ORR_nonres

p_data_test_anno <- HeatmapAnnotation(`Primary cancer site` = p_data_ordered_ORR$primary_cancer_site,
                                      TMB = p_data_ordered_ORR$`TMB (high and low)`,
                                      MSI = p_data_ordered_ORR$`MSI/MSS`,
                                      AS = p_data_ordered_ORR$AS_bi,
                                      HPV = p_data_ordered_ORR$`HPV (positive/negative)`,
                                      `CR/PR` = p_data_ordered_ORR$CR_PR,
                                      `PD-L1` = p_data_ordered_ORR$`PD-L1 IHC: positive/negative`,
                                      `Major mutational signature` = p_data_ordered_ORR$Short_names,
                                      col = list(`Primary cancer site` = col_anno,
                                                 HPV = col_anno,
                                                 MSI = col_anno,
                                                 AS = col_anno,
                                                 TMB = col_anno,
                                                 `PD-L1` = col_anno,
                                                 `CR/PR` = col_anno,
                                                 `Major mutational signature` = col_anno))

p_data_test_anno_ORR_res <- HeatmapAnnotation(`Primary cancer site` = p_data_ordered_ORR_res$primary_cancer_site,
                                              TMB = p_data_ordered_ORR_res$`TMB (high and low)`,
                                              MSI = p_data_ordered_ORR_res$`MSI/MSS`,
                                              AS = p_data_ordered_ORR_res$AS_bi,
                                              HPV = p_data_ordered_ORR_res$`HPV (positive/negative)`,
                                              `CR/PR` = p_data_ordered_ORR_res$CR_PR,
                                              `PD-L1` = p_data_ordered_ORR_res$`PD-L1 IHC: positive/negative`,
                                              `Major mutational signature` = p_data_ordered_ORR_res$Short_names,
                                              show_annotation_name = TRUE,
                                              col = list(`Primary cancer site` = col_anno,
                                                         HPV = col_anno,
                                                         MSI = col_anno,
                                                         AS = col_anno,
                                                         TMB = col_anno,
                                                         `PD-L1` = col_anno,
                                                         `CR/PR` = col_anno,
                                                         `Major mutational signature` = col_anno))

p_data_test_anno_ORR_nonres <- HeatmapAnnotation(`Primary cancer site` = p_data_ordered_ORR_nonres$primary_cancer_site,
                                                 TMB = p_data_ordered_ORR_nonres$`TMB (high and low)`,
                                                 MSI = p_data_ordered_ORR_nonres$`MSI/MSS`,
                                                 AS = p_data_ordered_ORR_nonres$AS_bi,
                                                 HPV = p_data_ordered_ORR_nonres$`HPV (positive/negative)`,
                                                 `CR/PR` = p_data_ordered_ORR_nonres$CR_PR,
                                                 `PD-L1` = p_data_ordered_ORR_nonres$`PD-L1 IHC: positive/negative`,
                                                 `Major mutational signature` = p_data_ordered_ORR_nonres$Short_names,
                                                 show_annotation_name = FALSE,
                                                 col = list(`Primary cancer site` = col_anno,
                                                            HPV = col_anno,
                                                            MSI = col_anno,
                                                            AS = col_anno,
                                                            TMB = col_anno,
                                                            `PD-L1` = col_anno,
                                                            `CR/PR` = col_anno,
                                                            `Major mutational signature` = col_anno))


oncoprint_pathway <- oncoPrint(mat_for_onco_baseline_pathways_ORR_nonres_test,
                               alter_fun = alter_fun, 
                               col = col.oncoprint,
                               show_column_names = FALSE,
                               remove_empty_columns = FALSE,
                               show_row_names = FALSE,
                               column_order = onco_colnames_ORR_nonres,
                               bottom_annotation = p_data_test_anno_ORR_nonres,
                               top_annotation = NULL,
                               column_title = "Patients not experiencing an objective response",
                               show_heatmap_legend = FALSE) +
  oncoPrint(mat_for_onco_baseline_pathways_ORR_res_test,
            alter_fun = alter_fun, 
            col = col.oncoprint,
            show_column_names = FALSE,
            show_row_names = TRUE,
            remove_empty_columns = FALSE,
            column_order = onco_colnames_ORR_res,
            bottom_annotation = p_data_test_anno_ORR_res,
            top_annotation = NULL,
            column_title = "Patients experiencing an objective response",
            show_heatmap_legend = FALSE)

oncoprint_pathway