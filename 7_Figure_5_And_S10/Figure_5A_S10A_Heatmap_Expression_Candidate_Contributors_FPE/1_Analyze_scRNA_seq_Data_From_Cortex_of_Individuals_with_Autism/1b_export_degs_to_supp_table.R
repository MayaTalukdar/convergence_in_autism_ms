###################
#I/O
###################
.libPaths(c("/n/groups/walsh/indData/Maya/RLibs", .libPaths()))
library(tidyverse)
library(openxlsx)

#read in degs
idiopathic_exn_degs <- readRDS("/n/groups/walsh/indData/Maya/Finalized_Convergence_In_ASD/Cell_Lines/figures/analysis_of_mo_dias_15q_data/Output/DEG_Results/MAST/Case_vs_Control_Male/ASD/Neurons_DEG_res.RDS")
cnv_exn_degs <- readRDS("/n/groups/walsh/indData/Maya/Finalized_Convergence_In_ASD/Cell_Lines/figures/analysis_of_mo_dias_15q_data/Output/DEG_Results/MAST/Case_vs_Control_Male/CNV/Neurons_DEG_res.RDS")
case_exn_degs <- readRDS("/n/groups/walsh/indData/Maya/Finalized_Convergence_In_ASD/Cell_Lines/figures/analysis_of_mo_dias_15q_data/Output/DEG_Results/MAST/Case_vs_Control_Male/case/Neurons_DEG_res.RDS")

write.xlsx(
  list(
    Idiopathic = idiopathic_exn_degs,
    CNV_15q = cnv_exn_degs,
    Case = case_exn_degs
  ),
  file = "/n/groups/walsh/indData/Maya/Finalized_Convergence_In_ASD/Cell_Lines/figures/analysis_of_mo_dias_15q_data/Output/Supp_Tables/Mo_Dias_DEG_results.xlsx",
  rowNames = TRUE,
  overwrite = TRUE
)