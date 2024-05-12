if(!require(devtools)) install.packages("devtools")
devtools::install_github("kassambara/survminer")
library("survminer")
require("survival")
library(cytoreason.validator.apps.client)
library(cytoreason.ccm.pipeline)

#######################################################################3

#survival for CCM dataset
CCM_run <- as_ccm_fit("wf-e64142ff5d") # !!!!!!!! #
# choose your ds
ds <- CCM_run$datasets$GSE30219 # !!!!!!!! #
dataset_values <- get_eset_from_dataset("GSE30219_eRh4m",platformID(ds))  # !!!!!!!! #
dataset_ann <- cytoreason.curator.annotations::get_annotations(experimentID(ds))
colnames(dataset_values@phenoData@data)

time <- as.numeric(dataset_values@phenoData@data[["follow-up time (months):ch1"]])
dataset_values$overallsurvival_os_num <- dataset_values@phenoData@data[["Status"]]
dataset_values$overallsurvival_os_num[dataset_values$overallsurvival_os_num=="ALIVE"] <- 0
dataset_values$overallsurvival_os_num[dataset_values$overallsurvival_os_num=="DEAD"] <- 1

# conect the data from TME table:
data_TME <- data.frame(sample_id = ds$tme_classifier$tme_predict$sample_id, cluster = ds$tme_classifier$tme_predict$predict)
data_TME <- data_TME[!duplicated(data_TME), ]
dataset_values$TME <- ""
known_index <- dataset_values$sample_id %in% data_TME$sample_id
dataset_values$TME[known_index] <- as.character(data_TME$cluster)

# plot regular
os_time <- time[known_index]
os_data <- as.numeric((dataset_values$overallsurvival_os_num)[known_index])
os_TME <- dataset_values$TME[known_index]
fit_gdf <- survfit(Surv(os_time, os_data) ~ os_TME, data = dataset_values) 

ggsurvplot(fit_gdf,  legend = "right", legend.title = "TMEs", 
           pval = TRUE, pval.method=T)+ ggtitle(paste0(experimentID(ds)," TMEs"))


# check possible annotations:
groups_data <- curator_sample_annotation(experiment_id = experimentID(ds))


# plot LUAD only
LUAD_index <- groups_data$group__LUAD_vs_LUSC == "B"
LUAD_time <- time[LUAD_index]
LUAD_data <- as.numeric((dataset_values$overallsurvival_os_num)[LUAD_index])
LUAD_TME <- dataset_values$TME[LUAD_index]
fit_gdf <- survfit(Surv(LUAD_time, LUAD_data) ~ LUAD_TME, data = dataset_values) 
ggsurvplot(fit_gdf,  legend = "right", legend.title = "TMEs", 
           pval = TRUE, pval.method=T)+ ggtitle(paste0(experimentID(ds)," TMEs", "\n", "N=", length(LUAD_TME)))



# plot LUSC only
LUSC_index <- groups_data$group__LUAD_vs_LUSC == "A"
LUSC_time <- time[LUSC_index]
LUSC_data <- as.numeric((dataset_values$overallsurvival_os_num)[LUSC_index])
LUSC_TME <- dataset_values$TME[LUSC_index]
fit_gdf <- survfit(Surv(LUSC_time, LUSC_data) ~ LUSC_TME, data = dataset_values) 
ggsurvplot(fit_gdf,  legend = "right", legend.title = "TMEs", 
           pval = TRUE, pval.method=T)+ ggtitle(paste0(experimentID(ds)," TMEs","\n", "N=", length(LUSC_TME)))
