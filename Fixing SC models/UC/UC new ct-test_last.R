### env ###

# install all packages 
library(unixtools)
unixtools::set.tempdir("/scratch/tmp")

install.packages('cytoreason.cc.client', repos = c(CRRAN = "https://crran.cytoreason.com/R/latest"))
library(cytoreason.cc.client)
library(SingleCellExperiment)
library(cytoreason.io)
library(cytoreason.ccm.pipeline)
library('cytoreason')
library(pkgmaker)
library(cytoreason.datasource)
# library(cytoreason.single.cell)
library(plyr)

# load_all("~/capsule/code/service-single-cell-preprocessing/")
load_all("~/capsule/code/service-single-cell/")


### ~ define output directory and image ~ ###
out_path <- "~/capsule/code/service-single-cell/notes/output/"
IMAGE <- "eu.gcr.io/cytoreason/ci-cytoreason.single.cell-package:MF_32_SC_gene_diff_wrapper_for_CCM_output_latest"

### update the last UC object ##################################
# get the last object 
uc_all = readRDS(get_task_outputs("wf-f6843f0de9", "0"))

# change the samples to avoid technical replicates
# healthy subjects have 2 EPI/LP each (without biological differences between them)
# N52, N58 and N111 have technical replicates as well
uc_samples = uc_all$biosample_id

# healthy: N8.LPA, N8.LPB -> N8.LP
# UC: "N52.LPB1a" "N52.LPB1b" "N52.LPB2a" "N52.LPB2b" -> N52.LPB
clean_reps <- function(samples, condition) {
  cleaned_samples <- ifelse(condition == "Healthy",
                           gsub("([AB]).*", "", samples),
                           gsub("([AB]).*", "\\1", samples))
  return(cleaned_samples)
}

new_sample = clean_reps(uc_all$biosample_id, uc_all$Type)

# comment: in the upload process all the dots chnage to underscores
# we will chanage it in the seurat object (seems to be optional)

# Define a function to replace dot with underscore and add underscore between "Epi"/"LP" and "A"/"B"
replace_characters <- function(string) {
  # Replace "." with "_"
  string <- gsub("\\.", "_", string)
  # Add "_" between "Epi"/"LP" and "A"/"B"
  string <- gsub("(Epi|LP)([A-Za-z])", "\\1_\\2", string)
  return(string)
}

# Apply the function to each string in the vector
new_sample_f <- sapply(new_sample, replace_characters)

# the function above change the sample ID from "N7.EpiA" to "N7_Epi_A"
uc_all$new_sample = new_sample_f

# adding subject ID to the seurat object MD 
subjectid = sub("\\..*", "", uc_all$biosample_id)
uc_all$"Subject ID" = subjectid  

save_to_cyto_cc(uc_all, tag = "uc_all_with_corrected_samples_and_subjects")
# Cyto-CC workflow: wf-3a4a083971
# Cyto-CC workflow: wf-b31b63f9c0 - with subject id
# Cyto-CC workflow: wf-084b8bf1a1

sc_fit <- as_sc_pipeline_run(data = "wf-084b8bf1a1")
seurat = as.Seurat(sc_fit, majority_vote = F)
wf_export <- sc_api_save_data(dataset_name = "p00_sc_uc_colon",
                              data = sc_fit,
                              tags = list(list(name = "stage", value = "upload_to_BQ"),
                                          list(name = "dataset", value = "uc_colon"),
                                          list(name="comment", value = "uc_all_with_corrected_samples")))

# [Tue Feb 20 09:01:12 2024] save_data_sc_pipeline - wf-96c3519949 -> failed (when setting directly the object)
# [Tue Feb 20 15:23:12 2024] save_data_sc_pipeline - wf-771786781e 
# [Sun Feb 25 09:36:55 2024] save_data_sc_pipeline - wf-0ea001cb25 -> with subject id
# [Sun Feb 25 14:19:13 2024] save_data_sc_pipeline - wf-86610a80e3



#### update celltypes names 140324 #######################################
# get the last seurat object
seurat = readRDS(get_task_outputs("wf-084b8bf1a1", task_id = 0))

celltypes = unique(seurat$biologist_annotation)

# get the last ontology 
colon.o = read.csv("~/capsule/scratch/ontology_colon_140324.csv")
seurat$biologist_annotation[seurat$biologist_annotation == "pericyte cell"] <- "pericyte"

library(dplyr)
colon.o.tmp = colon.o[colon.o$label %in% unique(seurat$biologist_annotation),]
updated_b_annotaions = mapvalues(seurat$biologist_annotation,
                                 from=colon.o.tmp$label,
                                 to=colon.o.tmp$Long.display.name)
seurat$biologist_annotation = updated_b_annotaions

save_to_cyto_cc(seurat, tag = "uc_without_technical_rep_and_updated_cell_names")
# Cyto-CC workflow: wf-723640b717

sc_fit <- as_sc_pipeline_run(data = "wf-723640b717")
seurat = as.Seurat(sc_fit, majority_vote = F)
wf_export <- sc_api_save_data(dataset_name = "p00_sc_uc_colon",
                              data = sc_fit,
                              dataset_id = "uc_smillie",
                              image = "eu.gcr.io/cytoreason/ci-cytoreason.single.cell.preprocessing-package:master_latest",
                              tags = list(list(name = "stage", value = "upload_to_BQ"),
                                          list(name = "dataset", value = "uc_colon"),
                                          list(name="comment", value = "uc_all_without_rep_correct_names")))
# [Sat Mar 16 16:12:59 2024] save_data_sc_pipeline - wf-47af6b80db
  
  
#### update the metadata #################################################
# get the last metadata
smillie.csv = read.csv("capsule/code/sliced_Smillie study Annotations - SC for UC.csv", row.names = 1)

# update the metadata to the new samples 
adjust_subject_ids <- function(subject_id, condition) {
  adjusted_ids <- character(length(subject_id))
  
  for (i in seq_along(subject_id)) {
    if (condition[i] == "healthy") {
      # Remove the trailing "_A" or "_B" for healthy condition
      adjusted_ids[i] <- sub("_(A|B)$", "", subject_id[i])
    } else if (condition[i] == "ulcerative colitis" && grepl("_(A|B)[1-9]", subject_id[i])) {
      # Remove only the values after "_A" or "_B" for disease condition, if characters exist after "_A" or "_B"
      adjusted_ids[i] <- sub("(A|B).*", "\\1", subject_id[i])
    } else {
      # For unknown conditions or if no characters exist after "_A" or "_B", keep the original subject ID
      adjusted_ids[i] <- subject_id[i]
    }
  }
  
  return(adjusted_ids)
}

new_sample_id = adjust_subject_ids(smillie.csv$sample_id, smillie.csv$condition)
smillie.csv$sample_id = new_sample_id
write.csv(smillie.csv, "~/capsule/code/smillie_study_anntation_corrected.csv")
smillie.csv = read.csv("~/capsule/code/smillie_study_anntation_corrected.csv", row.names = 1)

# smillie.csv["N51_Epi_A","tissue_comment"] = "epithelium"
# smillie.csv["N51_LP_A","tissue_comment"] = "lamina propria"
rownames(smillie.csv) = smillie.csv$sample_id


## rerun ct-test and csde per layer -----
config_test <- sc_api_configuration(
  asset_id = "wf-47af6b80db",
  cell_annotation = "biologist_annotation",
  experiment_id = "uc_colon_smillie",
  platform_id = "sc-rnaseq",
  biosample_id = "new_sample",
  sample_annotations = read.csv("~/capsule/code/smillie_study_anntation_corrected.csv", row.names = 1),
  term_metadata = list(UC = list(DZ_vs_HC = list(contrast_type = "disease_vs_control")),
                       UC_E = list(DZ.Epi_vs_HC.Epi = list(contrast_type = "disease_vs_control:epithelial_cells")),
                       UC_L = list(DZ.LP_vs_HC.LP = list(contrast_type = "disease_vs_control:lamina_propria_cells")),
                       
                       I_vs_NI = list(I_vs_NI = list(contrast_type = "inflamed_vs_non_inflamed")),
                       I_vs_NI_E = list(I.Epi_vs_NI.Epi = list(contrast_type = "inflamed_vs_non_inflamed:epithelial_cells")),
                       I_vs_NI_L = list(I.LP_vs_NI.LP = list(contrast_type = "inflamed_vs_non_inflamed:lamina_propria_cells")),
                       
                       UC_inf = list(DZ_inf_vs_HC = list(contrast_type = "inflamed_vs_healthy")),
                       UC_inf_E = list(DZ_inf.Epi_vs_HC.Epi = list(contrast_type = "inflamed_vs_healthy:epithelial_cells")),
                       UC_inf_L = list(DZ_inf.LP_vs_HC.LP = list(contrast_type = "inflamed_vs_healthy:lamina_propria_cells"))),
  
  analysis_model = list(UC = list(group = list(condition = c(HC = "healthy", DZ = "ulcerative colitis")), covariates = c("sub_tissue")),
                        UC_E = list(group = list(condition = c(HC = "healthy", DZ = "ulcerative colitis"), sub_tissue = "Epi")),
                        UC_L = list(group = list(condition = c(HC = "healthy", DZ = "ulcerative colitis"), sub_tissue = "LP")),
                        
                        
                        I_vs_NI = list(group = list(sample_classification = c(NI = "Non inflamed", I = "Inflamed")), covariates = c("sub_tissue")),
                        I_vs_NI_E = list(group = list(sample_classification = c(NI = "Non inflamed", I = "Inflamed"), sub_tissue = "Epi")),
                        I_vs_NI_L = list(group = list(sample_classification = c(NI = "Non inflamed", I = "Inflamed"), sub_tissue = "LP")),
                        
                        
                        UC_inf = list(group = list(sample_classification = c(HC = "Normal", DZ_inf = "Inflamed")), covariates = c("sub_tissue")),
                        UC_inf_E = list(group = list(sample_classification = c(HC = "Normal", DZ_inf = "Inflamed"), sub_tissue = "Epi")),
                        UC_inf_L = list(group = list(sample_classification = c(HC = "Normal", DZ_inf = "Inflamed"), sub_tissue = "LP"))
  ),
  parameters = list(model=list(services = c("ct_test", "gx_diff"), ct_test=list(), gx_diff=list(min_sample_size=15)))
)


#### check the configuration locally ----
sc_validate_dataset_configuration(config_test)

#### check that a valid dataset can be loaded from the config ----
pb_object <- prepare_single_cell_dataset(config_test)

# check a particular analysis model (if failing)
obj <- prepare_single_cell_dataset(config_test, validate = FALSE)
str(build_analysis_model_metadata(config_test$analysis_model$UC, obj))

sc_api_run_analysis("UC-CO",
                    config_test,
                    tags = list(message = "UC with colon tissue only no replication layer",
                                test_run = 4)) # ,image = IMAGE
# [Sun Mar 10 16:08:21 2024] run_analysis - wf-339724564f -> L instead of LP 
# [Sun Mar 10 16:13:40 2024] run_analysis - wf-b6ea0962a2
# [Tue Mar 19 12:18:41 2024] run_analysis - wf-7a96dd7824
# [Tue Mar 19 12:27:24 2024] run_analysis - wf-16add6d79a (force_execution = TRUE) failed
# [Tue Mar 19 16:20:36 2024] run_analysis - wf-094eb8ed56 failed
# [Tue Mar 19 16:41:32 2024] run_analysis - wf-f81aa8bced fully cashed 
# [Wed Mar 20 15:17:14 2024] run_analysis - wf-006a92c6ed -> correction of i ni -> FAILED DUE TO MISTAKE 
# [Wed Mar 20 16:07:56 2024] run_analysis - wf-ab6c641c10
# [Thu Mar 21 09:34:00 2024] run_analysis - wf-aaeda51e5f -> correct seurat object- updated names 

command <- paste("python /app/cytobigquery/exec_service.py single_cell", wf_id = "wf-f81aa8bced",
                 dataset_name = "p00_ulcerative_colitis", "--verbose -sc -ti=eu.gcr.io/cytoreason/cd-py-bigquery:develop_latest -tm=124Mi")
res <-
  run_command_dist(
    command,
    image = "eu.gcr.io/cytoreason/cd-py-bigquery:develop_latest",
    memory_request = "25Gi",
    force_execution = FALSE,
    replace_image_tags = TRUE
  )
# 110324 1033 Cyto-CC workflow: wf-ee11f5b78d
# 170324 Cyto-CC workflow: wf-858ef77aee
# 190324 1357 Cyto-CC workflow: wf-00ba5e9a66
# 190324 1720 Cyto-CC workflow: wf-179f4299f9
# 190324 1848 Cyto-CC workflow: wf-8a44e7a67f



########################### upload multiple datsets

make_upload_task <- function(wf_id, dataset_name, image, ..., 
                             memory_request = "500Mi",
                             task_name = "bq-upload",
                             data_access = get_workflow(wf_id)[["data_access"]]){
  if( !test_string(dataset_name, pattern = "^p[0-9xy]+_") ){
    message("* using BQ dataset: ... ", appendLF = FALSE)
    dataset_name <- .make_bq_dataset_name(dataset_name, data_access)
    message(dataset_name)
    
  }
  assert_string(dataset_name, pattern = "^p[0-9xy]+_")
  assert_string(wf_id)
  
  make_task(
    sprintf("python /app/cytobigquery/exec_service.py single_cell %s %s --verbose -sc -ti=%s -tm=124Mi",
            # NOTE: the uploader service is working on the old assumption of 2 workflow IDs
            wf_id,
            dataset_name,
            image),
    # CC config
    inputs = list(),
    outdir = "output/",
    image = image,
    memory_request = memory_request,
    task_name = task_name,
    ...
  )
  
}
# the first is UC smilie
combine_ds <- sc_api_save_analysis(c('wf-7a96dd7824','wf-7445c55bf9','wf-3cb167f2a0'), dataset_name = "p00_ulcerative_colitis")






