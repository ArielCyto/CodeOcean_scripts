install.packages('cytoreason.cc.client', repos = c(CRRAN = "https://crran.cytoreason.com/R/latest"))
library(cytoreason.cc.client)
library(plyr)
library(cytoreason)
library(cytoreason.io)
unixtools::set.tempdir("/scratch")

# usethis::edit_r_environ()
load_all("/root/capsule/code/service-single-cell/")



# check the columns in the original file
final_sc_object = readRDS(get_task_outputs('wf-52c69af287',"0")['seurat_object.rds']) # save_data_sc_pipeline - wf-9722097131
unique(final_sc_object$biosample_id)
final_sc_object$sample_id <- final_sc_object$biosample_id
unique(final_sc_object$sample_id)


config <- sc_api_configuration(
  # asset_id [mandatory]: workflow ID of an exported single-cell dataset as
  # output by `sc_api_save_data()`
  asset_id = "wf-9722097131", # "wf-6eb117c0c8", ########## "wf-ac5c91fe6a", # wf-ac5c91fe6a, wf-5861c60f41 -> add curator subject column
  cell_annotation = "new_annotations", # or: "final_2023"
  experiment_id = "GSE136831-", #"kong_ileum",
  # pseudobulk data
  biosample_id = "library_identity", # , Sample ID
  
  sample_annotations = read_data("~/capsule/code/service-single-cell/GSE136831_metadata_curator.csv"),
  
  # sample_annotations [optional]:
  # sample_annotations = list(
  #   condition = list(
  #     disease_identity = c(control = "Control",   #names in the config: control, healthy, idiopathic pulmonary fibrosis
  #                          "idiopathic pulmonary fibrosis" = "IPF")
  #   ),
  #   # tissue: values must be from the BTO ontology used in the Curator.
  #   tissue = list(
  #     tissue = c(
  #       "lung" = "lung"
  #     )
  #   ),
  #   # directly use column `donor_id` as `subject_id`
  #   subject_id = "subject_identity"
  # ),
  term_metadata = list(
    IPF = list(
      DZ_vs_HC = list(
        contrast_type = "disease_vs_control"
      )
    )
  ),
  analysis_model = list(
    IPF = list(
      group = list(
        # this will fit the same term as use-case 2, but label it as "DZ_vs_HC"
        condition = c(HC = "control", DZ = "idiopathic pulmonary fibrosis")
      )
    )
  ),
  parameters = list(
    model=list(
      services = c("gx_diff", "ct_test"),
      gx_diff=list(min_sample_size=15),
      ct_test=list()
    )
  )
)


# check the configuration locally
sc_validate_dataset_configuration(config)
# check that a valid dataset can be loaded from the config
pb_object <- prepare_single_cell_dataset(config)
# check a particular analysis model (if failing)

obj <- prepare_single_cell_dataset(config, validate = FALSE)
pData(obj)

str(build_analysis_model_metadata(config$analysis_model$IPF, obj))


# launch piepline
sc_api_run_analysis("IPF-LU", config, tags = list(message = "ipf lung with metadata, Ariel fixes change ds name"))
# [Sun Dec 10 08:49:08 2023] run_analysis - wf-187177d310
# [Sun Jan 21 17:42:18 2024] run_analysis - wf-042614961e -> probably updated code (0.01 bug fix)

# [Sun Feb  4 11:12:07 2024] run_analysis - wf-797abbe115, add curator metadata based on gsm biosample column
# [Tue Feb  6 12:37:44 2024] run_analysis - wf-6d3684065e, add manual curator metadata, based on sample_id (library_id)
# [Mon Feb 12 08:57:49 2024] run_analysis - wf-5bc2dca1d2, add manual curator metadata, change dataset to GSE_try while waiting for curator vs manual fix
# [Wed Apr 10 12:19:53 2024] run_analysis - wf-7edac959bd - Ariel fixes
# [Tue Apr 16 10:26:31 2024] run_analysis - wf-f8a1763546 - Ariel fixes ver2
# [Thu Apr 18 12:29:00 2024] run_analysis - wf-c0b3d4b370 - Ariel fixes ver3 - the "-" name

model_key = "IPF-LU"

command <- paste("python /app/cytobigquery/exec_service.py single_cell", wf_id = "wf-c0b3d4b370", dataset_name = "p00_idiopathic_pulmonary_fibrosis", "--verbose -sc -ti=eu.gcr.io/cytoreason/cd-py-bigquery:MF-67_latest -tm=124Mi")
res <-
  run_command_dist(
    command,
    image = "eu.gcr.io/cytoreason/cd-py-bigquery:MF-67_latest", #"eu.gcr.io/cytoreason/cd-py-bigquery:DI-640_latest",
    memory_request = "25Gi",
    force_execution = FALSE,
    replace_image_tags = TRUE
  )
res
# wf-01f7b34e0b -> input: wf-187177d310
# wf-f7faff76c7 -> input: wf-042614961e
# Cyto-CC workflow: wf-9c1eb08635
# Cyto-CC workflow: wf-238e222937, wf-5bc2dca1d2 change dataset to GSE_try while waiting for curator vs manual fix
# Cyto-CC workflow: wf-65cec382ca - Ariel fixes
# Cyto-CC workflow: wf-53758f3ba9 - Ariel fixes ver2
# Cyto-CC workflow: wf-7db7dc0851 - Ariel fixes ver3 - the "-" name