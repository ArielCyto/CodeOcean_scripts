
try_add_metadata_curator = function(){
  add_biosampleID = function(GSE136831_curator){
    library(Seurat)
    library(cytoreason.cc.client)
    library(qs)
    
    final_sc_object = readRDS(get_task_outputs('wf-b5c80ba585',"0")['seurat_object.rds'])
    
    # map columns Sample ID from curator to Subject_Identity (=subject id in curator)
    rownames(GSE136831_curator) = GSE136831_curator$subject_id
    final_sc_object$biosample_id = GSE136831_curator[unname(final_sc_object$Subject_Identity), 'samples_id']
    return(final_sc_object)
  }
  
  write.csv(final_sc_object@meta.data, 'metadata/ipf_metadata.csv')
  
  # how many unique in library_identity? why avital chose it for ct test?
  metadata = read.csv('metadata/ipf_metadata.csv')
  
  # filter for manual input for ct test
  small_metadata = metadata[, c('Subject_Identity', 'Disease_Identity', 'Library_Identity', 'LibraryIDs', "Sex", "Age", "Race", "Smoked", "Condition", "Tissue" )]
  # remove duplicate rows
  small_metadata <- distinct(small_metadata) 
  write.csv(small_metadata, 'metadata/ipf_meta_small.csv', row.names = F)
  
  length(unique(metadata$Subject_Identity)) # 60
  length(unique(metadata$Library_Identity)) # 83
  length(unique(GSE136831_curator$samples_id)) # 78
  length(unique(GSE136831_curator$subject_id)) # 78
  length(unique(metadata$biosample_id)) # 60 = 78- 18 COPD
  
  
  tags <- list(list(name="owner", value="Chaya Barbolin"),
               list(name="model", value="IPF sc p00"),
               list(name="data", value="add biosample_id from curator"),
               list(name="data_access", value="p00"))
  wf <- run_function_dist(add_biosampleID, 
                          GSE136831_curator=GSE136831_curator, 
                          tags = tags, 
                          image ="eu.gcr.io/cytoreason/ci-cytoreason.single.cell-package:sc_pipeline_e2e_0.0.7.3", memory_request = "80Gi")
  wf
  # Cyto-CC workflow: wf-9db04f14a4
  
  sc_fit <- as_sc_pipeline_run(data = "wf-9db04f14a4")
  wf_export <- sc_api_save_data(dataset_name = "p00_sc_idiopathic_pulmonary_fibrosis",
                                data = sc_fit)
  # [Sat Feb  3 00:03:34 2024] save_data_sc_pipeline - wf-5861c60f41
  
}





GSE136831_curator = read.csv('~/capsule/code/service-single-cell/GSE136831_metadata_curator.csv')
length(unique(GSE136831_curator$Subject.ID)) # 60
length(unique(GSE136831_curator$Sample.ID)) # 83

config <- sc_api_configuration(
  # asset_id [mandatory]: workflow ID of an exported single-cell dataset as
  # output by `sc_api_save_data()`
  asset_id = "wf-9722097131", # "wf-6eb117c0c8", ########## "wf-ac5c91fe6a", # wf-ac5c91fe6a, wf-5861c60f41 -> add curator subject column
  cell_annotation = "new_annotations", # or: "final_2023"
  experiment_id = "GSE136831_try", #"kong_ileum",
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
str(build_analysis_model_metadata(config$analysis_model$IPF, obj))


# launch piepline
sc_api_run_analysis("IPF-LU", config, tags = list(message = "ipf lung with metadata, Ariel fixes"))
# [Sun Dec 10 08:49:08 2023] run_analysis - wf-187177d310
# [Sun Jan 21 17:42:18 2024] run_analysis - wf-042614961e -> probably updated code (0.01 bug fix)

# [Sun Feb  4 11:12:07 2024] run_analysis - wf-797abbe115, add curator metadata based on gsm biosample column
# [Tue Feb  6 12:37:44 2024] run_analysis - wf-6d3684065e, add manual curator metadata, based on sample_id (library_id)
# [Mon Feb 12 08:57:49 2024] run_analysis - wf-5bc2dca1d2, add manual curator metadata, change dataset to GSE_try while waiting for curator vs manual fix
# [Wed Apr 10 12:19:53 2024] run_analysis - wf-7edac959bd - Ariel fixes
# [Tue Apr 16 10:26:31 2024] run_analysis - wf-f8a1763546 - Ariel fixes ver2

model_key = "IPF-LU"

command <- paste("python /app/cytobigquery/exec_service.py single_cell", wf_id = "wf-f8a1763546", dataset_name = "p00_idiopathic_pulmonary_fibrosis", "--verbose -sc -ti=eu.gcr.io/cytoreason/cd-py-bigquery:MF-67_latest -tm=124Mi")
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