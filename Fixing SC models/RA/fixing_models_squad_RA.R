
library(plyr)
library(cytoreason)
library(cytoreason.io)

# usethis::edit_r_environ()
load_all("/root/capsule/code/service-single-cell/")

install.packages('cytoreason.cc.client', repos = c(CRRAN = "https://crran.cytoreason.com/R/latest"))
library(cytoreason.cc.client)

#library(deaR)

# verify DOID
DOID=get_ontology_disease_terms()
grep('rheumatoid arthritis', DOID, value = T)
grep('osteoarthritis', DOID, value = T)


config <- sc_api_configuration(
  # asset_id [mandatory]: workflow ID of an exported single-cell dataset as 
  # output by `sc_api_save_data()`
  asset_id = "wf-f20e189f52", #  - final dataset+annotation
  cell_annotation = "biologist_annotation", # 
  # experiment_id [optional]: experiment identifier.
  # If not provided, the default is to use the asset_id.
  # It does not need to be a GSE ID. 
  # Sample annotations are looked up in the Currator using this key, but the 
  # pipeline will not fail if it does not find any
  experiment_id = "P01_AMP2",
  
  # platform_id [optional]: platform ID, if not provided it will be "sc-rnaseq" 
  # Note however that the platform ID in the text files used by the platform API will  
  # always be set to "sc-rnaseq".
  platform_id = "sc-rnaseq", # Chromium Next GEM Chip G
  
  # biosample_id [optional]: a string that indicates the name of the cell annotation
  # variable (in the cell metadata), which holds the bio-sample identifier, 
  # i.e. the key that is used to aggregate the single-cell data into a 
  # pseudobulk data
  biosample_id = "sample",
  
  # sample_annotations [optional]:
  sample_annotations = read_data("~/capsule/code/service-single-cell/RA SC AMP2 study annotations.csv"),
  term_metadata = list(
    
    RA_OA = list(
      RA_vs_OA = list(
        contrast_type = "disease_vs_disease" # disease_vs_disease disease_vs_control
      )
    )
  ),
  analysis_model = list(
    
    RA_OA = list(
      group = list(
        condition = c(OA = "osteoarthritis", RA = "rheumatoid arthritis") # osteoarthritis
      )#,covariates = 'tissue'
    )
  )#,
  # parameters = list(
  #   model=list(
  #     services = c("gx_diff", "ct_test"),
  #     gx_diff=list(min_sample_size=15),
  #     ct_test=list()
  #   )
  # )
)



# check the configuration locally
sc_validate_dataset_configuration(config)
# check that a valid dataset can be loaded from the config
pb_object <- prepare_single_cell_dataset(config)

# check a particular analysis model (if failing)
obj1 <- prepare_single_cell_dataset(config, validate = FALSE)
str(build_analysis_model_metadata(config$analysis_model$RA_OA, obj1))




# launch piepline
sc_api_run_analysis("RA-SYN", config, tags = list(message = "RA cttest and CSDE, fixed names 170424"))

# [Sun Jan 14 10:25:31 2024] run_analysis - wf-55fef32e7f, with final ct annotation
# [Sun Jan 21 08:34:08 2024] run_analysis - wf-fa13aceb1c, final ct, cttest+CSDE, with default image
# [Thu Jan 25 08:12:24 2024] run_analysis - wf-4d7607b773
# [Sun Feb  4 09:37:31 2024] run_analysis - wf-06daaa2564 - Ariel
# [Mon Feb  5 09:30:15 2024] run_analysis - wf-46d6580298 - Ariel change the names of group A and B
# [Mon Feb  5 12:32:08 2024] run_analysis - wf-3ca70026e3
# [Wed Apr 17 09:26:47 2024] run_analysis - wf-f513108a39 - Ariel fixes 1704

library(SingleCellExperiment)
library(cytoreason.io)
library(cytoreason.ccm.pipeline)
library('cytoreason')
library(pkgmaker)
library(cytoreason.datasource)
library(naturalsort)
library(cytoreason.single.cell)

library(cytoreason)
library(devtools)

library(cytoreason.cc.client)


wf_id='wf-f513108a39'
dataset_name = 'p01_ra_synovium' # 



command <- paste("python /app/cytobigquery/exec_service.py single_cell", wf_id, dataset_name, "--verbose -sc -ti=eu.gcr.io/cytoreason/cd-py-bigquery:SUP-4069_latest -tm=124Mi")
res <-
  run_command_dist(
    command,
    image = "eu.gcr.io/cytoreason/cd-py-bigquery:SUP-4069_latest",#"eu.gcr.io/cytoreason/cd-py-bigquery:DI-640_latest",
    memory_request = "25Gi",
    force_execution = FALSE,
    replace_image_tags = TRUE
  )
# Cyto-CC workflow: wf-8fae6aeccb
# Cyto-CC workflow: wf-e1c1479942
# Cyto-CC workflow: wf-3dc3dd8ca9 - Ariel
# Cyto-CC workflow: wf-93bb51c57f # 4069_latest
# Cyto-CC workflow: wf-672fb9d5a9 - Ariel fixes for RA 1704