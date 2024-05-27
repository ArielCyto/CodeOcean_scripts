library(plyr)
library(cytoreason)
library(cytoreason.io)
# usethis::edit_r_environ()
load_all("/root/capsule/code/service-single-cell/")
# install.packages('cytoreason.cc.client', repos = c(CRRAN = "https://crran.cytoreason.com/R/latest"))
library(cytoreason.cc.client)
#remotes::install_github("s-u/unixtools")
unixtools::set.tempdir("/scratch")


config <- sc_api_configuration(
  asset_id = 'wf-e684206586' ,# "wf-a8e3e424af",
  cell_annotation = "cr_annotation",
  experiment_id = "E-MTAB-8142-",
  # platform_id [optional]: platform ID, if not provided it will be "sc-rnaseq" 
  # Note however that the platform ID in the text files used by the platform API will  
  # always be set to "sc-rnaseq".
  platform_id = "10x 3' v2", 
  # biosample_id [optional]: a string that indicates the name of the cell annotation
  # variable (in the cell metadata), which holds the bio-sample identifier, 
  # i.e. the key that is used to aggregate the single-cell data into a 
  # pseudobulk data
  biosample_id = 'sample_pull', #'sample_pull',
  sample_annotations = read_data("~/capsule/code/service-single-cell/E-MTAB-8142 study annotations_new 2603.csv"),
  
  term_metadata = list(
    AD = list(
      DZ_vs_HC = list(
        contrast_type = "disease_vs_control"
      )
    ),
    L_vs_NL = list(
      # will only export this contrast term within this comparison
      L_vs_NL = list(contrast_type = "active_vs_non.active")
    )
    ,
    L_vs_HC = list(
      # will only export this contrast term within this comparison
      L_vs_HC = list(contrast_type = "active_vs_healthy")
    )
  ),
  analysis_model = list(
    AD = list(
      group = list(
        condition = c(HC = "healthy", DZ = "atopic dermatitis")
      ),covariates = 'sex'
    ),
    L_vs_NL = list(
      group = list(
        sample_classification = c(NL = "Non Lesion", L = "Lesion")
      ), pairing='donor_id'
      #, covariates = 'tissue'
    ) ,
    L_vs_HC = list(
      group = list(
        sample_classification = c(HC = "Normal", L = "Lesion")
      ),covariates = 'sex'
    )
  ),
  parameters = list(
    model=list(
      services = c("gx_diff"),
      gx_diff=list() #,
      # ct_test=list()
    )
  )
)


# check the configuration locally
sc_validate_dataset_configuration(config)
# check that a valid dataset can be loaded from the config
pb_object <- prepare_single_cell_dataset(config)

# check a particular analysis model (if failing)
obj1 <- prepare_single_cell_dataset(config, validate = FALSE)
str(build_analysis_model_metadata(config$analysis_model$L_vs_NL, obj1))
str(build_analysis_model_metadata(config$analysis_model$AD, obj1))- 

pData(obj1)
sc_api_run_analysis("AD-SK", config, tags = list(message = "AD Reynolds after cell names fix and sample_pull column", test_run = 1),
                    image = 'eu.gcr.io/cytoreason/ci-cytoreason.single.cell-package@sha256:f219909ddf492eb8707107cbf4d17f7ee1feaa0a922bcfa4f7d7bbb3f06cc878')
# [Thu Feb 22 12:55:34 2024] run_analysis - wf-2b1a6a0373
# [Sun Feb 25 14:35:39 2024] run_analysis - wf-2ffa3bae8b
# [Mon Feb 26 12:46:15 2024] run_analysis - wf-d37adee38a
# [Thu Feb 29 11:40:53 2024] run_analysis - wf-5367833e45 - sample_pull gx_diff + cttest
# [Sun Mar 10 10:41:35 2024] run_analysis - wf-e02ab0ad4b - sample_pull gx_diff only - failed
# [Sun Mar 10 11:03:25 2024] run_analysis - wf-fd303fcf7e - sample_pull gx_diff only with the previous replacemnt image
# [Thu Mar 28 10:59:52 2024] run_analysis - wf-84e233def4

# upload to BQ
wf_id='wf-fd303fcf7e'
dataset_name = 'p00_atopic' 
IMAGE <- 'eu.gcr.io/cytoreason/cd-py-bigquery:develop_latest'

command <- paste("python /app/cytobigquery/exec_service.py single_cell", wf_id, dataset_name, "--verbose -sc -ti=eu.gcr.io/cytoreason/cd-py-bigquery:develop_latest -tm=124Mi")
res <-
  run_command_dist(
    command,
    image = IMAGE,
    memory_request = "25Gi",
    force_execution = FALSE,
    replace_image_tags = TRUE
  )
# Cyto-CC workflow: wf-9b7f6896fa - Version 11 Tableau
# Cyto-CC workflow: wf-5c9a9ae371 - Version 12 Tableau
# Cyto-CC workflow: wf-0890df293e - Version 13 Tableau
# Cyto-CC workflow: wf-a6de9bcc82 - Version 14 Tableau
