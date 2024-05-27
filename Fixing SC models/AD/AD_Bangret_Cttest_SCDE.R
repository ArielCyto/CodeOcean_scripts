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
  asset_id = "wf-41a26790c9", #"wf-6211027953",
  cell_annotation = "final_annotation",
  experiment_id = "Bangret_AD",
  # platform_id [optional]: platform ID, if not provided it will be "sc-rnaseq" 
  # Note however that the platform ID in the text files used by the platform API will  
  # always be set to "sc-rnaseq".
  platform_id = "10x 3' v2", 
  # biosample_id [optional]: a string that indicates the name of the cell annotation
  # variable (in the cell metadata), which holds the bio-sample identifier, 
  # i.e. the key that is used to aggregate the single-cell data into a 
  # pseudobulk data
  biosample_id = 'sample_id',
  sample_annotations = read_data("~/capsule/code/service-single-cell/Bangret_metadata_080524.csv"),
  
  term_metadata = list(
    AD = list(
      DZ_vs_HC = list(
        contrast_type = "disease_vs_control"
      )
    ),
    AD1y_vs_AD16w = list(
      # will only export this contrast term within this comparison
      AD1y_vs_AD16w = list(contrast_type = "post_vs_post")
    )
    ,
    AD16w_vs_DZ = list(
      # will only export this contrast term within this comparison
      AD16w_vs_DZ = list(contrast_type = "post_vs_pre")
    ),
    AD1y_vs_DZ= list(
      # will only export this contrast term within this comparison
      AD1y_vs_DZ = list(contrast_type = "post_vs_pre")
    )
  ),
  analysis_model = list(
    AD = list(
      group = list(
        type = c(HC = "HC", DZ = "AD_UT")
      )
    ),
    AD1y_vs_AD16w = list(
      group = list(
        type = c(AD16w = "AD_16w",AD1y = "AD_1y")
      )
    ) ,
    AD16w_vs_DZ = list(
      group = list(
        type = c(DZ = "AD_UT",AD16w = "AD_16w")
      )
    ),
    AD1y_vs_DZ = list(
      group = list(
        type = c(DZ = "AD_UT", AD1y = "AD_1y")
      )
    )
  ),
  parameters = list(
    model=list(
      services = c("gx_diff", "ct_test"),
      gx_diff=list(),
      ct_test=list()
    )
  )
)


# check the configuration locally
sc_validate_dataset_configuration(config)
# check that a valid dataset can be loaded from the config
pb_object <- prepare_single_cell_dataset(config)

# check a particular analysis model (if failing)
obj1 <- prepare_single_cell_dataset(config, validate = FALSE)
pData(obj1)


sc_api_run_analysis("AD-SK", config, tags = list(message = "AD Bangret final with new metadata table", test_run = 3),
                    image = "eu.gcr.io/cytoreason/ci-cytoreason.single.cell-package:develop_latest")
# [Tue Mar 12 17:26:27 2024] run_analysis - wf-66f77c345c
# [Wed Mar 20 12:18:56 2024] run_analysis - wf-e5578fe1a4
# [Wed Mar 27 15:49:44 2024] run_analysis - wf-47cd96ba5c - new cell names
# [Wed May  8 15:16:35 2024] run_analysis - wf-6fb4cb5525 - add post_vs_pre, drug
# [Thu May  9 06:27:30 2024] run_analysis - wf-c08f099c86 - add post_vs_pre, drug, dosage, time
# [Thu May  9 08:06:42 2024] run_analysis - wf-fca22830bc - change directions
# [Thu May  9 14:06:02 2024] run_analysis - wf-f460ec08aa - rechekc


# upload to BQ
# wf_id='wf-66f77c345c'
# dataset_name = 'p00_atopic' 
# IMAGE <- 'eu.gcr.io/cytoreason/cd-py-bigquery:develop_latest'
# 
# command <- paste("python /app/cytobigquery/exec_service.py single_cell", wf_id, dataset_name, "--verbose -sc -ti=eu.gcr.io/cytoreason/cd-py-bigquery:develop_latest -tm=124Mi")
# res <-
#   run_command_dist(
#     command,
#     image = IMAGE,
#     memory_request = "25Gi",
#     force_execution = FALSE,
#     replace_image_tags = TRUE
#   )
# # Cyto-CC workflow: wf-08bdd92ff7
# 

# take the 0 mission wf to p01 capsule
combine_ds <- sc_api_save_analysis(c('wf-84e233def4','wf-f460ec08aa'), dataset_name = "p01_atopic_dermatitis")
# Cyto-CC workflow: wf-2f7ec94d78
# Cyto-CC workflow: wf-b4c2bb4d58
# Cyto-CC workflow: wf-85cd4bb5eb
# Cyto-CC workflow: wf-d75de590fc