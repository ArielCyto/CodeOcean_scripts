library(cytoreason.ccm.pipeline)
yuvals_config_CDi <- read_asset(get_task_inputs('wf-58591b8912',0)['config'])
# change the internal wf (yuvals_config_CDi$asset_id) so it will have the fixed names
yuvals_config_CDi$asset_id <- "wf-12bcd0761b"

# check the configuration locally
sc_validate_dataset_configuration(yuvals_config_CDi)
# check that a valid dataset can be loaded from the config
pb_object <- prepare_single_cell_dataset(yuvals_config_CDi)

# check a particular analysis model (if failing)
obj1 <- prepare_single_cell_dataset(yuvals_config_CDi, validate = FALSE)
str(build_analysis_model_metadata(yuvals_config_CDi$analysis_model$CDi, obj1))

sc_api_run_analysis("CD-IL", yuvals_config_CDi, tags = list(message = "CD ileum_kong_fixed_070524"))
# [Sun May  5 13:22:56 2024] run_analysis - wf-58b89c6bd1
# [Tue May  7 10:51:49 2024] run_analysis - wf-0e1b7bab54

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


wf_id='wf-0e1b7bab54'
dataset_name = 'p00_cd_ileum' # 



command <- paste("python /app/cytobigquery/exec_service.py single_cell", wf_id, dataset_name, "--verbose -sc -ti=eu.gcr.io/cytoreason/cd-py-bigquery:develop_latest -tm=124Mi")
res <-
  run_command_dist(
    command,
    image = "eu.gcr.io/cytoreason/cd-py-bigquery:SUP-4069_latest",#"eu.gcr.io/cytoreason/cd-py-bigquery:develop_latest",
    memory_request = "25Gi",
    force_execution = FALSE,
    replace_image_tags = TRUE
  )
# Cyto-CC workflow: wf-6403792b0a should be version 52
# Cyto-CC workflow: wf-d9d4202ad9 should be version 54