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
# out_path <- "~/capsule/scratch/sc_pipeline/ct_test_csde/"
# IMAGE <- "eu.gcr.io/cytoreason/ci-cytoreason.single.cell-package:MF_32_SC_gene_diff_wrapper_for_CCM_output_latest"

# term metadata
krzak = read.csv("~/capsule/code/service-single-cell/Krzak study Annotations - SC for CDi - updated.csv")
krzak$Sample.Classification[krzak$Condition == "healthy"] = "Normal"
krzak$Sample.Classification_untreated = krzak$Sample.Classification
krzak$Sample.Classification_untreated[krzak$Drug != ""] = ""
krzak$sample_id = krzak$Subject.ID

krzak$condition_untreated = krzak$Condition
krzak$condition_untreated[krzak$Drug != ""] = ""

write.csv(krzak, "~/capsule/code/service-single-cell/Krzak study Annotations - SC for CDi - updated.csv")


## run ct-test and csde -----
config <- sc_api_configuration(
  asset_id = "wf-66ecf2dcef",
  cell_annotation = "biologist_annotation",
  experiment_id = "CDi_krzak",
  platform_id = "sc-rnaseq",
  biosample_id = "sanger_sample_id",
  sample_annotations = read.csv("~/capsule/code/service-single-cell/Krzak study Annotations - SC for CDi - updated.csv", na.strings=''),
  term_metadata = list(CDi = list(DZ_vs_HC = list(contrast_type = "disease_vs_control")),
                       CDi_inf = list(inflamed_vs_HC = list(contrast_type = "inflamed_vs_healthy"))
                       
  ),
  
  analysis_model = list(CDi = list(group = list(condition_untreated = c(HC = "healthy", DZ = "Crohn's disease"))),
                        CDi_inf = list(group = list(sample_classification_untreated = c(HC = "Normal", inflamed = "Inflamed")))
  ),
  parameters = list(model=list(services = c("ct_test", "gx_diff"), ct_test=list(), gx_diff=list(min_sample_size=15))))

#### check the configuration locally ----
sc_validate_dataset_configuration(config)

#### check that a valid dataset can be loaded from the config ----
pb_object <- prepare_single_cell_dataset(config)

# check a particular analysis model (if failing)
obj <- prepare_single_cell_dataset(config, validate = FALSE)
str(build_analysis_model_metadata(config$analysis_model$CDi, obj))

sc_api_run_analysis("CD-IL", config, tags = list(message = "CD ileum_krzak_fixed_050524")) # ,image = IMAGE
# [Tue Mar 26 21:41:39 2024] run_analysis - wf-623d24045b failed
# [Wed Mar 27 07:33:27 2024] run_analysis - wf-bac5686d49 V
# [Wed Mar 27 11:29:03 2024] run_analysis - wf-7c23710c30
# [Sun May  5 13:22:43 2024] run_analysis - wf-7645848d04 - Ariel fixes









command <- paste("python /app/cytobigquery/exec_service.py single_cell", wf_id = "wf-7c23710c30",
                 dataset_name ="p00_cd_ileum", "--verbose -sc -ti=eu.gcr.io/cytoreason/cd-py-bigquery:develop_latest -tm=124Mi")
res <-
  run_command_dist(
    command,
    image = "eu.gcr.io/cytoreason/cd-py-bigquery:develop_latest",
    memory_request = "25Gi",
    force_execution = FALSE,
    replace_image_tags = TRUE
  )

# Cyto-CC workflow: wf-c5fc53e068 -> succeed but not the right dataset
# Cyto-CC workflow: wf-f1db42e82f

# wf-7c23710c30



# dataset_id = "cdi_krzak",

# multiple datasets
# the first is CDI kong, second cdi krzak
combine_ds <- sc_api_save_analysis(c('wf-0e1b7bab54','wf-7645848d04'), dataset_name = "p00_cd_ileum")
# Cyto-CC workflow: wf-15c68a992f 06.05.24, will be version 53
# Cyto-CC workflow: wf-dffd301159 07.05.24, will be version 55
