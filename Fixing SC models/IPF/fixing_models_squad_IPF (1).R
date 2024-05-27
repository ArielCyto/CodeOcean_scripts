# fix names of cells

#env 
library(cytoreason)
load_all("capsule/code/service-single-cell-preprocessing/")

library(devtools)
library(Seurat)
library(ggplot2)
library(future)
library(RColorBrewer)
# library(SingleR)
library(SeuratDisk)
library(cytoreason.cc.client)
library(tidyr)
library(plyr)
library(AnnotationDbi)
library(ggplot2)

library("SeuratDisk")
if (!require("annotables")){
  devtools::install_github("stephenturner/annotables", upgrade = FALSE)
}
library("annotables")
if (!require("unixtools")){
  remotes::install_github("s-u/unixtools", upgrade = FALSE)
}
unixtools::set.tempdir("/scratch")

#defile output directory
out_path <- "/scratch/sc_pipeline/"

#the main image used in the pipeline
IMAGE <-"eu.gcr.io/cytoreason/ci-cytoreason.single.cell.preprocessing-package:master_latest"











save_to_cyto_cc <- function(obj, filename = "seurat_object.rds", data_access = cytoreason.cc.client:::default_data_access(), tag = "none"){
  # saving the output file to cyto-cc
  run_command_dist(
    command = sprintf('cp /cyto_cc/inputs/obj.rds "output/%s"',
                      filename),
    image = "alpine",
    inputs = list(
      input_local_obj(obj, 'obj.rds')
    ),
    memory_request = '2Gi',
    force_execution = FALSE,
    replace_image_tags = FALSE,
    data_access = data_access,
    tags = list(list(name = "save", value = "dataset"),
                list(name = "name", value = tag))
  )
  
}

##preparing IPF's dataset for DA analysis

#splitting the object
tmp <- get_task_outputs("wf-9db04f14a4","0")
data <- readRDS(tmp)

# see all optional cell names
data.frame(unique(data$new_annotations))

data$new_annotations[which(data$new_annotations %in% "CD1a-positive Langerhans cell")]<- "CD1a+ Langerhans cell"
data$new_annotations[which(data$new_annotations %in% "CD4-positive, alpha-beta memory T cell")]<- "CD4+ memory T cell"
data$new_annotations[which(data$new_annotations %in% "Inflammatory macrophage")]<- "inflammatory macrophage"
data$new_annotations[which(data$new_annotations %in% "Th17")]<- "T-helper 17 cell"
data$new_annotations[which(data$new_annotations %in% "mature conventional dendritic cell")]<- "migratory conventional dendritic cell"
data$new_annotations[which(data$new_annotations %in% "elicited macrophage")]<- "monocyte derived macrophage"

data.frame(unique(data$new_annotations))

save_to_cyto_cc(obj = data, tag = "IPF_with_fixed_names_150424")
  # Cyto-CC workflow: wf-52c69af287

sc_fit <- as_sc_pipeline_run(
  data = "wf-52c69af287"
)

wf_export <- sc_api_save_data(
  data = sc_fit,dataset_id = 'GSE136831', image= IMAGE
)
# [Wed Apr 10 10:20:45 2024] save_data_sc_pipeline - wf-6eb117c0c8
# [Tue Apr 16 08:31:36 2024] save_data_sc_pipeline - wf-9722097131


# Now take this wf and use it for cttest
