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

## old object for Alopecia
#splitting the object
tmp <- get_task_outputs("wf-6eee71afda","0")
data <- readRDS(tmp)

View(data.frame(unique(data$final_annotation)))

# see all optional cell names
data$final_annotation[which(data$final_annotation %in% "B cell, CD19-positive")]<- "CD19+ B cell"
data$final_annotation[which(data$final_annotation %in% "central memory CD8-positive, alpha-beta T cell")]<- "central memory CD8+ T cell"
data$final_annotation[which(data$final_annotation %in% "central memory CD4-positive, alpha-beta T cell")]<- "central memory CD4+ T cell"

data$final_annotation[which(data$final_annotation %in% "endothelial cell of lymphatic vessel")]<- "lymphatic endothelial cell"
data$final_annotation[which(data$final_annotation %in% "hair follicle dermal papilla cell")]<- "hair follicle dermal papilla cell"
data$final_annotation[which(data$final_annotation %in% "mature NK T cell")]<- "natural killer T cell"

data$final_annotation[which(data$final_annotation %in% "mature conventional dendritic cell")]<- "migratory conventional dendritic cell"
data$final_annotation[which(data$final_annotation %in% "resident memory CD4-positive, alpha-beta T cell")]<- "resident memory CD4+ T cell"
data$final_annotation[which(data$final_annotation %in% "resident memory CD8-positive, alpha-beta T cell")]<- "resident memory CD8+ T cell"

View(data.frame(unique(data$final_annotation)))


View(data.frame(unique(data$final_annotation))) #data is already fixed :)

save_to_cyto_cc(obj = data, tag = "ALO_with_fixed_names_240424")
# Cyto-CC workflow: wf-9be824fbb3

sc_fit <- as_sc_pipeline_run(
  data = "wf-9be824fbb3"
)

wf_export <- sc_api_save_data(
  data = sc_fit,dataset_id = 'GSE212447', image= IMAGE
)
# [Wed Apr 24 10:31:43 2024] save_data_sc_pipeline - wf-5e233d5aa1

# Now take this wf and use it for cttest
