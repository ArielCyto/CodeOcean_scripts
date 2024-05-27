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

## old object for CDi Krzak
#splitting the object
tmp <- get_task_outputs("wf-66cad48543","0")
data <- readRDS(tmp)
View(data.frame(unique(data$biologist_annotation)))
# BQ upload: wf-66ecf2dcef

## old object for CDi kong
#splitting the object
tmp <- get_task_outputs("wf-f720070d5c","0")
data <- readRDS(tmp)
View(data.frame(unique(data$biologist_annotation)))
# BQ upload: wf-12bcd0761b



# Yuval wf
tmp <- get_task_outputs("wf-f720070d5c","0")
data <- readRDS(tmp)
View(data.frame(unique(data$biologist_annotation)))

# find Yuval's cttest
library(cytoreason.ccm.pipeline)
kong_cdi_config <- get_task_inputs(read_asset('wf-58591b8912'))

