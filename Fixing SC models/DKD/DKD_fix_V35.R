### env ###
install.packages(c("cytoreason.deconvolution", "cytoreason.datasource", "cytoreason.dbi", "cytoreason.integration",'curl','rlang'))
install.packages('unixtools', repos = 'http://www.rforge.net/')
library(unixtools)
unixtools::set.tempdir("/scratch")
library(curl)

install.packages("googledrive", repos = "https://cloud.r-project.org")# to update googledrive to >= 2.0.0
# restart #
library(rlang)
library(googledrive)

setwd("~/capsule/code/analysis-p01-public-sphere/notes/sphere-migration")
pkgmaker::irequire('future')

install.packages("bigrquery", repos = "https://cloud.r-project.org")
library("bigrquery")
library(cytoreason.cc.client)
devtools::load_all()
MODEL_METADATA_FILE("sphere-model-metadata.csv") #press twice!
devtools::load_all("~/capsule/code/analysis-p01-public-sphere/") 
library(stringr)

### generate term list ### 
model_metadata_dkd = list("DKD",
                          "DKD_Advanced",
                          DKD_Advanced = list("term"="advanced_vs_HC", "contrast_type"="disease_vs_control:advanced"),
                          #  "DKD_ANCA_Associated_Vasculitis",
                          "DKD_Chronic_Kidney_Disease",
                          # "DKD_Cortex",
                          # DKD_Cortex = list("term"="DZ_cortex_vs_HC_cortex", "contrast_type"="disease_vs_control:cortex"), 
                          "DKD_diabetes_mellitus",
                          "DKD_DZ_tissue",
                          "DKD_Early",
                          DKD_Early = list("term"="early_vs_HC", "contrast_type"="disease_vs_control:early"),
                          "DKD_Focal_Seg_Glome",
                          "DKD_Glomerulus",
                          DKD_Glomerulus = list("term"="DZ_glome_vs_HC_glome", "contrast_type"="disease_vs_control:glomerulus"), 
                          "DKD_HC_tissue",
                          "DKD_IgA_Glome",
                          "DKD_Kidney",
                          DKD_Kidney = list("term"="DZ_kidney_vs_HC_kidney", "contrast_type"="disease_vs_control:kidney"),
                          "DKD_Kidney_Cancer",
                          "DKD_Lipoid_Nephrosis",
                          "DKD_Membranous_Glome",
                          "DKD_Rapidly_Prog_Glome",
                          "DKD_Renal_Hypertension",
                          "DKD_stage",
                          "DKD_Systemic_Lupus_Erythematosus",
                          "DKD_Thin_Membrane_Disease",
                          "DKD_Tubule",
                          DKD_Tubule = list("term"="DZ_tubule_vs_HC_tubule", "contrast_type"="disease_vs_control:tubule"))


### generate term_md table ### 
devtools::load_all("~/capsule/code/cytoreason.ccm.pipeline/")
# term_md_3 <- build_designModelTermMetadata("dm://p00_diabetic_kidney:ccm:3",
#                                            drop = NA, 
#                                            fill_default_variables = TRUE,
#                                            model_metadata = model_metadata_dkd, 
#                                            group_level_variables = "dkd_160524 (1).csv" # Use the config link from google instead of a local config file if possible.
# )

term_md_5 <- build_designModelTermMetadata("dm://p00_diabetic_kidney:ccm:5",
                                           drop = NA, 
                                           fill_default_variables = TRUE,
                                           model_metadata = model_metadata_dkd, 
                                           group_level_variables = "dkd_160524 (1) (1).csv" # Use the config link from google instead of a local config file if possible.
)

### check the term_md and keep the csvs ### 
assert_true(all(is.na(term_md_5[, "term_error"])))
write.csv(term_md_5, paste0(getwd(),"/","DKD_06062024_unfiltered.csv"))
term_md_BC_filtered <- filter_small_meta_term_id(term_md_5, threshold = 2L)
validate_term_metadata(term_md_BC_filtered)
write.csv(term_md_BC_filtered, paste0(getwd(),"/","DKD_06062024_filtered.csv"))

### upload new term md version to BQ
devtools::load_all() # branch: issue/SUP-3758-export-TMEs
update_sphere_disease_model_entry(term_metadata = term_md_BC_filtered,
                                  image = "feature_SUP_3758_augment_contrast_term_metadata@0.63.0",
                                  disease_model = "DKD-KD",
                                  tags = list(message = "export DKD test 5"),
                                  group_level_variables = "dkd_160524 (1) (1).csv") # the csv of V5

# https://cyto-cc.cytoreason.com/workflow/wf-9669fd59e8 - V3
# https://cyto-cc.cytoreason.com/workflow/wf-99aac7266f - V5
