config <- sc_api_configuration(
  # asset_id [mandatory]: workflow ID of an exported single-cell dataset as 
  # output by `sc_api_save_data()`
  asset_id = "wf-1b14933a5b",
  cell_annotation = "final_annotation", # or: "final_2023"
  # experiment_id [optional]: experiment identifier.
  # If not provided, the default is to use the asset_id.
  # It does not need to be a GSE ID. 
  # Sample annotations are looked up in the Currator using this key, but the 
  # pipeline will not fail if it does not find any
  experiment_id = "GSE183276",
  
  # platform_id [optional]: platform ID, if not provided it will be "sc-rnaseq" 
  # Note however that the platform ID in the text files used by the platform API will  
  # always be set to "sc-rnaseq".
  platform_id = "GPL24676",
  
  # biosample_id [optional]: a string that indicates the name of the cell annotation
  # variable (in the cell metadata), which holds the bio-sample identifier, 
  # i.e. the key that is used to aggregate the single-cell data into a 
  # pseudobulk data
  biosample_id = "sample_id",
  # biosample_id = "donor_id",
  
  # sample_annotations [optional]:
  sample_annotations = list(
    # condition: values must be from the DOID ontology used in the Curator.
    # Here we can remap invalid values in the original data into valid ontology terms.
    # Values that are not remapped remain untouched.
    # Format: list(<original_variable_name> = c(new_value1 = old_value1, new_value2 = old_value2))
    condition = list(
      condition = c(control = "control",
                    "diabetic kidney disease" = "diabetic kidney disease")
    ),
    # tissue: values must be from the BTO ontology used in the Curator.
    # Here we can remap invalid values in the original data into valid ontology terms 
    # See example on 'condition' for format description
    tissue = list(
      tissue = c(
        "kidney" = "kidney"
      )
    ),
    # 
    # sample_classification = list(
    #   sample_classification = c("Lesion" = "Lesion",
    #   Normal = "Normal")
    # ),
    # directly use column `donor_id` as `subject_id`
    # subject_id = "orig.ident" # use this if the last horse fell ***
    subject_id = "subject_id"
    # # a custom variable as a cross-product of 2 other variables
    # custom_variable = list(
    #   type = c("NonI", "Infl", "Heal"),
    #   disease = c(healthy = "normal",
    #               "CD" = "Crohn disease")
    # )
  ),
  term_metadata = list(
    DKD = list(
      DZ_vs_HC = list(
        contrast_type = "disease_vs_control"
      )
    )
  ),
  analysis_model = list(
    DKD = list(
      group = list(
        # this will fit the same term as use-case 2, but label it as "DZ_vs_HC"
        condition = c(HC = "control", DZ = "diabetic kidney disease")
      )#,
      # pairing = NULL,
      #covariates = c("assay","layer")
      #covariates = "layer"
    )
  ),
  parameters = list(
    model=list(
      services = c("ct_test","gx_diff"),
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
str(build_analysis_model_metadata(config$analysis_model$DKD, obj))

# launch piepline
sc_api_run_analysis("DKD-KD", config, tags = list(message = "Ariel names fixes 240424"))
# [Wed Apr 24 09:24:08 2024] run_analysis - wf-880653748b

model_key = "DKD-KD"

command <- paste("python /app/cytobigquery/exec_service.py single_cell", wf_id = "wf-880653748b", dataset_name = "p00_diabetic_kidney", "--verbose -sc -ti=eu.gcr.io/cytoreason/cd-py-bigquery:develop_latest -tm=124Mi")
res <-
  run_command_dist(
    command,
    image = "eu.gcr.io/cytoreason/cd-py-bigquery:develop_latest", #"eu.gcr.io/cytoreason/cd-py-bigquery:DI-640_latest",
    memory_request = "25Gi",
    force_execution = FALSE,
    replace_image_tags = TRUE
  )
res
# Cyto-CC workflow: wf-42f1f56850
