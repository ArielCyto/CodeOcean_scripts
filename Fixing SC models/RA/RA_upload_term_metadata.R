task_image = 'eu.gcr.io/cytoreason/cd-py-bigquery:SUP-4069_latest'
#####################

memory_request = '124Mi'
# ARGS (positional) #
#####################
service = 'term_metadata'
dataset = "p01_ra_synovium"
text_wf_id ="wf-700eb1cff3"  #"wf-e74f542b22" #change to the wf you got when running update_sphere_disease_model_entry (e.g. 'wf-dd8723885e')

###################
# ARGS (optional) #
###################
verbose = '--verbose'
ti = sprintf('-ti=%s', task_image)
tm = sprintf('-tm=%s', memory_request)###########
# command #
###########
command <- sprintf('python /app/cytobigquery/exec_service.py %s %s %s %s %s %s',
                   service, text_wf_id, dataset, verbose, ti, tm)

############
# workflow #
############
tags=list(list(name="de-process", value='bigquery-upload'),
          list(name="service", value=service),
          list(name="export_workflow", value=text_wf_id),
          list(name="target_dataset", value=dataset)
)
task_env_vars = list(list("name"="DE_PROCESS","value"="BigQuery_upload")
)
wf <- run_command_dist(command,
                       outdir = "./output/",
                       image = task_image ,
                       task_env_vars = task_env_vars,
                       tags=tags,
                       memory_request = memory_request)

# Cyto-CC workflow: wf-eb7a4a7915 old
# Cyto-CC workflow: wf-b677c574d1 new