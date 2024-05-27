clustering_object <- readRDS(get_task_outputs('wf-a6d9afce3d', '0')['output.rds'])
# generate_new_column
clustering_object$sample_pull <- clustering_object$SampleID
# combine
clustering_object$sample_pull[clustering_object$sample_pull == 'SKN8090561' | clustering_object$sample_pull == 'SKN8090565'] <- 'SKN8090561_SKN8090565'
clustering_object$sample_pull[clustering_object$sample_pull == 'SKN8090560' | clustering_object$sample_pull == 'SKN8090564'] <- 'SKN8090560_SKN8090564'
clustering_object$sample_pull[clustering_object$sample_pull == 'SKN8090563' | clustering_object$sample_pull == 'SKN8090567'] <- 'SKN8090563_SKN8090567'
clustering_object$sample_pull[clustering_object$sample_pull == 'SKN8090562' | clustering_object$sample_pull == 'SKN8090566'] <- 'SKN8090562_SKN8090566'


clustering_object$sample_pull[clustering_object$sample_pull == 'SKN8090530' | clustering_object$sample_pull == 'SKN8090531'] <- 'SKN8090530_SKN8090531'
clustering_object$sample_pull[clustering_object$sample_pull == 'SKN8090528' | clustering_object$sample_pull == 'SKN8090529'] <- 'SKN8090528_SKN8090529'
clustering_object$sample_pull[clustering_object$sample_pull == 'SKN8090526' | clustering_object$sample_pull == 'SKN8090527'] <- 'SKN8090526_SKN8090527'
clustering_object$sample_pull[clustering_object$sample_pull == 'SKN8090524' | clustering_object$sample_pull == 'SKN8090525'] <- 'SKN8090524_SKN8090525'


clustering_object$sample_pull[clustering_object$sample_pull == 'SKN8090540' | clustering_object$sample_pull == 'SKN8090541'] <- 'SKN8090540_SKN8090541'
clustering_object$sample_pull[clustering_object$sample_pull == 'SKN8090542' | clustering_object$sample_pull == 'SKN8090543'] <- 'SKN8090542_SKN8090543'
clustering_object$sample_pull[clustering_object$sample_pull == 'SKN8090536' | clustering_object$sample_pull == 'SKN8090537'] <- 'SKN8090536_SKN8090537'
clustering_object$sample_pull[clustering_object$sample_pull == 'SKN8090538' | clustering_object$sample_pull == 'SKN8090539'] <- 'SKN8090538_SKN8090539'


clustering_object$sample_pull[clustering_object$sample_pull == 'SKN8090552' | clustering_object$sample_pull == 'SKN8090553'] <- 'SKN8090552_SKN8090553'
clustering_object$sample_pull[clustering_object$sample_pull == 'SKN8090554' | clustering_object$sample_pull == 'SKN8090555'] <- 'SKN8090554_SKN8090555'
clustering_object$sample_pull[clustering_object$sample_pull == 'SKN8090548' | clustering_object$sample_pull == 'SKN8090549'] <- 'SKN8090548_SKN8090549'
clustering_object$sample_pull[clustering_object$sample_pull == 'SKN8090550' | clustering_object$sample_pull == 'SKN8090551'] <- 'SKN8090550_SKN8090551'


clustering_object$sample_pull[clustering_object$sample_pull == '4820STDY7388991' | clustering_object$sample_pull == '4820STDY7388992'
                              | clustering_object$sample_pull == '4820STDY7388993' | clustering_object$sample_pull == '4820STDY7388994'] <- '4820STDY7388991_4820STDY7388992_4820STDY7388993_4820STDY7388994'

clustering_object$sample_pull[clustering_object$sample_pull == '4820STDY7388995' | clustering_object$sample_pull == '4820STDY7388996'
                              | clustering_object$sample_pull == '4820STDY7388997' | clustering_object$sample_pull == '4820STDY7388998'] <- '4820STDY7388995_4820STDY7388996_4820STDY7388997_4820STDY7388998'

clustering_object$sample_pull[clustering_object$sample_pull == '4820STDY7388999' | clustering_object$sample_pull == '4820STDY7389000'
                              | clustering_object$sample_pull == '4820STDY7389001' | clustering_object$sample_pull == '4820STDY7389002'] <- '4820STDY7388999_4820STDY7389000_4820STDY7389001_4820STDY7389002'

clustering_object$sample_pull[clustering_object$sample_pull == '4820STDY7389003' | clustering_object$sample_pull == '4820STDY7389004'
                              | clustering_object$sample_pull == '4820STDY7389005' | clustering_object$sample_pull == '4820STDY7389006'] <- '4820STDY7389003_4820STDY7389004_4820STDY7389005_4820STDY7389006'

clustering_object$sample_pull[clustering_object$sample_pull == '4820STDY7389007' | clustering_object$sample_pull == '4820STDY7389008'
                              | clustering_object$sample_pull == '4820STDY7389009' | clustering_object$sample_pull == '4820STDY7389010'] <- '4820STDY7389007_4820STDY7389008_4820STDY7389009_4820STDY7389010'

clustering_object$sample_pull[clustering_object$sample_pull == '4820STDY7389011' | clustering_object$sample_pull == '4820STDY7389012'
                              | clustering_object$sample_pull == '4820STDY7389013' | clustering_object$sample_pull == '4820STDY7389014'] <- '4820STDY7389011_4820STDY7389012_4820STDY7389013_4820STDY7389014'


clustering_object$sample_pull[clustering_object$sample_pull == 'SKN8104899' | clustering_object$sample_pull == 'SKN8104900'
                              | clustering_object$sample_pull == 'SKN8104901' | clustering_object$sample_pull == 'SKN8104902'] <- 'SKN8104899_SKN8104900_SKN8104901_SKN8104902'

clustering_object$sample_pull[clustering_object$sample_pull == 'SKN8104894' | clustering_object$sample_pull == 'SKN8104895'
                              | clustering_object$sample_pull == 'SKN8104896' | clustering_object$sample_pull == 'SKN8104897'] <- 'SKN8104894_SKN8104895_SKN8104896_SKN8104897'

clustering_object$sample_pull[clustering_object$sample_pull == 'SKN8105197' | clustering_object$sample_pull == 'SKN8105198'
                              | clustering_object$sample_pull == 'SKN8105199' | clustering_object$sample_pull == 'SKN8105200'] <- 'SKN8105197_SKN8105198_SKN8105199_SKN8105200'

clustering_object$sample_pull[clustering_object$sample_pull == 'SKN8105192' | clustering_object$sample_pull == 'SKN8105193'
                              | clustering_object$sample_pull == 'SKN8105194' | clustering_object$sample_pull == 'SKN8105195'] <- 'SKN8105192_SKN8105193_SKN8105194_SKN8105195'


unique(clustering_object$sample_pull)

wf <- run_function_dist(function(obj){return(obj)},
                        obj=clustering_object) 
wf$workflow_id # wf-ee6d57f696

sc_fit <- as_sc_pipeline_run(
  data = "wf-ee6d57f696", 
  annotate = list(
    singler = NULL,
    scanvi = NULL,
    deepsort = NULL
  )
)

wf_export <- sc_api_save_data(dataset_name = "p00_sc_atopic_skin",
                              data = sc_fit #, image = IMAGE
)
# wf-a8e3e424af