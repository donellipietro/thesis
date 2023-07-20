generate_data_combinations <- function(N, 
                                       test.options, test.options_names, 
                                       mesh.directory, mesh.options,
                                       data.directory,
                                       generate_X,
                                       FORCE_GENERATION = FALSE) {
  
  # Generating all possible combinations
  test.options <- data.frame(as.matrix(expand.grid(test.options)))
  test.options_names <- as.matrix(expand.grid(test.options_names))
  
  # Generating data
  for(i in 1:dim(test.options)[1]){
    
    # Names
    test.options$mesh.path[i] <- paste(mesh.directory, test.options$mesh.name[i], "/", sep = "")
    test.options$test.name[i] <- paste(test.options_names[i,], collapse = "_")
    test.options$data.name[i] <- paste("data_", test.options$test.name[i], ".RData", sep = "")
    test.options$data.path[i] <- paste(data.directory, test.options$data.name[i], sep = "")
    
    cat(paste("\n", test.options$test.name[i], ": ", sep = ''))
    
    # Options
    if(!file.exists(test.options$data.path[i]) | FORCE_GENERATION){
      load(mesh.options[[test.options$mesh.name[i]]]$generate_X)
      generate_data(N,
                    test.options$mesh.path[i], test.options$data.path[i],
                    generate_X, test.options[i,],
                    mesh.options[[test.options$mesh.name[i]]]$mesh_finer.path,
                    mesh.options[[test.options$mesh.name[i]]]$mesh.area_refine)
      cat("Generated!")
    } else{
      cat("It already exists")
    }
    
  }
  
  return(test.options)
  
}

save(generate_data_combinations, file = "utils/functions/generate_data_combinations.RData")