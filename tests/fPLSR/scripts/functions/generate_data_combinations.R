generate_data_combinations <- function(N, 
                                       test.options, test.options_names,
                                       options.global,
                                       data.directory,
                                       FORCE_GENERATION = FALSE) {
  
  # Generating all possible combinations
  test.options <- data.frame(as.matrix(expand.grid(test.options)))
  test.options_names <- as.matrix(expand.grid(test.options_names))
  
  # Generating data
  for(i in 1:dim(test.options)[1]){
    
    # Names
    test.options$test.name[i] <- paste(test.options_names[i,], collapse = "_")
    test.options$data.name[i] <- paste("data_", test.options$test.name[i], ".RData", sep = "")
    test.options$data.path[i] <- paste(data.directory, test.options$data.name[i], sep = "")
    
    cat(paste("\n", test.options$test.name[i], ": ", sep = ''))
    
    # Options
    if(!file.exists(test.options$data.path[i]) | FORCE_GENERATION){
      generate_data(N, test.options$data.path[i], test.options[i,], options.global)
      cat("Generated!\n")
    } else{
      cat("It already exists")
    }
    
  }
  
  return(test.options)
  
}

save(generate_data_combinations, file = "scripts/functions/generate_data_combinations.RData")