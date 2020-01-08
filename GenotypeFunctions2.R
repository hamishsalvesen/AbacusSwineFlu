edit_genes <- function(animals, edit_prop = 0.0, edit_type, edit_trials, embryo_trials, fail_rate, death_rate, sel_prop) { 
  # gene editing workflow
  #     0.0. Subset individuals that would have to be edited
  # Check and see if we have recessives to edit.
  ToEdit <- filter(popdata, genoA == "AA" | "Aa" | "aA")
  Others <- filter (popdata , genoA == "aa")
  
  #     0.1. Sort individuals based on merit to identify the top animals to edit
  SPFpopEdit <- arrange(SPFpop, desc(merit)) ## selection criteria for pigs to be edited. May need to be applied to males or females
  
  # For each recessive to be edited:
  #     1. Select the top edit_prop proportion of animals, currently at least 1 animal always will be edited
  #     2. Do the edit for PP, Pp and pP genotypes
  #     3. Check to see if the edit succeeded
  #     4. Update the animal's genotype
  #     5. Update the edit_status list
  #     6. Update the ge_trials (number of tries until successful edit)
  
  # check that homozygotes do not already supply all the required animals for selection
  if (length(homozygotes$ID) < ceiling(length(animals$ID) * sel_prop)) {
    n_edit <- ceiling(length(animals$ID) * edit_prop)
    if (n_edit >= 1) {
      cat("Attempting to edit ", n_edit, " animals edited with edit_prop = ", edit_prop, "and edit_type = ", edit_type)
    }
    ## 0.1. Sort individuals based on merit ##
    others <- others[order(others$merit, decreasing = TRUE),]
    # add columns to allow for recording editing information
    others$edit_status <- 0
    others$ge_trials <- others$et_trials <- NA
    
    ## 1. Select the top edit_prop proportion of animals.
    edits <- 0
    for (animal in 1:length(others$ID)) {
      if(edits < n_edit){
        # assuming that both heterozygotes and undesirable homozygotes need to be edited (and desirable homozygotes do not):
        # if (others$geno[animal] == "PP" | others$geno[animal] == "Pp" | animals$geno[animal] == "pP") {  - Don't need this anymore, as homozygotes are already removed
        # 2. Do the edit for PP and Pp genotypes
        if (edit_trials > 0) {
          edits <- edits+1
          # 3a. (i) if edit_trials > 0 then only a fixed number of trials will be carried out. 
          # If there is no success before the final trial then the editing process fails.
          outcomes <- rbinom(edit_trials, 1, (1 - fail_rate))
          # check if edit was successful
          if (1 %in% outcomes) {
            # 4. Update the animal's genotype - for now assume that we can edit both in the same editing step
            #animals$geno[animal] <- ifelse((animals$geno[animal] == "Pp"), "PP", "pp")
            others$geno[animal] <-  "pp"
            # 5. Update the edit_status 
            others$edit_status[animal] <- 1
            # 6. Update the animal's edit count with the time of the first successful edit
            others$ge_trials[animal] <- which(outcomes != 0)[1]
            
            # 3b. Was the successfully edited embryo carried to term?
            if (embryo_trials > 0){
              # 3b. (i) If the embryo died then we need to update the fate. If edit_trials > 0 then only a fixed number of trials
              #         will be carried out. If there is no success before the final trial then the editing process
              #         fails.
              outcomes <- rbinom(edit_trials, 1, (1 - death_rate))
              if (any(outcomes)){
                others$et_trials[animal] <- which(outcomes != 0)[1]
              } else {
                # If ET fails, animal is dead/not available for selection
                others$et_trials[animal] <- embryo_trials
                others$fate[animal] <- 0
              } 
            }
          } else {
            others$edit_status[animal] <- -1
          }
        }
        else { print ("edit_trials should never be 0, skipping editing step!") }
        #}
      }
    }
    # print how many animals were edited AND survived the ET
    #print(paste(sum(others$edit_status == 1 & others$fate == 1), "animals successfully edited", sep = " "))
    # record the number of edits attempted, and the number of successful edits to an output list
    edits <- sum(others$edit_status == 1) + sum(others$edit_status == -1)
    success <- sum(others$edit_status == 1 & others$fate == 1)
    animals <- rbind.fill(homozygotes, others)
    return(list(animals=animals, edits = edits, success = success))
  }
  else { cat("Enough resistant animals in the population for selection, no individuals edited")
    # If we have enough homozygotes available for to meet the selection 
    edits <- 0
    success <- 0
    animals$edit_status <- NA
    animals$et_trials <- NA
    animals$ge_trials <- NA
    animals <- animals
    return(list(animals = animals, edits = edits, success = success))
  }
}