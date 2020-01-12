edit_genes <- function(animals) { 
  # gene editing workflow
  #     0.0. Subset individuals that would have to be edited
  # Check and see if we have recessives to edit.
  A_ToEdit <- filter(SPFpop, genoA %in% c("A/A", "A/a", "a/A"))
  A_Others <- filter (SPFpop , genoA == "a/a")
  
  ## Change to piglets for function format
  if (nrow(A_ToEdit) > 1){
# Split_A_ToEdit <- strsplit(A_ToEdit$genoA, 2) 
### need to Strspilt, unnest, add seperated columns
  
A_ToEdit <- separate(A_ToEdit, genoA, into = c("Allele1", "Allele2"), sep = "/", remove = TRUE) %>% mutate(Allele1 = as.character(Allele1), Allele2 = as.character(Allele2))
  
#     0.1. Sort individuals based on merit to identify the top animals to edit
A_ToEdit <- arrange(A_ToEdit, desc(merit)) #%>% top_frac(0.5) ## selection criteria for pigs to be edited. May need to be applied to males or females
  ##filter out top_frac if only certain animals are to be edited? filter sires/dams or piglets born?

###need to have 'if' for only editing A alleles and keping a in the table

A_ToEditA1 <- A_ToEdit  %>% filter(Allele1 == "A") ## need to select at the same time as genes are in the same pig vessel
A_ToEditA1$Allele2 <- NULL ### individual alleles are edited
A_ToEditA2 <- A_ToEdit  %>% filter(Allele2 == "A") ## need to select at the same time as genes are in the same pig vessel
A_ToEditA2$Allele1 <- NULL

  ### filter or select pigs eligible for editing ## filter for piglets at -3 then can put in after CreatePiglets
  ##rbind A_others with A_toEdit tables  

  ## Select proportion of piglets for relevant editing efficiency. #Can do separately as alleles are independently edited?
A_ToEditA1 <- sample_frac(A_ToEditA1, 0.7, replace = FALSE)
A_ToEditA2 <- sample_frac(A_ToEditA2, 0.7, replace = FALSE) ## random selection of alleles to edit. Essentially editing efficiency proxy
## need to save the other 30% and remerge below  

A_ToEditA1$Allele1 <- tolower(A_ToEditA1$Allele1) ## use pipe for putting through and avoid splitting table?
A_ToEditA2$Allele2 <- tolower(A_ToEditA2$Allele2) ### apply to and retain in retain in data frame
### need to merge the column of Allele2 data to allele1

## want to superimpose new alleles list on top of the A_ToEdit table
## Add new alleles data over the old alleles in the table. Bind others with edited allele data
  
## Need join for allele1 and allele2 to be in the same table for appropriate pigs
### check that the right alleles are merged back to the relevant piglet!! 
SPFpopEdit <- full_join(A_ToEditA1, A_ToEditA2, by = "ID") 
## want to join allele2 column to the rest of the table. 

SPFpopEdit <- SPFpopEdit %>% unite(col = genoA, c(Allele1, Allele2), sep = "/", remove = TRUE)

} ## end of if function
  
##join for the original population. anti_join?
  
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