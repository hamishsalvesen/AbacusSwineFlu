

               ####### AbacusBio 2018/19 ########
### Simulating IAV the economic/geneitc cost of IAV Gene Editing in Pigs #####

               
### LIBRARIES ###               
library (plyr)
library (dplyr)
library (purrr)
library (tidyverse)

################## PARAMETERS/VALUES #########################
AgeDist <- c(0.092,0.055,0.07,0.09,0.115,0.147,0.189,0.242) ###
indexSD <- 10
littersize <- 12 ### Increase within functions for simulated 'embryo flushing'
AgeFirstMate <- 8
FarrowInt <- 5
rem <- AgeFirstMate - FarrowInt #### rem = 3. Use %% to select for breeding every 5 months ######

###############################

BasePop <- function(n, indexSD, AgeDist) {
  NucleusA <- setNames(data.frame(matrix(nrow = n, ncol =11)), c("ID", "gen", "herd", "sex", "merit", "genoA", "genoB", "age", "fate", "IAV", "EditTrials"))
  
  NucleusA$ID <- rep(1:n)
  NucleusA$gen <- 0
  NucleusA$herd <- sample(c("A", "B", "T"), size = n, replace = TRUE, prob = c(0.5, 0.25, 0.25))
  NucleusA$sex <- sample (c("M", "F"), n, replace = TRUE, prob = c(0.3, 0.7))
  NucleusA$merit <- rnorm(n, mean = 0, sd = indexSD) %>% round(digits = 0) 
  NucleusA$genoA <- "A/A"
  NucleusA$genoB <- "B/B"
  NucleusA$age <- sample(rep(1:length(AgeDist) + 8), size = n, replace = TRUE, prob = AgeDist) # +8 so that piglets are ready for breeding more quickly ##3
  NucleusA$fate <- 1
  NucleusA$IAV <- ifelse(NucleusA$genoA == "a/a", "0", "1") ### 1 = susceptible. Therefore $10 cost
  NucleusA$EditTrials <- 0
  
  
  return (NucleusA)
}


CreatePiglets <- function (sires, dams, indexSD, genNo, label, littersize){
  nP <- length(dams$ID)*littersize ###
  piglets <- setNames(data.frame(matrix(nrow = nP, ncol =13)), c("ID", "gen", "herd", "sex", "merit", "genoA", "genoB", "age", "fate", "sire", "dam", "IAV", "EditTrials"))
  
  siresOrdered <- sires %>% filter(age >= AgeFirstMate) %>% arrange(desc(merit))
  siresOrdered2 <- sample(1:length(siresOrdered$ID), length(dams$ID), replace = TRUE)
  siresOrdered2 <- sires[c(siresOrdered2), ]  ####### Creates a list of sires and repeats for the length of the dams vector. Hard selection on males causes repeat use
  
  # split the sire and dam genotypes into the separate alleles - A 
  sgA <- siresOrdered2$genoA
  s1A <- substr(sgA, 1, 1)
  s2A <- substr(sgA, 3, 3)
  dgA <- dams$genoA
  d1A <- substr(dgA, 1, 1)
  d2A <- substr(dgA, 3, 3)
  # Sample which allele of the PRRS resistant genotype is inherited from each parent
  sireAlleleA <- sample(c(1,2),nP,replace = TRUE)
  damAlleleA <- sample(c(1,2),nP,replace = TRUE)
  
  # split the sire and dam genotypes into the separate alleles - B
  sgB <- siresOrdered2$genoB
  s1B <- substr(sgB, 1, 1)
  s2B <- substr(sgB, 3, 3)
  dgB <- dams$genoB
  d1B <- substr(dgB, 1, 1)
  d2B <- substr(dgB, 3, 3)
  # Sample which allele of the PRRS resistant genotype is inherited from each parent
  sireAlleleB <- sample(c(1,2),nP,replace = TRUE)
  damAlleleB <- sample(c(1,2),nP,replace = TRUE)
  
  piglets$ID <- paste0(label, seq(1:nP)) ## label in function for identifier ###
  piglets$gen <- genNo ##generation associated to loop
  piglets$herd <- siresOrdered2$herd ## Male herd can be interpreted as A, AB, ABT lineage. 
  piglets$sex <- sample (c("M", "F"), nP, replace = TRUE, prob = c(0.5, 0.5))
  piglets$merit <- 0.5 * (siresOrdered2$merit + dams$merit) + rnorm(nP, 0, indexSD) %>% round(digits = 1) ## average merit of sire and dam with random variation. 
  piglets$genoA <- paste(ifelse(sireAlleleA==1,s1A,s2A),ifelse(damAlleleA==1,d1A,d2A), sep = "/") ##pasting the alleles selected from boar and dam
  piglets$genoB <- paste(ifelse(sireAlleleB==1,s1B,s2B),ifelse(damAlleleB==1,d1B,d2B), sep = "/")
  piglets$age <- -3 ## -3 fo creation of piglet and 3 months gestation
  piglets$fate <- 1 
  piglets$sire <- siresOrdered2$ID 
  piglets$dam <- dams$ID
  piglets$IAV <- ifelse(piglets$genoA == "a/a", "0", "1") ## $10 is cost of IAV/pig
  piglets$EditTrials <- 0
  
  return(piglets)
}

### AGEING FUNCTION ###

ageing <- function (popdata){
  popdata$age <- popdata$age + 1
  
  alive_pigs <- sum(nrow(popdata))
  
  if (nrow(popdata[popdata$sex == "M" & popdata$age > 38, ]) > 0) {popdata[popdata$sex == "M" & popdata$age > 38,]$fate <- 0} 
  if (nrow(popdata[popdata$sex == "F" & popdata$age > 42, ]) > 0) {popdata[popdata$sex == "F" & popdata$age > 42,]$fate <- 0} ### 42 for females to allow for weaing after piglets creation at 38 months
 
  popdata_mort <- popdata %>% filter(age > 1)     ### No killing off piglets before they are weaned # time to allow for genomic assesment after birth
  popdata_mort$mort <- runif(nrow(popdata_mort)) ### uniform distribution between 1 and 0 randomly assigned
  popdata_mort <- top_frac(popdata_mort, 0.975, mort) ### removes 2.5% of population over the age of 1 
  popdata_mort$mort <- NULL
  popdata_piglets <- popdata %>% filter(age <= 1) 
  popdata <- rbind (popdata_piglets, popdata_mort) ## rebind piglets and undead pigs again 
  
  popdataSows <- filter(popdata, sex == "F" & age >= AgeFirstMate)
  popdataSows$gest <- ((popdataSows$age - AgeFirstMate) / FarrowInt)+1 ### splits into gestation groups
  popdataSows$gest <- floor(popdataSows$gest)
  popdataSows <- split(popdataSows, popdataSows$gest) ## forms lists with the differetn gestations
 
  popdataSows <- lapply(popdataSows, function(x) top_frac(x, 0.5, merit))
  popdataSows <- unlist(data.frame(popdataSows))
  # lapply, applys to each element of the list. unlist function. 
  
  ## does not show culled pigs for looking at just removes the ones dieing
  
  popdata <- filter(popdata, fate == 1) ## only alive pigs retained

  dead_pigs <- alive_pigs - sum(nrow(popdata)) ### is this counting correctly??
  
  print(paste0 ("Died:", dead_pigs))

  return (popdata)
}

        ##### MI gene editing function ###########

edit_genes <- function(popdata, Edit_Efficiency, Embryo_Survival) { 

  A_ToEdit <- filter(popdata, genoA %in% c("A/A", "A/a", "a/A"))
  A_Others <- filter (popdata , genoA == "a/a")
  
  if (nrow(A_ToEdit) >= 1){
    
    A_ToEdit <- separate(A_ToEdit, genoA, into = c("A_Allele1", "A_Allele2"), sep = "/", remove = TRUE) %>% mutate(A_Allele1 = as.character(A_Allele1), A_Allele2 = as.character(A_Allele2))
    
    #  SELECTION OF PIGS ON MERIT ## would need to be on sires and dams the selection, not on the progeny
    # A_ToEdit <- arrange(A_ToEdit, desc(merit)) #%>% top_frac(0.5) # filter out top_frac if only certain animals are to be edited? filter sires/dams, can't select on piglets as we don't know merit!?

    A_ToEditA1 <- A_ToEdit  %>% filter(A_Allele1 == "A") %>% sample_frac(Edit_Efficiency, replace = FALSE) ## Split columns as alleles are edited separately
    A_ToEditA1$A_Allele1 <- tolower(A_ToEditA1$A_Allele1) 
    A_ToEditA1$A_Allele2 <- NULL
    A_ToEditA2 <- A_ToEdit  %>% filter(A_Allele2 == "A") %>% sample_frac(Edit_Efficiency, replace = FALSE) 
    A_ToEditA2$A_Allele2 <- tolower(A_ToEditA2$A_Allele2) 
    A_ToEditA2$A_Allele1 <- NULL
    
    A_ToEditJoinA1 <- A_ToEditA1 %>% dplyr::select(ID, A_Allele1) ### ID as identifier and allele only selected from dataframes
    A_ToEditJoinA2 <- A_ToEditA2 %>% dplyr::select(ID, A_Allele2)
    
    A_ToEdit$A_Allele1 <- NULL 
    A_ToEdit$A_Allele2 <- NULL
    
    ## Bind Alleles back together by ID match. replace NA (non-succesfully edited alleles) with A. ### Both alleles put back together in a single column
    A_ToEdit <- A_ToEdit %>% left_join (A_ToEditJoinA1, by = "ID") %>% left_join(A_ToEditJoinA2, by = "ID") 
    A_ToEdit [c("A_Allele1", "A_Allele2")][is.na(A_ToEdit[c("A_Allele1", "A_Allele2")])] <- "A" 
    A_ToEdit <- A_ToEdit %>% unite(col = genoA, c(A_Allele1, A_Allele2), sep = "/", remove = TRUE) 

    popdata <- rbind(A_ToEdit, A_Others) ### resistant piglets are remerged with edit attempt piglets
   } ## end of if function for A ########
  
  ###### B editing #######
  
  B_ToEdit <- filter(popdata, genoB %in% c("B/B", "B/b", "b/B"))
  B_Others <- filter (popdata , genoB == "b/b")
  
  if (nrow(B_ToEdit) >= 1){
    
    B_ToEdit <- separate(B_ToEdit, genoB, into = c("B_Allele1", "B_Allele2"), sep = "/", remove = TRUE) %>% mutate(B_Allele1 = as.character(B_Allele1), B_Allele2 = as.character(B_Allele2))
    
    #  SELECTION OF PIGS ON MERIT ## would need to be on sires and dams the selection, not on the progeny
    # A_ToEdit <- arrange(A_ToEdit, desc(merit)) #%>% top_frac(0.5) # filter out top_frac if only certain animals are to be edited? filter sires/dams, can't select on piglets as we don't know merit!?
    
    B_ToEditB1 <- B_ToEdit  %>% filter(B_Allele1 == "B") %>% sample_frac(Edit_Efficiency, replace = FALSE) ## Split columns as alleles are edited separately
    B_ToEditB1$B_Allele1 <- tolower(B_ToEditB1$B_Allele1) 
    B_ToEditB1$B_Allele2 <- NULL
    B_ToEditB2 <- B_ToEdit  %>% filter(B_Allele2 == "B") %>% sample_frac(Edit_Efficiency, replace = FALSE) 
    B_ToEditB2$B_Allele2 <- tolower(B_ToEditB2$B_Allele2) 
    B_ToEditB2$B_Allele1 <- NULL
    
    B_ToEditJoinB1 <- B_ToEditB1 %>% dplyr::select(ID, B_Allele1) ### ID as identifier and allele only selected from dataframes
    B_ToEditJoinB2 <- B_ToEditB2 %>% dplyr::select(ID, B_Allele2)
    
    B_ToEdit$B_Allele1 <- NULL 
    B_ToEdit$B_Allele2 <- NULL
    
    ## Bind Alleles back together by ID match. replace NA (non-succesfully edited alleles) with A. ### Both alleles put back together in a single column
    B_ToEdit <- B_ToEdit %>% left_join (B_ToEditJoinB1, by = "ID") %>% left_join(B_ToEditJoinB2, by = "ID") 
    B_ToEdit [c("B_Allele1", "B_Allele2")][is.na(B_ToEdit[c("B_Allele1", "B_Allele2")])] <- "B" 
    B_ToEdit <- B_ToEdit %>% unite(col = genoB, c(B_Allele1, B_Allele2), sep = "/", remove = TRUE) 
    
 
    popdata <- rbind(B_ToEdit, B_Others) ### resistant piglets are remerged with edit attempt piglets
  } ## end of if function  
  
  if (nrow(A_ToEdit) >=1 | nrow(B_ToEdit) >=1) {
  popdata$EditTrials <- 1
  ### add 1 to EditTrials when a piglet (zygote) undergoes editing #same cost independent of number of allele targets
  
  popdata$mort <- runif(nrow(popdata)) ### New column. Random, even 0-1 distribution to all piglets
  popdata <- top_frac(popdata, Embryo_Survival, mort) ### selects top 60% of mort column in popdata. Set in function for embryo survival 
  popdata$mort <- NULL 
  }
  
  return (popdata)  } 


########################################## SPF POPULATION CREATION ################################

SPFpop <- BasePop(9000, indexSD, AgeDist)  #### creates base population of piggys 

  ProdPop <- data.frame(matrix(ncol = 13, nrow = 0))
  colnames(ProdPop) <- c("ID", "gen", "herd", "sex", "merit", "genoA", "genoB", "age", "fate", "sire", "dam", "IAV", "EditTrials")
  
  Prod_Females <- data.frame(matrix(ncol = 13, nrow = 0))
  colnames(Prod_Females) <- c("ID", "gen", "herd", "sex", "merit", "genoA", "genoB", "age", "fate", "sire", "dam", "IAV", "EditTrials")
  
  MultPop <- data.frame(matrix(ncol = 13, nrow = 0))
  colnames(MultPop) <- c("ID", "gen", "herd", "sex", "merit", "genoA", "genoB", "age", "fate", "sire", "dam", "IAV", "EditTrials")
  
  Mult_Females <- data.frame(matrix(ncol = 13, nrow = 0))
  colnames(Mult_Females) <- c("ID", "gen", "herd", "sex", "merit", "genoA", "genoB", "age", "fate", "sire", "dam", "IAV", "EditTrials")
  
  BWpop <- data.frame(matrix(ncol = 13, nrow = 0))
  colnames(BWpop) <- c("ID", "gen", "herd", "sex", "merit", "genoA", "genoB", "age", "fate", "sire", "dam", "IAV", "EditTrials")
  
  CommercialPop <- data.frame(matrix(ncol = 13, nrow = 0))
  colnames(BWpop) <- c("ID", "gen", "herd", "sex", "merit", "genoA", "genoB", "age", "fate", "sire", "dam", "IAV", "EditTrials")
  
  for (g in -60:0){

      print(paste0("Gen:", g, sep = " "))
    
    ### selecting 85% of gilts to conceive. Need to somehow be in next cycle for breeding and not removed from breeders if didn't conceive?
    #SPF_Gilts <- SPFpop %>% filter(age = AgeFirstMate) 
    #SPF_Gilts$Conceived <- runif(nrow(SPF_Gilts))   ### uniform distribution between 1 and 0 randomly assigned
    #SPF_Gilts <- top_frac(SPF_Gilts, 0.85, Conceived) ### 
    #SPF_Gilts$Conceived <- NULL ### remove pigs that don't conceive from herd? culling measure often used.... 
    ### remove = sign from other Fem selection. create new SPF$Conception col. random distribution for 85%, assign to 1 in new col. if 1 then add to fem breed 
    #BreedSPFA_Fem <- rbind.fill(BreedSPFA_Fem, SPF_Gilts) %>% top_frac(0.25, merit)
      
    BreedSPFA_Male <- SPFpop %>% filter (sex == "M" & herd == "A") %>% filter(age >= AgeFirstMate) %>% top_frac(0.05, merit) 
    BreedSPFA_Fem <- filter(SPFpop, sex == "F" & herd == "A") %>% filter(age >= AgeFirstMate & age %% FarrowInt == rem)  %>% top_frac(0.25, merit) ### Kept from prod tier. Can't breed twice!
    
    newGenSPFNucA <- CreatePiglets(BreedSPFA_Male, BreedSPFA_Fem, indexSD, g, paste0("SPFA",g,"_"), littersize)
    
    SPFNucA_Fem <- filter(newGenSPFNucA, sex == "F")  %>% top_frac(0.2, merit) 
    SPFNucA_Male <- filter(newGenSPFNucA, sex == "M") %>% top_frac (0.05, merit)
    
    ##### separate into B herd ###########
    BreedSPFB_Male <- SPFpop %>% filter (sex == "M" & herd == "B") %>% filter(age >= AgeFirstMate) %>% top_frac(0.05, merit) 
    BreedSPFB_Fem <- filter(SPFpop, sex == "F" & herd == "B") %>% filter(age >= AgeFirstMate & age %% FarrowInt == rem) %>% top_frac(0.25, merit) ### remove these pigs from prod tier. Will leave some pigs that can be bred from!
    
    newGenSPFNucB <- CreatePiglets(BreedSPFB_Male, BreedSPFB_Fem, indexSD, g, paste0("SPFB",g,"_"), littersize)
    
    SPFNucB_Fem <- filter(newGenSPFNucB, sex == "F")  %>% top_frac(0.2, merit) #### maybe need ageDist to make sure pigs of all ages are selected ### 
    SPFNucB_Male <- filter(newGenSPFNucB, sex == "M") %>% top_frac (0.05, merit)
    
    ##### separate into T herd ###########
    BreedSPFT_Male <- SPFpop %>% filter (sex == "M" & herd == "T") %>% filter(age >= AgeFirstMate) %>% top_frac(0.05, merit) 
    BreedSPFT_Fem <- filter(SPFpop, sex == "F" & herd == "T") %>% filter(age >= AgeFirstMate & age %% FarrowInt == rem ) %>% top_frac(0.25, merit) ### remove these pigs from prod tier. Will leave some pigs that can be bred from!
    
    newGenSPFNucT <- CreatePiglets(BreedSPFT_Male, BreedSPFT_Fem, indexSD, g, paste0("SPFT",g,"_"), littersize)

    SPFNucT_Fem <- filter(newGenSPFNucT, sex == "F")  %>% top_frac(0.2, merit) #### maybe need ageDist to make sure pigs of all ages are selected ### 
    SPFNucT_Male <- filter(newGenSPFNucT, sex == "M") %>% top_frac (0.05, merit)
    
    ### new SPF pop piglets rbind at the end of the code
    
    #########
    
    SPF_Breeders <- rbind.fill(BreedSPFA_Fem, BreedSPFB_Fem, BreedSPFT_Fem) ### animals for SPF breeding. Make table to exclude from ProdPop breeding cycle
    ### will rejoin SPFpop at end of code incase they are no longer in top % for mating and can be used in below tiers  
    ### males are not removed as they can inseminate multiple tiers through an AI program
    
    
    if (g == -60) { 
      SPFpopAll <- data.frame(SPFpop)
    } else  (SPFpopAll <- rbind.fill(SPFpopAll, SPFpop))
    
    
if (g >= -50){ 
  
  SPFpop <-  anti_join(SPFpop, SPF_Breeders, by = "ID") 

  SPF_ProdPop_Fem <- SPFpop %>% filter(sex == "F" & herd == "A") %>% filter(age >= AgeFirstMate & age %% FarrowInt == rem) %>% top_frac(0.25, merit) 
  SPFpop <-  anti_join(SPFpop, SPF_ProdPop_Fem, by = "ID") ## SPF pop loses the animals transferred to prod pop. Not reintegrated and remain in ProdPop until culled 
  
  SPF_ProdPop_Males <- SPFpop %>% filter(sex == "M" & herd == "A") %>% filter(age >= AgeFirstMate) %>% top_frac(0.1, merit) 
  #SPFpop <-  anti_join(SPFpop, SPF_ProdPop_Males, by = "ID") ## removes prod pop males from the SPFpop ### don't need as AI is performed on SPF pops

  ProdPop <- rbind.fill(ProdPop, SPF_ProdPop_Fem) 
  
  Prod_Males <- rbind.fill(ProdPop, SPF_ProdPop_Males) %>% filter(sex == "M") %>% filter(age >= AgeFirstMate) %>% top_frac(0.15, merit) # selects for top males from either ProdPop and SPF Nucleus ## AI for breeding 
  Prod_Females <- ProdPop %>% filter(sex == "F") %>% filter(age >= AgeFirstMate & age %% FarrowInt == rem) %>% top_frac(0.5, merit) %>% mutate(ID = as.character(ID))### already over age first mate, ensure they are correct farrowing interval. 
  ### coerce ID to be a character vector for multpop antijoin below
  
  NewProdPop <- CreatePiglets(Prod_Males, Prod_Females, indexSD, g, paste0("Prod",g,"_"), littersize)
  
  NewProd_Females <- NewProdPop %>% filter(sex == "F") %>% top_frac(0.2, merit) 
  NewProd_Males <- NewProdPop %>% filter(sex == "M") %>% top_frac(0.05, merit) ### Retained males and female piglets for ProdPop
  
  ProdPop <- rbind.fill(ProdPop, NewProd_Females, NewProd_Males) ### combine all animals to stay in ProdPop
  
  if (g >= -50) { 
    ProdPopAll <- data.frame(ProdPop)
  } else  (ProdPopAll <- rbind.fill(ProdPopAll, ProdPop))
  
  }
  
if (g >= -35){
  
  ### A or B females as can be bred within mult tier or from prod?? ## no transfer of Mult pigs back into breeding data frame so only A females used, piglets will be B herd. 
  ### Need to make sure Prod_Females are put back into ProdPop at the end of generation! 
  
  ProdPop <- anti_join(ProdPop, Prod_Females, by = "ID") ## removes pigs bred in ProdPop from being bred again in same generation 
  
  Prod_MultPop_Females <- ProdPop %>% filter(sex == "F" & herd == "A" | herd == "B") %>% filter(age >= AgeFirstMate & age %% FarrowInt == rem) %>% top_frac(0.5, merit)  
  ProdPop <- anti_join(ProdPop, Prod_MultPop_Females, by = "ID") ## removes pigs transferred down to multpop from ProdPop permanently
  
  SPF_MultPop_Males <- SPFpop %>% filter(sex == "M" & herd == "B") %>% filter(age >= AgeFirstMate) %>% top_frac(0.1, merit) 
  SPFpop <-  anti_join(SPFpop, SPF_MultPop_Males, by = "ID") ## removes MultPop boars from SPF pop 
  
  MultPop <- rbind.fill(MultPop, Prod_MultPop_Females, SPF_MultPop_Males) # All pigs used in MultPop are now permanently there and not in other tiers
  
  Mult_Males <- MultPop %>% filter(sex == "M") %>% filter(age >= AgeFirstMate)
  Mult_Females <- MultPop %>% filter(sex == "F") %>% filter (age >= AgeFirstMate & age %% FarrowInt == rem) %>% top_frac(0.5, merit)
  
  NewMultPop <- CreatePiglets(Mult_Males, Mult_Females, indexSD, g, paste0("Mult",g,"_"), littersize)
  
  NewMult_Females <- NewMultPop %>% filter(sex == "F") %>% top_frac(0.25, merit) 
  NewMult_Males <- NewMultPop %>% filter(sex == "M") %>% top_frac(0.05, merit) ### retain piglets in MultPop for flowing down.
  
  MultPop <- rbind.fill(MultPop, NewMult_Females, NewMult_Males) %>% mutate(ID = as.character(ID))
  ### coerce ID to be a character vector for BW antijoin below 

  if (g >= -35) { 
    MultPopAll <- data.frame(MultPop)
  } else  (MultPopAll <- rbind.fill(MultPopAll, MultPop))
  
  } 
    
    ###############################
  
  if (g >= -25){
    
    MultPop <-  anti_join(MultPop, Mult_Females, by = "ID") #### keep MultPop females that have bred out of BW pop
    
    Mult_BW_Fem <- MultPop %>% filter(sex == "F" & herd == "B" | herd == "A") %>% filter(age >= AgeFirstMate & age %% FarrowInt == rem) %>% top_frac(0.7, merit) 
    MultPop <- anti_join(MultPop, Mult_BW_Fem, by = "ID") ## removes pigs transferred down to BWpop from MultPop permanently

    SPF_BWpop_Males <- SPFpop %>% filter(sex == "M" & herd == "T") %>% filter(age >= AgeFirstMate) %>% top_frac(0.1, merit)
    SPFpop <-  anti_join(SPFpop, SPF_MultPop_Males, by = "ID") ## no returning pigs as not an AI program
    
    BWpop <- rbind.fill(BWpop, Mult_BW_Fem, SPF_BWpop_Males) 
    
    BW_Fem_Breed <- BWpop %>% filter(sex == "F") %>% filter(age >= AgeFirstMate & age %% FarrowInt == rem) 
    BW_Males <- BWpop %>% filter(sex == "M") %>% filter(age >= AgeFirstMate)

    BW_piglets <- CreatePiglets(SPF_BWpop_Males, BW_Fem_Breed, indexSD, g, paste0("BW",g,"_"), littersize) 

    NewBW_Females <- BW_piglets %>% filter(sex == "F") %>% top_frac(0.1, merit) 

    BWpop <- rbind.fill(BWpop, NewBW_Females)
 
    if (g >= -25) { 
      BWpopAll <- data.frame(BWpop)
    } else  (BWpopAll <- rbind.fill(BWpopAll, BWpop))
    
     }

    
#### POPULATION CONTROLS ######

    SPFpop <- rbind.fill(SPFpop, SPF_Breeders, SPFNucA_Male, SPFNucA_Fem, SPFNucB_Male, SPFNucB_Fem, SPFNucT_Male, SPFNucT_Fem) %>% distinct() 
    ## All SPF pigs. SPF pigs moved to another tier are already excluded by antijoins
    ## put breeders back with SPFpop, if they are not used for SPF pop breeding next round they may be used for PN. Need to select breeders from same pool as before ##
    
    if (sum(SPFpop$age > 8) > 9000) {
      
      NextBreeders <- SPFpop %>% filter (sex =="F" & age > 9 & age %% FarrowInt == 2) ### filters all pigs that will breed in the next loop to be included for selection 
      
      SPFpopA_females <- SPFpop %>% filter(sex == "F" & herd == "A") %>% filter(age >= 8) %>% top_n(4500, merit)
      SPFpopB_females <- SPFpop %>% filter(sex == "F" & herd == "B") %>% filter(age >= 8) %>% top_n(2250, merit) 
      SPFpopT_females <- SPFpop %>% filter(sex == "F" & herd == "T") %>% filter(age >= 8) %>% top_n(2250, merit)
      
      SPFpop_males <- SPFpop %>% filter(sex == "M") %>% filter(age >= 8) %>% top_n(500, merit) ### changed for forward simm to maintain male herd distribution 
      #SPFpopB_males <- SPFpop %>% filter(sex == "M" & herd == "B") %>% filter(age >= 8) %>% top_n(500, merit) 
      #SPFpopT_males <- SPFpop %>% filter(sex == "M", herd == "T") %>% filter(age >= 8) %>% top_n(500, merit) 
      
      #SPFpopA_piglets <- SPFpop %>% filter(herd == "A" & sex =="M") %>% filter(age < 4 & age >= 1) %>% top_n(1000, merit) 
      #SPFpopB_piglets <- SPFpop %>% filter(herd == "B" & sex =="M") %>% filter(age < 4 & age >= 1) %>% top_n(1000, merit) 
      #SPFpopT_piglets <- SPFpop %>% filter(herd == "T" & sex =="M") %>% filter(age < 4 & age >= 1) %>% top_n(1000, merit)
      
      topSPFpop <-  rbind.fill(NextBreeders, SPFpopA_females, SPFpopB_females, SPFpopT_females, SPFpop_males)

      SPFpop <- semi_join(topSPFpop, SPFpop, by = "ID") %>% distinct() %>% mutate(ID = as.character(ID))### joins only the selected SPFpop number and removes any pigs duplicated. 
      
      SPFpopA_females <- NULL ### remove tables to prevent accumulation of data slowing code
      SPFpopB_females <- NULL
      SPFpopT_females <- NULL
      SPFpopA_males <- NULL
      SPFpopB_males <- NULL
      SPFpopT_males <- NULL
      SPFpopA_piglets <- NULL
      SPFpopB_piglets <- NULL
      SPFpopT_piglets <- NULL
      topSPFpop <- NULL
      SPF_Breeders <- NULL 
    }
    
    SPFpop <- data.frame(ageing(SPFpop))  
    
    print(paste0("SPF Males:", sum(SPFpop$sex == "M"), sep = " "))
    print(paste0 ("SPF Females:", sum(SPFpop$sex =="F"), sep = " "))
    
    ########################
  
 if (sum(ProdPop$age > 8) > 33000){
   
   NextBreedersProd <- ProdPop %>% filter (sex =="F" & age > 9 & age %% FarrowInt == 2)
   
   ProdPopFem <- ProdPop %>% filter(sex == "F") %>% filter(age >= AgeFirstMate) %>% top_n(33000, merit) 
  ProdPopMales <- ProdPop %>% filter (sex == "M") %>% filter(age > AgeFirstMate) %>% top_n(800, merit)
  ProdPopPiglets <- ProdPop %>% filter(sex == "F") %>% filter(age < AgeFirstMate & age >= 1) %>% top_n(15000, merit)
  ProdPop <- rbind.fill(ProdPopFem, ProdPopMales, ProdPopPiglets, NextBreedersProd) 
 }
    
  ProdPop <- rbind.fill(ProdPop, Prod_Females) %>% distinct() #
   
  ProdPop <- data.frame(ageing(ProdPop)) 
  
  print(paste0("Prod Males:", sum(ProdPop$sex == "M"), sep = " "))
  print(paste0 ("Prod Females:", sum(ProdPop$sex =="F"), sep = " "))
  
  #############
  
  if (sum(MultPop$age > 8) > 80000) {  
  MultPopFem <- MultPop %>% filter(sex == "F") %>% filter(age >= AgeFirstMate) %>% top_n(80000, merit)
  MultPopMales <- MultPop %>% filter(sex == "M") %>% filter(age >= AgeFirstMate) %>% top_n(10000, merit)
  MultPopPiglets <- MultPop %>% filter(sex == "F") %>% filter(age < AgeFirstMate & age >= 1) %>% top_n(20000, merit)
  MultPop <- rbind.fill(MultPopFem, MultPopMales, MultPopPiglets)
  }
  
  MultPop <- rbind.fill(Mult_Females, MultPop) %>% distinct()

  MultPop <- ageing(MultPop)
  
  print(paste0("Mult Males:", sum(MultPop$sex == "M"), sep = " "))
  print(paste0 ("Mult Females:", sum(MultPop$sex =="F"), sep = " "))
  
  if (sum(BWpop$age > 8) > 150000) {  
    
  BWpopFem <- BWpop %>% filter(sex == "F") %>% filter(age >= AgeFirstMate) %>% top_n(100000, merit)
  BWpopMales <- BWpop %>% filter(sex == "M") %>% filter(age >= AgeFirstMate) %>% top_n(10000, merit)
  BWpopPiglets <- BWpop %>% filter(sex == "F") %>% filter(age < AgeFirstMate & age >= 1) %>% top_n(20000, merit) 
  BWpop <- rbind.fill(BWpopFem, BWpopMales, BWpopPiglets) 
  }
  
  BWpop <- ageing(BWpop)
  
  print(paste0("BW Males:", sum(BWpop$sex == "M"), sep = " "))
  print(paste0 ("BW Females:", sum(BWpop$sex =="F"), sep = " "))

  Allpop <- rbind.fill(SPFpop, ProdPop, MultPop, BWpop) 

  if (g == -60) { 
    Allpop2 <- data.frame(Allpop)
  } else  (Allpop2 <- rbind.fill(Allpop2, Allpop))
  
  
  } ### end of g loop

  BurnIn <- list(SPFpop = SPFpop, ProdPop = ProdPop, MultPop = MultPop, BWpop = BWpop)
  
  
  ##### write.csv(BurnInPop, paste0("BurnInPop", Sys.Date(), ".csv")) - file for reiterating usage of the same burn in population
  
 
  


  

 

 
