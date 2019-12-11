####### functions and running code ########

library (plyr)
library (dplyr)
library (purrr)
library (tidyverse)
library (genetics)


################## PARAMETERS/VALUES #########################
######## setting parameters/values required for code ########## 

AgeDist <- c(0.242,0.189,0.147,0.115,0.09,0.07,0.055,0.092)
indexSD <- 10 
littersize <- 12 ### consider function for changing the littersize and having a fertilisation proportion ###
AgeFirstMate <- 8
FarrowInt <- 5
rem <- AgeFirstMate - FarrowInt #### females to breed. rem = 3 which is selects them for cycles to breed in %% ######

#SPFNucMalesSelect <- 0.02
#SPFNucFemalesSelect <- 0.1 

genotypesA_susceptible <- c("A/A") ## introduce after gene editing
allelesA_susceptible <- c("A")

genotypesB_susceptible <- c("B/B") ## introduce after gene editing
allelesB_susceptible <- c("B")

#### ForwardSim Genotyping needs alleles for resistance ###

genotypesA <- c("A/A", "A/a", "a/A", "aa") ## introduce after gene editing
allelesA <- c("A", "a")

genotypesB <- c("B/B", "B/b", "b/B", "bb") ## introduce after gene editing
allelesB <- c("B", "b")

                            ###############################
                                   ####################

BasePop <- function(n, indexSD, AgeDist) {
  NucleusA <- setNames(data.frame(matrix(nrow = n, ncol =9)), c("ID", "gen", "herd", "sex", "merit", "genoA", "genoB", "age", "fate"))
  
  NucleusA$ID <- rep(1:n)
  NucleusA$gen <- 0
  NucleusA$herd <- sample(c("A", "B", "T"), size = n, replace = TRUE, prob = c(0.5, 0.25, 0.25))
  NucleusA$sex <- sample (c("M", "F"), n, replace = TRUE, prob = c(0.3, 0.7))
  NucleusA$merit <- rnorm(n, mean = 0, sd = indexSD)
  NucleusA$genoA <- sample(expectedGenotypes(as.genotype(genotypesA_susceptible), alleles=allelesA_susceptible, ploidy=2, sort=TRUE,haplotype=FALSE), n, replace = TRUE)
  NucleusA$genoB <- sample(expectedGenotypes(as.genotype(genotypesB_susceptible), alleles=allelesB_susceptible, ploidy=2, sort=TRUE,haplotype=FALSE), n, replace = TRUE)
  NucleusA$age <- sample(rep(1:length(AgeDist)), size = n, replace = TRUE, prob = AgeDist)
  NucleusA$fate <- 1
  
  return (NucleusA)
}

CreatePiglets <- function (sires, dams, indexSD, genNo, label, littersize){
  nP <- length(dams$ID)*littersize ###
  piglets <- setNames(data.frame(matrix(nrow = nP, ncol =11)), c("ID", "gen", "herd", "sex", "merit", "genoA", "genoB", "age", "fate", "sire", "dam"))
  
  siresOrdered <- sires %>% filter(age >= AgeFirstMate) %>% arrange(desc(merit)) #### need to assign sires to the relevant dam 
  siresOrdered2 <- sample(1:length(siresOrdered$ID), length(dams$ID), replace = TRUE)
  siresOrdered2 <- sires[c(siresOrdered2), ]  ############ don't understand this code well #########
  
  piglets$ID <- paste0(label, seq(1:nP)) ## label refers to generation of creation as identifier, doesn't need to be additive ###
  piglets$gen <- genNo 
  piglets$herd <- siresOrdered2$herd ## showing male herd is enough to know lineage? ###
  piglets$sex <- sample (c("M", "F"), nP, replace = TRUE, prob = c(0.5, 0.5))
  piglets$merit <- 0.5 * (siresOrdered2$merit + dams$merit) + rnorm(nP, 0, indexSD) ##### is this variation a realistic meiotic figure? #### are these the actual dam and sire merits? ### how to check...?##
  piglets$genoA <- sample(expectedGenotypes(as.genotype(genotypesA_susceptible), alleles=allelesA_susceptible, ploidy=2, sort=TRUE,haplotype=FALSE), nP, replace = TRUE)
  piglets$genoB <- sample(expectedGenotypes(as.genotype(genotypesB_susceptible), alleles=allelesB_susceptible, ploidy=2, sort=TRUE,haplotype=FALSE), nP, replace = TRUE)
  piglets$age <- -3         
  piglets$fate <- 1 
  piglets$sire <- siresOrdered2$ID ### doesn't appear to always work... ###
  piglets$dam <- dams$ID
  
  return(piglets)
}


ageing <- function (popdata){
  popdata$age <- popdata$age + 1
  
  alive_pigs <- sum(nrow(popdata))
  
  if (nrow(popdata[popdata$sex == "M" & popdata$age > 38, ]) > 0) {popdata[popdata$sex == "M" & popdata$age > 38,]$fate <- 0} ### go to previous code to see why it won't work!!!
  if (nrow(popdata[popdata$sex == "F" & popdata$age > 40, ]) > 0) {popdata[popdata$sex == "F" & popdata$age > 40,]$fate <- 0} ### as to help with understanding ###
 
  popdata_mort <- popdata %>% filter(age > 1)     ### don't want to be killing off piglets before they are weaned ## have littersize variability for realism ###
  popdata_mort$mort <- runif(nrow(popdata_mort)) ### uniform distribution between 1 and 0 randomly assigned
  popdata_mort <- top_frac(popdata_mort, 0.975, mort) ### removes 2.5% of population over the age of 1 
  popdata_mort$mort <- NULL
  popdata_piglets <- popdata %>% filter(age <= 1) ### extract piglets excluded from mortality ###
  popdata <- rbind (popdata_piglets, popdata_mort) ## combine piglets and undead pigs again ### was losing >1 before so population stagnated ###
  
  ## does not show culled pigs for looking at just removes the ones dieing
  
  popdata <- filter(popdata, fate == 1)   # remove rows where popdata fate == 0

  dead_pigs <- alive_pigs - sum(nrow(popdata))
  
  print(paste0 ("Died:", dead_pigs))

  return (popdata)
}




                                          ########################
########################################## SPF POPULATION CREATION ################################


SPFpop <- BasePop(10000, indexSD, AgeDist)   #### creates base population of piggys ### start with ~10,000 to have enough for flowing down ### 

  ProdPop <- data.frame(matrix(ncol = 11, nrow = 0))
  colnames(ProdPop) <- c("ID", "gen", "herd", "sex", "merit", "genoA", "genoB", "age", "fate", "sire", "dam")
  
  MultPop <- data.frame(matrix(ncol = 11, nrow = 0))
  colnames(MultPop) <- c("ID", "gen", "herd", "sex", "merit", "genoA", "genoB", "age", "fate", "sire", "dam")
  
  BWpop <- data.frame(matrix(ncol = 11, nrow = 0))
  colnames(BWpop) <- c("ID", "gen", "herd", "sex", "merit", "genoA", "genoB", "age", "fate", "sire", "dam")
  
  
  for (g in -60:0){
  print(paste0("Gen:", g, sep = " "))
    
  ##### separate into A herd ########### 
    
  BreedSPFA_Male <- filter(SPFpop, sex == "M", herd == "A", age >= AgeFirstMate)
  BreedSPFA_Fem <- filter(SPFpop, sex == "F", herd == "A", age >= AgeFirstMate, age %% FarrowInt == rem)
  
  newGenSPFNucA <- CreatePiglets(BreedSPFA_Male, BreedSPFA_Fem, indexSD, g, paste0("SPF",g,"_"), littersize)
                                   
  SPFNucA_Fem <- filter(newGenSPFNucA, sex == "F")  %>% top_frac(0.1, merit) #### maybe need ageDist to make sure pigs of all ages are selected ### 
  SPFNucA_Male <- filter(newGenSPFNucA, sex == "M") %>% top_frac (0.02, merit)
  
  ##### separate into B herd ###########
  
  BreedSPFB_Male <- filter(SPFpop, sex == "M", herd == "B", age >= AgeFirstMate)
  BreedSPFB_Fem <- filter(SPFpop, sex == "F", herd == "B", age >= AgeFirstMate, age %% FarrowInt == rem)
  
  newGenSPFNucB <- CreatePiglets(BreedSPFB_Male, BreedSPFB_Fem, indexSD, g, paste0("SPF",g,"_"), littersize)
  
  SPFNucB_Fem <- filter(newGenSPFNucB, sex == "F") %>% top_frac(0.1, merit)
  SPFNucB_Male <- filter(newGenSPFNucB, sex == "M") %>% top_frac (0.02, merit)
  
  ##### separate into T herd ###########
  
  BreedSPFT_Male <- filter(SPFpop, sex == "M", herd == "T", age >= AgeFirstMate)
  BreedSPFT_Fem <- filter(SPFpop, sex == "F", herd == "T", age >= AgeFirstMate, age %% FarrowInt == rem)
  
  newGenSPFNucT <- CreatePiglets(BreedSPFT_Male, BreedSPFT_Fem, indexSD, g, paste0("SPF",g,"_"), littersize)
  
  SPFNucT_Fem <- filter(newGenSPFNucT, sex == "F") %>% top_frac(0.1, merit)
  SPFNucT_Male <- filter(newGenSPFNucT, sex == "M") %>% top_frac (0.02, merit)
  
  SPFpop <- rbind.fill(SPFpop, SPFNucA_Fem, SPFNucA_Male, SPFNucB_Fem, SPFNucB_Male, SPFNucT_Fem, SPFNucT_Male) 
  ### put this at the end of the code and the remove breeders will stay out? ###
  
  SPFpop <- data.frame(ageing(SPFpop)) 
  
  
  ### need to remove pigs from the list when transferring to a new population #####
  ### SPF_A pigs that are filtered into the production nucleus still breed in the SPF... Can't be both! ####
 
  print(paste0("SPF Males:", sum(SPFpop$sex == "M"), sep = " "))
  print(paste0 ("SPF Females:", sum(SPFpop$sex =="F"), sep = " "))
  
if (g >= -40){
  
  ## subset function may be 
  SPF_ProdPop_Fem <- SPFpop %>% filter(sex == "F" & herd == "A") %>% top_frac(0.25, merit) #adds top 25% of SPF_A females to herd. leaves in SPF herd for multiple breeding... ### ### Take not top merit but middling aniamls
  SPF_ProdPop_Males <- SPFpop %>% filter(sex == "M" & herd == "A") %>% top_frac(0.1, merit)  #adds top 10% of SPF_A males to herd. Create a new herd and rename both herds to remove theses males
  
  ProdPop <- rbind.fill(ProdPop, SPF_ProdPop_Fem, SPF_ProdPop_Males) ### need ProdPop as some are retained
  
  Prod_Males <- ProdPop %>% filter(sex == "M", age >= AgeFirstMate)
  Prod_Females <- ProdPop %>% filter(sex == "F", age >= AgeFirstMate, age %% FarrowInt == rem)
  
  NewProdPop <- CreatePiglets(Prod_Males, Prod_Females, indexSD, g, paste0("Prod",g,"_"), littersize)
  
  AddProd_Females <- NewProdPop %>% filter(sex == "F") %>% top_frac(0.2, merit) ## takes the top 20% of new females created, all that is required. ####
  
  ProdPop <- rbind.fill(ProdPop, AddProd_Females) ### may be confusing with adding the above rows to data frame in wrong order? ###
  
  ProdPop <- data.frame(ageing(ProdPop)) 
  
  print(paste0("Prod Males:", sum(ProdPop$sex == "M"), sep = " "))
  print(paste0 ("Prod Females:", sum(ProdPop$sex =="F"), sep = " "))
  
}
 

if (g >= -25){
  
 MultPop_Fem <- ProdPop %>% filter(sex == "F" & herd == "A" | herd == "B", age >= AgeFirstMate) %>% top_frac(0.5, merit) ### need to be creaming off the percentage and not returning!!! ###
 MultPop_Males <- SPFpop %>% filter(sex == "M" & herd == "B", age >= AgeFirstMate) %>% top_frac(0.1, merit) ### again filtering the boars, but they still exist in the SPF population  ##
 MultPop <- rbind.fill(MultPop_Fem, MultPop_Males) #puts new Multpop with new piglets from SPF & Prod 
  
  Mult_Males <- MultPop %>% filter(sex == "M", age >= AgeFirstMate) ### are these required ??? #### maybe to keep the pigs that have been put into these herds.. ###
  Mult_Females <- MultPop %>% filter(sex == "F", age >= AgeFirstMate & age %% FarrowInt == rem) ## probably best to filter for age here, not realistic to redo it for each cycle ##
  
  NewMultPop <- CreatePiglets(MultPop_Males, Mult_Females, indexSD, g, paste0("Mult",g,"_"), littersize)
  
  MultPop <- rbind.fill(NewMultPop, MultPop) 
  
  MultPop <- ageing(MultPop)
  
  print(paste0("Mult Males:", sum(MultPop$sex == "M"), sep = " "))
  print(paste0 ("Mult Females:", sum(MultPop$sex =="F"), sep = " "))
  
} ###############################
  
  if (g >= -15){
    
    BW_Fem <- MultPop %>% filter(sex == "F", herd =="A" | herd == "B" | herd =="T", age >= AgeFirstMate) %>% top_frac(0.7, merit) ### need to be creaming off the percentage and not returning!!! ### Take not top merit but middling aniamls
    BW_Males <- SPFpop %>% filter(sex == "M" & herd == "T", age >= AgeFirstMate) %>% top_frac(0.1, merit) ### again filtering the boars, but they still exist in the SPF population  #
    BWpop <- rbind.fill(BWpop, BW_Fem, BW_Males) 
    
    BW_Fem <- BWpop %>% filter(sex == "F", age >= AgeFirstMate & age %% FarrowInt == rem) ## probably best to filter for age here, not realistic to redo it for each cycle ##
    BW_Males <- BWpop %>% filter(sex == "M", age >= AgeFirstMate) 
    
    NewBWpop <- CreatePiglets(BW_Males, BW_Fem, indexSD, g, paste0("BW",g,"_"), littersize)
    
    BWPop <- rbind(NewBWpop, BWpop)
    
    BWpop <- ageing(BWpop)
    
    print(paste0("BW Males:", sum(BWpop$sex == "M"), sep = " "))
    print(paste0 ("BW Females:", sum(BWpop$sex =="F"), sep = " "))
  
    
  }
    Allpop <- rbind (SPFpop, ProdPop, MultPop, BWpop)
  }
  return(Allpop) 
  }

#}  ## BurnIn close brace 


  
  
  

 

  
 
