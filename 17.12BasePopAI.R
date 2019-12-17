####### functions and running code ########

library (plyr)
library (dplyr)
library (purrr)
library (tidyverse)
library (genetics)


################## PARAMETERS/VALUES #########################
######## setting parameters/values required for code ########## 
AgeDist <- c(0.092,0.055,0.07,0.09,0.115,0.147,0.189,0.242) ###
indexSD <- 10 
littersize <- 12 ### consider function for changing the littersize and having a fertilisation proportion ###
AgeFirstMate <- 8
FarrowInt <- 5
rem <- AgeFirstMate - FarrowInt #### females to breed. rem = 3 which is selects them for cycles to breed in %% ######

genotypesA_susceptible <- c("AA") ## introduce after gene editing
allelesA_susceptible <- c("A")

genotypesB_susceptible <- c("BB") ## introduce after gene editing
allelesB_susceptible <- c("B")

#### ForwardSim Genotyping needs alleles for resistance ###

genotypesA <- c("AA", "Aa", "aA", "aa") ## introduce after gene editing
allelesA <- c("A", "a")

genotypesB <- c("BB", "Bb", "bB", "bb") ## introduce after gene editing
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
  #NucleusA$genoA <- sample(expectedGenotypes(as.genotype(genotypesA, sep = ""), alleles=allelesA, ploidy=2, sort=TRUE,haplotype=FALSE), n, replace = TRUE)
  NucleusA$genoA <- "Rr"
  NucleusA$genoB <- "Rr"
  #NucleusA$genoB <- sample(expectedGenotypes(as.genotype(genotypesB, sep = ""), alleles=allelesB, ploidy=2, sort=TRUE,haplotype=FALSE), n, replace = TRUE)
  NucleusA$age <- sample(rep(1:length(AgeDist) + 8), size = n, replace = TRUE, prob = AgeDist) # +8 so that piglets are ready for breeding more quickly ##3
  NucleusA$fate <- 1
  
  return (NucleusA)
}


CreatePiglets <- function (sires, dams, indexSD, genNo, label, littersize){
  nP <- length(dams$ID)*littersize ###
  piglets <- setNames(data.frame(matrix(nrow = nP, ncol =11)), c("ID", "gen", "herd", "sex", "merit", "genoA", "genoB", "age", "fate", "sire", "dam"))
  
  siresOrdered <- sires %>% filter(age >= AgeFirstMate) %>% arrange(desc(merit)) #### need to assign sires to the relevant dam 
  siresOrdered2 <- sample(1:length(siresOrdered$ID), length(dams$ID), replace = TRUE)
  siresOrdered2 <- sires[c(siresOrdered2), ]  ############ don't understand this code well #########
  
  ### creating substring allows selecting of each character ###
  
  # split the sire and dam genotypes into the separate alleles - A 
  ### will always select the allele on the left. not a true shuffling or random meiotic pattern? how to include recombination factors???? 
  sgA <- siresOrdered2$genoA
  s1A <- substr(sgA, 1, 1)
  s2A <- substr(sgA, 2, 2)
  dgA <- dams$genoA
  d1A <- substr(dgA, 1, 1)
  d2A <- substr(dgA, 2, 2)
  # Sample which allele of the PRRS resistant genotype is inherited from each parent
  sireAlleleA <- sample(c(1,2),nP,replace = TRUE)
  damAlleleA <- sample(c(1,2),nP,replace = TRUE)
  
  # split the sire and dam genotypes into the separate alleles - B
  sgB <- siresOrdered2$genoB
  s1B <- substr(sgB, 1, 1)
  s2B <- substr(sgB, 2, 2)
  dgB <- dams$genoB
  d1B <- substr(dgB, 1, 1)
  d2B <- substr(dgB, 2, 2)
  # Sample which allele of the PRRS resistant genotype is inherited from each parent
  sireAlleleB <- sample(c(1,2),nP,replace = TRUE)
  damAlleleB <- sample(c(1,2),nP,replace = TRUE)
  
  piglets$ID <- paste0(label, seq(1:nP)) ## label refers to generation of creation as identifier, doesn't need to be additive ###
  piglets$gen <- genNo 
  piglets$herd <- siresOrdered2$herd ## showing male herd is enough to know lineage? ###
  piglets$sex <- sample (c("M", "F"), nP, replace = TRUE, prob = c(0.5, 0.5))
  piglets$merit <- 0.5 * (siresOrdered2$merit + dams$merit) + rnorm(nP, 0, indexSD) ##### is this variation a realistic meiotic figure? #### are these the actual dam and sire merits? ### how to check...?##
  piglets$genoA <- paste0(ifelse(sireAlleleA==1,s1A,s2A),ifelse(damAlleleA==1,d1A,d2A))
  piglets$genoB <- paste0(ifelse(sireAlleleB==1,s1B,s2B),ifelse(damAlleleB==1,d1B,d2B)) ###paste0 puts in letter not the logical result
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

SPFpop <- BasePop(9000, indexSD, AgeDist)  #### creates base population of piggys ### start with ~10,000 to have enough for flowing down ### 

  ProdPop <- data.frame(matrix(ncol = 11, nrow = 0))
  colnames(ProdPop) <- c("ID", "gen", "herd", "sex", "merit", "genoA", "genoB", "age", "fate", "sire", "dam")
  
  MultPop <- data.frame(matrix(ncol = 11, nrow = 0))
  colnames(MultPop) <- c("ID", "gen", "herd", "sex", "merit", "genoA", "genoB", "age", "fate", "sire", "dam")
  
  BWpop <- data.frame(matrix(ncol = 11, nrow = 0))
  colnames(BWpop) <- c("ID", "gen", "herd", "sex", "merit", "genoA", "genoB", "age", "fate", "sire", "dam")
  
  CommercialPop <- data.frame(matrix(ncol = 11, nrow = 0))
  colnames(BWpop) <- c("ID", "gen", "herd", "sex", "merit", "genoA", "genoB", "age", "fate", "sire", "dam")
  
  for (g in -60:0){
  print(paste0("Gen:", g, sep = " "))
    
  ##### separate into A herd ########### 
      
   ### select males from top percentage. Should leave other breeding age pigs for below tiers
  BreedSPFA_Male <- SPFpop %>% filter (sex == "M" & herd == "A", age >= AgeFirstMate) %>% top_frac(0.02, merit) 
  SPFA_Fem <- filter(SPFpop, sex == "F" & herd == "A", age >= AgeFirstMate)  %>% top_frac(0.2, merit) ### remove these pigs from prod tier. Will leave some pigs that can be bred from!
  BreedSPFA_Fem <- filter(SPFA_Fem, age %% FarrowInt == rem) #already filtered on age and merit
  
  newGenSPFNucA <- CreatePiglets(BreedSPFA_Male, BreedSPFA_Fem, indexSD, g, paste0("SPF",g,"_"), littersize)

  
  ##### separate into B herd ###########
  BreedSPFB_Male <- SPFpop %>% filter (sex == "M" & herd == "B", age >= AgeFirstMate) %>% top_frac(0.02, merit) 
  SPFB_Fem <- filter(SPFpop, sex == "F" & herd == "B" , age >= AgeFirstMate) %>% top_frac(0.2, merit) ### remove these pigs from prod tier. Will leave some pigs that can be bred from!
  BreedSPFB_Fem <- filter(SPFB_Fem, age %% FarrowInt == rem) #already filtered on age and merit
  
  newGenSPFNucB <- CreatePiglets(BreedSPFB_Male, BreedSPFB_Fem, indexSD, g, paste0("SPF",g,"_"), littersize)
  
  
  ##### separate into T herd ###########
  BreedSPFT_Male <- SPFpop %>% filter (sex == "M" & herd == "T", age >= AgeFirstMate) %>% top_frac(0.02, merit) 
  SPFT_Fem <- filter(SPFpop, sex == "F" & herd == "T", age >= AgeFirstMate ) %>% top_frac(0.2, merit) ### remove these pigs from prod tier. Will leave some pigs that can be bred from!
  BreedSPFT_Fem <- filter(SPFT_Fem, age %% FarrowInt == rem) #already filtered on age and merit
  
  newGenSPFNucT <- CreatePiglets(BreedSPFT_Male, BreedSPFT_Fem, indexSD, g, paste0("SPF",g,"_"), littersize)
  ### new SPF pop piglets join rows at the end of the code
 
                  #########
  SPF_Breeders <- rbind.fill(BreedSPFA_Fem, BreedSPFB_Fem, BreedSPFT_Fem) ### animals for SPF breeding
  ### will rejoin SPFpop at end of code incase they are no longer in top % for mating and can be used in below tiers  
  ### males are not removed as they can inseminate multiple tiers through an AI program
  
  
  print(paste0("SPF Males:", sum(SPFpop$sex == "M"), sep = " "))
  print(paste0 ("SPF Females:", sum(SPFpop$sex =="F"), sep = " "))
  
  
if (g >= -40){ #### selecting all the breeding sows above removes the breeders for Prod breeding. Need to select elsewhere ###
  
  SPFpop <-  anti_join(SPFpop, SPF_Breeders, by = "ID") ## removes breeding SPF breeding animals from SPF pop,
     # selection will exclude nucleus breeding animals ## keeps out for mult & BW tier 
  
  SPF_ProdPop_Fem <- SPFpop %>% filter(sex == "F", herd == "A", age >= AgeFirstMate) %>% top_frac(0.25, merit) #takes away top 25% for PN 
  ProdPop_Fem_Breed <- filter(SPF_ProdPop_Fem, age %% FarrowInt == rem) #already filtered on age and merit
  
  SPF_ProdPop_Males <- SPFpop %>% filter(sex == "M", herd == "A", age >= AgeFirstMate) %>% top_frac(0.1, merit) #takes top 10% available boars, SPF pop excluded.
  #SPFpop <-  anti_join(SPFpop, SPF_ProdPop_Males, by = "ID") ## removes prod pop males from the SPFpop #can't be reused in SPFpop
  
  SPFpop <-  anti_join(SPFpop, ProdPop_Fem_Breed, by = "ID") ## returns all the pigs in SPFpop that are not present in SPF_ProdNuc_Breeders_M. ##SPF pop loses the animals transferred
  ##put in print function to see how many females are being bred in each generation ##
  
  ProdPop <- rbind.fill(ProdPop, SPF_ProdPop_Fem, SPF_ProdPop_Males) ### NO animals ######
  
  Prod_Males <- ProdPop %>% filter(sex == "M")
  Prod_Females <- ProdPop %>% filter(sex == "F")
  
  NewProdPop <- CreatePiglets(Prod_Males, Prod_Females, indexSD, g, paste0("Prod",g,"_"), littersize)
  
  AddProd_Females <- NewProdPop %>% filter(sex == "F") %>% top_frac(0.3, merit) ## takes the top 20% of new females created, all that is required. #### increase if I want more PN breeding ###
  AddProd_Males <- NewProdPop %>% filter(sex == "M") %>% top_frac(0.1, merit) ## top 10%% of new males are added to the production population
  
  ProdPop <- rbind.fill(ProdPop, AddProd_Females, AddProd_Males) ### may be confusing with adding the above rows to data frame in wrong order? ###
  
  ## need culling function to remove low perfroming animals from breeding program. Set a population max figure to keep consistent. 
  ## Will reduce numbers but avoids low merit sows from early tiers being used until death?
  
  print(paste0("Prod Males:", sum(ProdPop$sex == "M"), sep = " "))
  print(paste0 ("Prod Females:", sum(ProdPop$sex =="F"), sep = " "))
}

  
if (g >= -25){
  
  ### A or B females as can be bred within mult tier or from prod?? ## no transfer of Mult pigs back into breeding data frame so only A females used, piglets will be B herd. 
 
  ProdPop <-  anti_join(ProdPop, ProdPop_Fem_Breed, by = "ID") ### removes breeding ProdPop animals before selecting for MultPop. Are they 
  
  Prod_MultPop_Females <- ProdPop %>% filter(sex == "F", herd == "A" | herd == "B", age >= AgeFirstMate) %>% top_frac(0.5, merit) #selects top 50% for mult tier from Prod. 
  MultPop_Fem_Breed <- filter(Prod_MultPop_Females, age %% FarrowInt == rem) #already filtered on age and merit
  
  SPF_MultPop_Males <- SPFpop %>% filter(sex == "M", herd == "B", age >= AgeFirstMate) %>% top_frac(0.1, merit) #takes top 10% available boars, SPF pop excluded.## breeders removed above
  #SPFpop <-  anti_join(SPFpop, SPF_MultPop_Males, by = "ID") ## removes mult pigs from SPF pop 
  
  MultPop <- rbind.fill(MultPop, MultPop_Fem_Breed, SPF_MultPop_Males) #puts new Multpop with new piglets from SPF & Prod. Are these females now permanently in this tier?
  
  Mult_Males <- MultPop %>% filter(sex == "M") ### are these required ??? #### maybe to keep the pigs that have been put into these herds.. ###
  Mult_Females <- MultPop %>% filter(sex == "F" & age %% FarrowInt == rem) ## probably best to filter for age here, not realistic to redo it for each cycle ##
  
  NewMultPop <- CreatePiglets(Mult_Males, Mult_Females, indexSD, g, paste0("Mult",g,"_"), littersize)
  
  MultPop <- rbind.fill(NewMultPop, MultPop)
  MultPop <- MultPop %>% filter(sex == "F")
  

  print(paste0("Mult Males:", sum(MultPop$sex == "M"), sep = " "))
  print(paste0 ("Mult Females:", sum(MultPop$sex =="F"), sep = " "))
  
} ###############################
  
  if (g >= -10){
    
    MultPop <-  anti_join(MultPop, MultPop_Fem_Breed, by = "ID") #### need to rejoin breeding animals?
    
    BW_Fem <- MultPop %>% filter(sex == "F", herd == "B" | herd =="A", age >= AgeFirstMate) %>% top_frac(0.7, merit) ### need to be creaming off the percentage and not returning!!! 
    ## no A herd as they must have been bred in Multt tier ### no T herd either as the piglets are transferred directly to commercial tier ###
   
    #### not selecting correctly here.... 
    Mult_BW_Fem_Breed <- BW_Fem %>% filter(age %% FarrowInt == rem) ## probably best to filter for age here, not realistic to redo it for each cycle ## No BW_Fem_Breed as filtered here
    
    SPF_BWpop_Males <- SPFpop %>% filter(sex == "M", herd == "T", age >= AgeFirstMate) %>% top_frac(0.1, merit) #take away top 2% for SPF breeding, select next 10% for PN 
    #SPFpop <-  anti_join(SPFpop, SPF_MultPop_Males, by = "ID") ## returns all the pigs in SPFpop that are not present in SPF_ProdNuc_Breeders_M
    
    BWpop <- rbind.fill(BWpop, Mult_BW_Fem_Breed) 
    
    BW_Fem_Breed <- BWpop %>% filter(sex == "F", age %% FarrowInt == rem)

    MultPop <- anti_join(MultPop, Mult_BW_Fem_Breed, by = "ID") ### remove BW bred females, leave rest in multPop
    
    MultPop <- rbind.fill(MultPop, MultPop_Fem_Breed) ### rejoin bred multpop females for next generation ###
    
    CommercialPop <- CreatePiglets(SPF_BWpop_Males, BW_Fem_Breed, indexSD, g, paste0("BW",g,"_"), littersize) #piglets will be T herd, females from B herd. 
    ### maybe transfer BW males and retain some BW females for breeding if they are of better merit than piglets born??? ###
    
    ## maybe need culling function to remove low perfroming animals. May also not need population max as these animals are passed to the commercial breeders for growing. 
    
    print(paste0("BW Males:", sum(BWpop$sex == "M"), sep = " "))
    print(paste0 ("BW Females:", sum(BWpop$sex =="F"), sep = " "))
    
    print(paste0 ("Commercial Pigs:", nrow(CommercialPop), sep = " "))
    
  }
  
  
  SPFpop <- rbind.fill(SPFpop, SPF_Breeders, newGenSPFNucT, newGenSPFNucA, newGenSPFNucB) ### includes all SPF pigs, no selection. Cull after ageing function
  ## All SPF pigs. SPF pigs moved to another tier are excluded by antijoin from SPF tables
  ### remerges all SPFpops, older pigs with lower merit will not get stuck in breeding population - all SPF pigs as a single herd!!! ###
  ## need to rebind as top 2% of boars from total SPF pop will not remain the same!
  ## put breeders back with SPFpop, if they are not used for SPF pop breeding next round they may be used for PN. Need to select breeders from same pool as before ##
  
  
  SPFpop <- data.frame(ageing(SPFpop))
  
  SPFpopA_females <- SPFpop %>% filter(sex == "F", herd == "A") %>% filter(age >= 8) %>% top_n(4500, merit) ### filter breeding animals o ensure they remain ### this code works!!!
  SPFpopB_females <- SPFpop %>% filter(sex == "F", herd == "B") %>% filter(age >= 8) %>% top_n(2500, merit) ### filter breeding animals o ensure they remain ### this code works!!!
  SPFpopT_females <- SPFpop %>% filter(sex == "F", herd == "T") %>% filter(age >= 8) %>% top_n(2500, merit) ### filter breeding animals o ensure they remain ### this code works!!!
  
  SPFpopA_males <- SPFpop %>% filter(sex == "M", herd == "A") %>% filter(age >= 8) %>% top_n(500, merit) #### may need to ensure not all old pigs are culled so that there are enough older animals for breeding!
  SPFpopB_males <- SPFpop %>% filter(sex == "M", herd == "B") %>% filter(age >= 8) %>% top_n(500, merit) #### may need to ensure not all old pigs are culled so that there are enough older animals for breeding!
  SPFpopT_males <- SPFpop %>% filter(sex == "M", herd == "T") %>% filter(age >= 8) %>% top_n(500, merit) #### may need to ensure not all old pigs are culled so that there are enough older animals for breeding!
  
  SPFpopA_piglets <- SPFpop %>% filter(herd == "A") %>% filter(age < 8) %>% top_n(1000, merit) #### culling piglets before maturity..
  SPFpopB_piglets <- SPFpop %>% filter(herd == "B") %>% filter(age < 8) %>% top_n(1000, merit)
  SPFpopT_piglets <- SPFpop %>% filter(herd == "T") %>% filter(age < 8) %>% top_n(1000, merit)
  
  ## SPF pop retainers need to consider herd distributions with such low numbers initially breeding
  
  topSPFpop <-  rbind.fill(SPFpopA_females, SPFpopB_females, SPFpopT_females, SPFpopA_males, SPFpopB_males, SPFpopT_males, SPFpopA_piglets, SPFpopB_piglets, SPFpopT_piglets) ## remove low performing males from SPF population 
  ### careful that older pigs are not killed as piglets will be higher merit but can't breed/ put in selection for if over 16 months???
  #### may need to ensure not all old pigs are culled so that there are enough older animals for breeding!
  #already only 1000 males selected. 
  
  
  SPFpop <- semi_join(topSPFpop, SPFpop, by = "ID") 
  
  ### makes SPFpop only contain the topSPFpop pigs
  ########
  ### have anti joins for SPF vs Prod, Prod Vs mult and mult Vs BW to ensure no cross over in pig populations???? #####3
  
  SPFpopA_females <- NULL ### remove tables to prevent accumulation of data slowing code
  SPFpopB_females <- NULL
  SPFpopT_females <- NULL
  SPFpopA_males <- NULL
  SPFpopB_males <- NULL
  SPFpopT_males <- NULL
  SPFpopA_piglets <- NULL
  SPFpopB_piglets <- NULL
  SPFpopT_piglets <- NULL
  topSPFpop <- NULL ### remove tables to prevent accumulation of data slowing code
  SPF_Breeders <- NULL
  
    ########################
  
  ProdPop <- data.frame(ageing(ProdPop)) 
  
  ProdPopBreeders <- ProdPop %>% filter(sex == "F") %>% filter(age >= AgeFirstMate) %>% top_n(33000, merit) ### filter breeding animals o ensure they remain ### this code works!!!
  ProdPopPiglets <- ProdPop %>% filter(sex == "F") %>% filter(age < AgeFirstMate) %>% top_n(15000, merit)
  ProdPopMales <- ProdPop %>% filter (sex == "M") %>% filter(age > AgeFirstMate) %>% top_n(800, merit)
  ProdPop <- rbind.fill(ProdPopBreeders, ProdPopPiglets) ### males for prod pop are taken each time from SPF_A. Should I set this to retain some males??
  
  #ProdPop <- anti_join (ProdPop, SPFpop, by = "ID") ## no animals should be able to exist in both tiers. Checks females, males are not transferred as it is an AI program. 
  #############
  
  MultPop <- ageing(MultPop)
  MultPopBreeders <- MultPop %>% filter(sex == "F") %>% filter(age >= AgeFirstMate) %>% top_n(80000, merit)
  MultPop <- semi_join(MultPopBreeders, MultPop)
  ### may need to include population limitation as for SPF and ProdPops
 # MultPop <- anti_join(MultPop, ProdPop, by = "ID")
  ### remerge multpop breeders so they can be used again??
  
  
  BWpop <- ageing(BWpop)
  ### remerge BWpop breeders so they can be used again??
  
  
  Allpop <- rbind.fill (SPFpop, ProdPop, MultPop, BWpop)
   
  # CommercialPop <- BWpop
  }

  ### have anti joins for SPF vs Prod, Prod Vs mult and mult Vs BW to ensure no cross over in pig populations???? #####3
  ## must still remove ProdNuc, Multpop and BWpop males from the SPFpop. Once they transtition they can't return 
  





  

 

  
 
