AgeDist <- c(0.092,0.055,0.07,0.09,0.115,0.147,0.189,0.242) ###
indexSD <- 10 
littersize <- 12 ### consider function for changing the littersize and having a fertilisation proportion ###
AgeFirstMate <- 8
FarrowInt <- 5
rem <- AgeFirstMate - FarrowInt #### females to breed. rem = 3 which is selects them for cycles to breed in %% ######

###############################
####################

BasePop <- function(n, indexSD, AgeDist) {
  NucleusA <- setNames(data.frame(matrix(nrow = n, ncol =9)), c("ID", "gen", "herd", "sex", "merit", "genoA", "genoB", "age", "fate"))
  
  NucleusA$ID <- rep(1:n)
  NucleusA$gen <- 0
  NucleusA$herd <- sample(c("A", "B", "T"), size = n, replace = TRUE, prob = c(0.5, 0.25, 0.25))
  NucleusA$sex <- sample (c("M", "F"), n, replace = TRUE, prob = c(0.3, 0.7))
  NucleusA$merit <- rnorm(n, mean = 0, sd = indexSD)
  NucleusA$genoA <- "Rr"
  NucleusA$genoB <- "Rr"
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