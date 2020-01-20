# FORWARD SIMULATION WITH EDITING

## Create BasePop - generations -70 to 0



#BurnInPop <- read.csv("BurnInPop2020-01-14.csv") ### source ("") ### setwd ("")

SPFpop <- BurnIn$SPFpop
Prodpop <- BurnIn$ProdPop
MultPop <- BurnIn$MultPop
BWpop <- BurnIn$BWpop 

### reset all populations to BurnIn read file

#### changing litter size is effectively increasing the number of donor females, flushing pigs. alter for different techniques. Think about!!!! 
### Doesn't matter if the recip pig isn't recorded, flushed female can have more piglets???  

for (g in 0:120){
  
  print(paste0("Gen:", g, sep = " "))
  
  ### bias selection towards resistant or hets in editing population?
  ### don't breed with A/A animals??? 
  ### add in selection for having a resistance allele? resistant > hets >susceptible
  ### select males from top percentage. Should leave other breeding age pigs for below tiers
  ### if for top_frac whereby if a/a not on merit, otherwise only top_frac of merit column
  
  if (nrow(filter(SPFpop, sex == "F" & herd == "A") %>% filter(age >= AgeFirstMate & age %% FarrowInt == rem)  %>% top_frac(0.5, merit)) >2) {
    
    BreedSPFA_Male <- SPFpop %>% filter (sex == "M" & herd == "A") %>% filter(age >= AgeFirstMate) # %>% top_frac(0.05, merit)
    
    if (nrow(filter(BreedSPFA_Male, genoA == "a/a")) >= 5) {
      BreedSPFA_Male <- BreedSPFA_Male %>% filter(genoA == "a/a") # %>% top_frac(0.5, merit)  
    } else if (nrow(filter(BreedSPFA_Male, genoA != "A/A")) > 10) {
      BreedSPFA_Male <- BreedSPFA_Male %>% filter(genoA != "A/A")  %>% top_frac(0.2, merit) 
    } else  (BreedSPFA_Male <- BreedSPFA_Male %>% top_frac(0.05, merit) )  
    
    BreedSPFA_Fem <- filter(SPFpop, sex == "F" & herd == "A") %>% filter(age >= AgeFirstMate & age %% FarrowInt == rem)  %>% top_frac(0.5, merit) ### remove these pigs from prod tier. Will leave some pigs that can be bred from!
    
    newGenSPFNucA <- CreatePiglets(BreedSPFA_Male, BreedSPFA_Fem, indexSD, g, paste0("SPFA",g,"_"), littersize)
    
    newGenSPFNucA <- edit_genes(newGenSPFNucA, 0.6, 0.6) ### editing of piglets - done to zygotes but is easier to code post creation
    
    SPFNucA_Fem <- filter(newGenSPFNucA, sex == "F")  %>% top_frac(0.25, merit)
    SPFNucA_Fem1 <- newGenSPFNucA %>% filter(genoA == "a/a" &  sex == "F") # %>% top_frac(1, merit)
    SPFNucA_Fem <- rbind.fill(SPFNucA_Fem, SPFNucA_Fem1) %>% distinct()
    
    SPFNucA_Male <- filter(newGenSPFNucA, sex == "M") %>% top_frac (0.05, merit)
    SPFNucA_Male1 <-  newGenSPFNucA %>% filter(genoA == "a/a" &  sex == "M") #%>% top_frac(1, merit)
    SPFNucA_Male <- rbind.fill(SPFNucA_Male, SPFNucA_Male1) %>% distinct()
    ### these pigs are going into the SPF population. Then the next iteration selects piglets on merit. 
    
  } ### only runs breeding for A if there are 2 females present. Keeps stopping with no sows!!! ####
  
  
  ################### 
  
  ##### separate into B herd ###########
  ### could avoid editing in B herd as there will be flow down from prodpop to multpop and hets can then breed with T hets...
  
  if (nrow(SPFpop %>% filter(sex == "F" & herd == "B") %>% filter(age >= AgeFirstMate & age %% FarrowInt == rem)) >2) {
    
    BreedSPFB_Male <- SPFpop %>% filter (sex == "M" & herd == "B") %>% filter(age >= AgeFirstMate) # %>% top_frac(0.05, merit) 
    
    if (nrow(filter(BreedSPFB_Male, genoA == "a/a")) >= 5) {
      BreedSPFB_Male <- BreedSPFB_Male %>% filter(genoA == "a/a") # %>% top_frac(0.1, merit) 
    } else if (nrow(filter(BreedSPFB_Male, genoA != "A/A")) > 10) {
      BreedSPFB_Male <- BreedSPFB_Male %>% filter(genoA != "A/A")  %>% top_frac(0.2, merit) 
    } else  (BreedSPFB_Male <- BreedSPFB_Male %>% top_frac(0.05, merit))  
    
    BreedSPFB_Fem <- filter(SPFpop, sex == "F" & herd == "B") %>% filter(age >= AgeFirstMate & age %% FarrowInt == rem) # %>% top_frac(0.25, merit) ### remove these pigs from prod tier. Will leave some pigs that can be bred from!
    
    if (nrow(filter(BreedSPFB_Fem, genoA == "a/a")) >= 3) {
      BreedSPFB_Fem <- BreedSPFB_Fem %>% filter(genoA == "a/a")  #%>% top_frac(0.5, merit) 
    } else if (nrow(filter(BreedSPFB_Fem, genoA != "A/A")) > 4) {
      BreedSPFB_Fem <- BreedSPFB_Fem %>% filter(genoA != "A/A") %>% top_frac(0.5, merit) 
    } else  (BreedSPFB_Fem <- BreedSPFB_Fem %>% top_frac(0.25, merit))
    
    newGenSPFNucB <- CreatePiglets(BreedSPFB_Male, BreedSPFB_Fem, indexSD, g, paste0("SPFB",g,"_"), littersize)
    
    newGenSPFNucB <- edit_genes(newGenSPFNucB, 0.6, 0.6)
    
    SPFNucB_Fem <- filter(newGenSPFNucB, sex == "F")  %>% top_frac(0.2, merit)
    SPFNucB_Fem1 <- newGenSPFNucB %>% filter(genoA == "a/a" &  sex == "F") # %>% top_frac(1, merit)
    SPFNucB_Fem <- rbind.fill(SPFNucB_Fem, SPFNucB_Fem1) %>% distinct()
    
    SPFNucB_Male <- filter(newGenSPFNucB, sex == "M") %>% top_frac (0.05, merit)
    SPFNucB_Male1 <-  newGenSPFNucB %>% filter(genoA == "a/a" &  sex == "M") #%>% top_frac(1, merit)
    SPFNucB_Male <- rbind.fill(SPFNucB_Male, SPFNucB_Male1) %>% distinct()
  }
  
  ########################################
  
  ##### separate into T herd ###########
  
  if (nrow(SPFpop %>% filter(sex == "F" & herd == "T") %>% filter(age >= AgeFirstMate & age %% FarrowInt == rem)) >2) {
    
    BreedSPFT_Male <- SPFpop %>% filter (sex == "M" & herd == "T") %>% filter(age >= AgeFirstMate) # %>% top_frac(0.05, merit) 
    
    if (nrow(filter(BreedSPFT_Male, genoA == "a/a")) >= 3) {
      BreedSPFT_Male <- BreedSPFT_Male %>% filter(genoA == "a/a")  %>% top_frac(0.5, merit) 
    } else if (nrow(filter(BreedSPFT_Male, genoA != "A/A")) > 4) {
      BreedSPFT_Male <- BreedSPFT_Male %>% filter(genoA != "A/A")  %>% top_frac(0.5, merit) 
    } else  (BreedSPFT_Male <- BreedSPFT_Male %>% top_frac(0.05, merit) )  
    
    BreedSPFT_Fem <- filter(SPFpop, sex == "F" & herd == "T") %>% filter(age >= AgeFirstMate & age %% FarrowInt == rem) # %>% top_frac(0.25, merit) ### remove these pigs from prod tier. Will leave some pigs that can be bred from!
    
    if (nrow(filter(BreedSPFT_Fem, genoA == "a/a")) >= 3) {
      BreedSPFT_Fem <- BreedSPFT_Fem %>% filter(genoA == "a/a") # %>% top_frac(0.5, merit) 
    } else if (nrow(filter(BreedSPFT_Fem, genoA != "A/A")) > 4) {
      BreedSPFT_Fem <- BreedSPFT_Fem %>% filter(genoA != "A/A")  %>% top_frac(0.2, merit) 
    } else  (BreedSPFT_Fem <- BreedSPFT_Fem %>% top_frac(0.25, merit))
    newGenSPFNucT <- CreatePiglets(BreedSPFT_Male, BreedSPFT_Fem, indexSD, g, paste0("SPFT",g,"_"), littersize)
    
    newGenSPFNucT <- edit_genes(newGenSPFNucT, 0.6, 0.6)
    
    
    SPFNucT_Fem <- filter(newGenSPFNucT, sex == "F")  %>% top_frac(0.2, merit)
    SPFNucT_Fem1 <- newGenSPFNucT %>% filter(genoA == "a/a" &  sex == "F") # %>% top_frac(1, merit)
    SPFNucT_Fem <- rbind.fill(SPFNucT_Fem, SPFNucT_Fem1) %>% distinct()
    
    SPFNucT_Male <- filter(newGenSPFNucT, sex == "M") %>% top_frac (0.1, merit)
    SPFNucT_Male1 <-  newGenSPFNucT %>% filter(genoA == "a/a" &  sex == "M") #%>% top_frac(1, merit)
    SPFNucT_Male <- rbind.fill(SPFNucT_Male, SPFNucT_Male1) %>% distinct()
    
  }
  ################################
  
  # if (g == 0) { 
  #   SPFpopAll <- data.frame(SPFpop)
  # } else  (SPFpopAll <- rbind.fill(SPFpopAll, SPFpop))
  
  SPF_Breeders <- rbind.fill(BreedSPFA_Fem, BreedSPFB_Fem, BreedSPFT_Fem) 
  
  ############# need to bias to preferential mating of resistance alleles in ProdPop 
  
  if (g >= 0){ 
    #### selecting all the breeding sows above removes the breeders for Prod breeding. Need to select elsewhere ###
    #### check to ensure that resistant pigs aren't all being taken by SPF pop. Bias Males towards resistant for AI in Prodpop
    # & g %% 2 == 0 ### attempt to not remove so many from SPFpop as ProdPop has enough males  
    
    SPFpop <-  anti_join(SPFpop, SPF_Breeders, by = "ID") ## removes breeding SPF breeding animals from SPF pop,
    
    if (nrow(SPFpop %>% filter(sex == "F" & herd == "A") %>% filter(age >= AgeFirstMate & age %% FarrowInt == rem) %>% top_frac(0.5, merit)) >2) {
      
      SPF_ProdPop_Fem <- SPFpop %>% filter(sex == "F" & herd == "A") %>% filter(age >= AgeFirstMate & age %% FarrowInt == rem) #%>% top_frac(0.25, merit)
      ### removing selection affects flow massively!! speeds code also
      
      if (nrow(filter(SPF_ProdPop_Fem, genoA == "a/a")) >= 20) {
        SPF_ProdPop_Fem <- SPF_ProdPop_Fem %>% filter(genoA == "a/a")  %>% top_frac(0.5, merit)
      } else if (nrow(filter(SPF_ProdPop_Fem, genoA != "A/A")) > 40) {
        SPF_ProdPop_Fem <- SPF_ProdPop_Fem %>% filter(genoA != "A/A")  %>% top_frac(0.25, merit) 
      } else  (SPF_ProdPop_Fem <- SPF_ProdPop_Fem %>% top_frac(0.25, merit))  
      
      
      
      SPFpop <-  anti_join(SPFpop, SPF_ProdPop_Fem, by = "ID") ## SPF pop loses the animals transferred to prod pop. Not reintegrated and cannot be used in SPF in future 
      
      SPF_ProdPop_Males <- SPFpop %>% filter(sex == "M" & herd == "A") %>% filter(age >= AgeFirstMate) # %>% top_frac(0.1, merit) 
      #takes top 10% available boars ### genotype?
      
      if (nrow(filter(SPF_ProdPop_Males %>% filter(age >= AgeFirstMate & genoA == "a/a"))) >= 10) {
        SPF_ProdPop_Males <- SPF_ProdPop_Males %>% filter(genoA == "a/a")  # %>% top_frac(0.5, merit) 
      } else if (nrow(filter(SPF_ProdPop_Males, genoA != "A/A")) > 20) {
        SPF_ProdPop_Males <- SPF_ProdPop_Males %>% filter(genoA != "A/A")  %>% top_frac(0.5, merit) 
      } else  (SPF_ProdPop_Males <- SPF_ProdPop_Males %>% top_frac(0.05, merit) )  
      
      #SPFpop <-  anti_join(SPFpop, SPF_ProdPop_Males, by = "ID") ## removes prod pop males from the SPFpop #can't be reused in SPFpop 
      ### don't need as AI is performed on SPF pops
      
      ProdPop <- rbind.fill(ProdPop, SPF_ProdPop_Fem) ### No males from SPF stored as used by AI only . Only PN bred males will be in PN ######
      
      Prod_Males <- rbind.fill(ProdPop, SPF_ProdPop_Males) %>% filter(sex == "M") %>% filter(age >= AgeFirstMate) #%>% top_frac(0.15, merit) # selects for males from ProdPop and SPF Nucleus ## never integrated with full ProdPop
      if (nrow(filter(Prod_Males %>% filter(age >= AgeFirstMate & genoA == "a/a"))) >= 10) {
        Prod_Males <- Prod_Males %>% filter(genoA == "a/a")  # %>% top_frac(0.5, merit) 
      } else if (nrow(filter(Prod_Males, genoA != "A/A")) > 20) {
        Prod_Males <- Prod_Males %>% filter(genoA != "A/A")  %>% top_frac(0.5, merit) 
      } else  (Prod_Males <- Prod_Males %>% top_frac(0.05, merit) )  
      
      Prod_Females <- ProdPop %>% filter(sex == "F") %>% filter(age >= AgeFirstMate & age %% FarrowInt == rem) %>% top_frac(0.5, merit) %>% mutate(ID = as.character(ID))### already over age first mate, ensure they are correct farrowing interval. 
      
      NewProdPop <- CreatePiglets(Prod_Males, Prod_Females, indexSD, g, paste0("Prod",g,"_"), littersize)
      
      NewProd_Females <- NewProdPop %>% filter(sex == "F") %>% top_frac(0.2, merit) ## takes the top 20% of new females created, all that is required. #### increase if I want more PN breeding ###
      NewProd_Females <- rbind.fill(NewProd_Females, (filter(NewProdPop, sex == "F" & genoA == "a/a")))
      
      NewProd_Males <- NewProdPop %>% filter(sex == "M") %>% top_frac(0.05, merit) ## top 10%% of new males are added to the production population
      NewProd_Males <- rbind.fill(NewProd_Males, (filter(NewProdPop, sex == "M" & genoA =="a/a")))
      
      ProdPop <- rbind.fill(ProdPop, NewProd_Females, NewProd_Males) ### combine all animals to stay in ProdPop ### prod_Females must be added in after MultPophas selected ##
      
      # if (g == 0) { 
      #    ProdPopAll <- data.frame(ProdPop)
      #  } else  (ProdPopAll <- rbind.fill(ProdPopAll, ProdPop))
      
    }
  }
  
  if (g >= 0){
    
    ### A or B females as can be bred within mult tier or from prod?? ## no transfer of Mult pigs back into breeding data frame so only A females used, piglets will be B herd. 
    ### Need to make sure Prod_Females are put back into ProdPop at the end of generation! 
    
    ## 'while' function may be needed here to no dwindle prodpop stocks
    
    ProdPop <- anti_join(ProdPop, Prod_Females, by = "ID") ## removes pigs bred in ProdPop from being bred again in same generation ## should only be top 50% of breeding animals
    
    Prod_MultPop_Females <- ProdPop %>% filter(sex == "F" & herd == "A" | herd == "B") %>% filter(age >= AgeFirstMate & age %% FarrowInt == rem) %>% top_frac(0.5, merit) #selects top 50% for mult tier from Prod. 
    ProdPop <- anti_join(ProdPop, Prod_MultPop_Females, by = "ID") ## removes pigs transferred down to multpop from ProdPop permanently
    
    SPF_MultPop_Males <- SPFpop %>% filter(sex == "M" & herd == "B") %>% filter(age >= AgeFirstMate) %>% top_frac(0.1, merit) #takes top 10% available boars, SPF pop excluded.## breeders removed above
    # if (nrow(filter(SPF_MultPop_Males %>% filter(age >= AgeFirstMate & genoA == "a/a"))) >= 50) {
    #   SPF_MultPop_Males <- SPF_MultPop_Males %>% filter(genoA == "a/a") # %>% top_frac(0.5, merit) 
    #  } else if (nrow(filter(SPF_MultPop_Males, genoA != "A/A")) > 100) {
    #   SPF_MultPop_Males <- SPF_MultPop_Males %>% filter(genoA != "A/A")  %>% top_frac(0.5, merit) 
    #  } else  (SPF_MultPop_Males <- SPF_MultPop_Males %>% top_frac(0.05, merit) )
    
    SPFpop <-  anti_join(SPFpop, SPF_MultPop_Males, by = "ID") ## removes mult pigs from SPF pop 
    
    MultPop <- rbind.fill(MultPop, Prod_MultPop_Females, SPF_MultPop_Males) #puts new Multpop with new piglets from SPF & Prod. Males are moved as not AI
    
    Mult_Males <- MultPop %>% filter(sex == "M") %>% filter(age >= AgeFirstMate) %>% top_frac(0.1, merit) ### are these required ??? #### maybe to keep the pigs that have been put into these herds.. ###
    Mult_Females <- MultPop %>% filter(sex == "F") %>% filter (age >= AgeFirstMate & age %% FarrowInt == rem) %>% top_frac(0.5, merit) ## probably best to filter for age here, not realistic to redo it for each cycle ##
    
    NewMultPop <- CreatePiglets(Mult_Females, Mult_Females, indexSD, g, paste0("Mult",g,"_"), littersize)
    
    NewMult_Females <- NewMultPop %>% filter(sex == "F") %>% top_frac(0.25, merit) ## takes the top 20% of new females created, all that is required. #### increase if I want more PN breeding ###
    NewMult_Females <- rbind.fill(NewMult_Females, (filter(NewMultPop, sex == "F" & genoA =="a/a")))
    
    NewMult_Males <- NewMultPop %>% filter(sex == "M" & genoA =="a/a") %>% top_frac(0.05, merit) ## top 10%% of new males are added to the production population
    
    MultPop <- rbind.fill(MultPop, NewMult_Females, NewMult_Males) %>% mutate(ID = as.character(ID))
    ### ### add in operator to coerce ID to be a character vector for BW antijoin below ###
    
    ## MultPop should retain some males or male piglets born to breed with the new females. Breeding in MultPop not to be AI only
    
    # if (g == 0) { 
    #  MultPopAll <- data.frame(MultPop)
    #  } else  (MultPopAll <- rbind.fill(MultPopAll, MultPop))
    
  } 
  
  ###############################
  
  if (g >= 0){
    
    MultPop <-  anti_join(MultPop, Mult_Females, by = "ID") #### keep mult pop females that have bred out of potential BW pop
    
    Mult_BW_Fem <- MultPop %>% filter(sex == "F" & herd == "B" | herd =="A") %>% filter(age >= AgeFirstMate & age %% FarrowInt == rem) %>% top_frac(0.7, merit) ### need to be creaming off the percentage and not returning!!! 
    MultPop <- anti_join(MultPop, Mult_BW_Fem, by = "ID") ## removes pigs transferred down to BWpop from MultPop permanently
    ## no A herd as they must have been bred in Mult tier ### no T herd either as the piglets are transferred directly to commercial tier?#
    
    SPF_BWpop_Males <- SPFpop %>% filter(sex == "M" & herd == "T") %>% filter(age >= AgeFirstMate) %>% top_frac(0.2, merit) #take away top 2% for SPF breeding, select next 10% for PN 
    
    SPFpop <-  anti_join(SPFpop, SPF_BWpop_Males, by = "ID") ##no returning pigs as not an AI program
    
    BWpop <- rbind.fill(BWpop, Mult_BW_Fem, SPF_BWpop_Males) ## put multipler females into BW population
    
    BW_Fem_Breed <- BWpop %>% filter(sex == "F") %>% filter(age >= AgeFirstMate & age %% FarrowInt == rem) 
    if (nrow(filter(BW_Fem_Breed %>% filter(age >= AgeFirstMate & genoA == "a/a"))) >= 50) {
      BW_Fem_Breed <- BW_Fem_Breed %>% filter(genoA == "a/a") # %>% top_frac(0.5, merit) 
    } else if (nrow(filter(BW_Fem_Breed, genoA != "A/A")) > 100) {
      BW_Fem_Breed <- BW_Fem_Breed %>% filter(genoA != "A/A")  %>% top_frac(0.5, merit) 
    } else  (BW_Fem_Breed <- BW_Fem_Breed %>% top_frac(0.05, merit))
    
    BW_Males <- BWpop %>% filter(sex == "M") %>% filter(age >= AgeFirstMate) %>% top_frac(0.5, merit) 
    
    CommercialPop <- CreatePiglets(BW_Males, BW_Fem_Breed, indexSD, g, paste0("BW",g,"_"), littersize) #piglets will be T herd, females from B herd. 
    ### maybe transfer BW males and retain some BW females for breeding if they are of better merit than piglets born??? ###
    
    NewBW_Females <- CommercialPop %>% filter(sex == "F" & genoA == "a/a") %>% top_frac(0.1, merit) ## takes the top 20% of new females created, all that is required. #### increase if I want more PN breeding ###
    
    BWpop <- rbind.fill(BWpop, NewBW_Females)
    
    ## maybe need culling function to remove low perfroming animals. May also not need population max as these animals are passed to the commercial breeders for growing.
    
    ## BWpop should retain some males or male piglets born to breed with the new females. Breeding in BWpop not to be AI only
    
    
    #######
    
    #BWpop <- rbind.fill(BWpop, BW_Fem_Breed) ###combines all piglets in or from Multpop for next generation selection
    
    #MultPop <- rbind.fill(MultPop, MultPop_Fem_Breed) ###puts breeders back into MultPop. Can't be selected above, needed in table for selection
    if (g == 0) { 
      BWpopAll <- data.frame(BWpop)
    } else  (BWpopAll <- rbind.fill(BWpopAll, BWpop))
    
    
    # if (g == 0) { 
    #    CommercialPopAll <- data.frame(CommercialPop)
    #  } else  (CommercialPopAll <- rbind.fill(CommercialPopAll, CommercialPop))  
  }
  
  
  SPFpop <- rbind.fill(SPFpop, SPF_Breeders, SPFNucA_Male, SPFNucA_Fem, SPFNucB_Male, SPFNucB_Fem, SPFNucT_Male, SPFNucT_Fem) %>% distinct() ### includes all SPF pigs, no selection. Cull after ageing function
  ## All SPF pigs. SPF pigs moved to another tier are excluded by antijoin from SPF tables
  ### remerges all SPFpops, older pigs with lower merit will not get stuck in breeding population - all SPF pigs as a single herd!!! ###
  ## need to rebind as top 2% of boars from total SPF pop will not remain the same!
  ## put breeders back with SPFpop, if they are not used for SPF pop breeding next round they may be used for PN. Need to select breeders from same pool as before ##
  
  if (sum(SPFpop$age > 8) > 9000) {
    
    NextBreeders <- SPFpop %>% filter (sex =="F" & age > 9 & age %% FarrowInt == 2) ### filters all pigs that will breed next loop to be included for selection
    
    #### removing females with neither allele for resistance from herd
    SPFpopA_females <- SPFpop %>% filter(sex == "F" & herd == "A") %>% filter(age >= 8) %>% filter (genoA != "A/A") %>% top_n(4500, merit) ### same as numbers for initial basepop
    SPFpopB_females <- SPFpop %>% filter(sex == "F" & herd == "B") %>% filter(age >= 8) %>% filter (genoA != "A/A") %>% top_n(2250, merit) ### not removing females available for breeding in the next generation??? ###
    SPFpopT_females <- SPFpop %>% filter(sex == "F" & herd == "T") %>% filter(age >= 8) %>% filter (genoA != "A/A") %>% top_n(2250, merit)
    
    
    SPFpop_males <- SPFpop %>% filter(sex == "M") %>% filter(age >= 8) %>% filter(genoA != "A/A") %>% top_n(500, merit) ### do we need these? increades gene pool for selection so maybe?
    SPFA_males <- SPFpop %>% filter(sex == "M" & herd == "A") %>% filter(age >= 8) %>% filter(genoA != "A/A") %>% top_n(100, merit)
    SPFB_males <- SPFpop %>% filter(sex == "M" & herd == "B") %>% filter(age >= 8) %>% filter(genoA != "A/A") %>% top_n(100, merit)
    SPFT_males <- SPFpop %>% filter(sex == "M" & herd == "T") %>% filter(age >= 8) %>% filter(genoA != "A/A") %>% top_n(100, merit)
    
    ### insert code to remove bottom frac of A/A pigs if getting too high population. not over 2000 piglets, 3000 breeders?
    
    SPFpop_males <- rbind.fill(SPFpop_males, SPFA_males, SPFB_males, SPFT_males)  %>% distinct ()
    SPFA_males <- NULL
    SPFB_males <- NULL
    SPFT_males <- NULL
    #SPFpopB_males <- SPFpop %>% filter(sex == "M" & herd == "B") %>% filter(age >= 8) %>% top_n(500, merit) 
    #SPFpopT_males <- SPFpop %>% filter(sex == "M", herd == "T") %>% filter(age >= 8) %>% top_n(500, merit) 
    
    #SPFpopA_piglets <- SPFpop %>% filter(herd == "A" & sex =="M") %>% filter(age < 4 & age >= 1) %>% top_n(1000, merit) 
    #SPFpopB_piglets <- SPFpop %>% filter(herd == "B" & sex =="M") %>% filter(age < 4 & age >= 1) %>% top_n(1000, merit) ### must be 1 to cull. allows for birth and genomic/pedigree selection
    #SPFpopT_piglets <- SPFpop %>% filter(herd == "T" & sex =="M") %>% filter(age < 4 & age >= 1) %>% top_n(1000, merit)
    
    topSPFpop <-  rbind.fill(NextBreeders, SPFpopA_females, SPFpopB_females, SPFpopT_females, SPFpop_males) 
    
    ##remove low performing males from SPF population 
    
    ### careful that older pigs are not killed as piglets will be higher merit but can't breed/ put in selection for if over 16 months???
    #### may need to ensure not all old pigs are culled so that there are enough older animals for breeding!
    #already only 1000 males selected.
    
    SPFpop <- semi_join(topSPFpop, SPFpop, by = "ID") %>% distinct() %>% mutate(ID = as.character(ID))### joins only the selected SPFpop number and removes any pigs duplicated. 
    
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
    SPF_Breeders <- NULL ### next generation breeders are a different cycle
  }
  
  SPFpop <- data.frame(ageing(SPFpop))  
  
  print(paste0("SPF Males:", sum(SPFpop$sex == "M"), sep = " "))  #### could be worth changing these to age > age first mate to see breeding animals???
  print(paste0 ("SPF Females:", sum(SPFpop$sex =="F"), sep = " "))
  
  print(paste0 ("SPFA_Resistance:", round(nrow(filter(SPFpop, herd == "A" & genoA == "a/a")) / nrow(filter(SPFpop, herd == "A")) *100), "%"))
  print(paste0 ("SPFB_Resistance:", round(nrow(filter(SPFpop, herd == "B" & genoA == "a/a")) / nrow(filter(SPFpop, herd == "B"))*100), "%"))
  print(paste0 ("SPFT_Resistance:", round(nrow(filter(SPFpop, herd == "T" & genoA == "a/a")) / nrow(filter(SPFpop, herd == "T")) *100), "%"))
  
  ########################
  
  ### need to have ProdPop breeders merged into daat frame before this selection step ### may be too many!!!
  if (sum(ProdPop$age > 8) > 33000){
    
    NextBreedersProd <- ProdPop %>% filter (sex =="F" & age > 9 & age %% FarrowInt == 2)
    
    ProdPopFem <- ProdPop %>% filter(sex == "F") %>% filter(age >= AgeFirstMate) %>% filter(genoA != "A/A") %>% top_n(33000, merit) ### add in selection for resistance geno here ### filter breeding animals o ensure they remain ### this code works!!!
    ProdPopMales <- ProdPop %>% filter (sex == "M") %>% filter(age > AgeFirstMate) %>% filter(genoA != "A/A") %>% top_n(800, merit)
    ProdPopPiglets <- ProdPop %>% filter(sex == "F") %>% filter(age < AgeFirstMate & age >= 1) %>% filter(genoA != "A/A") %>% top_n(15000, merit)
    ProdPop <- rbind.fill(ProdPopFem, ProdPopMales, ProdPopPiglets, NextBreedersProd) ### males for prod pop are taken each time from SPF_A. Should I set this to retain some males??
  }
  
  ProdPop <- rbind.fill(ProdPop, Prod_Females) %>% distinct() #### may need to make prod females table
  
  ProdPop <- data.frame(ageing(ProdPop)) 
  
  print(paste0("Prod Males:", sum(ProdPop$sex == "M"), sep = " "))
  print(paste0 ("Prod Females:", sum(ProdPop$sex =="F"), sep = " "))
  
  print(paste0 ("ProdPop:", round(nrow(filter(ProdPop, genoA == "a/a")) / nrow(filter(ProdPop)) *100), "%"))
  
  # ProdPop <- anti_join (ProdPop, SPFpop, by = "ID") %>% mutate(ID = as.character(ID)) ## no animals should be able to exist in both tiers. Checks females, males are not transferred as it is an AI program. 
  #############
  
  ## sows bred from already rejoined to data frame
  if (sum(MultPop$age > 8) > 80000) {  ### may be too many!!!
    MultPopFem <- MultPop %>% filter(sex == "F") %>% filter(age >= AgeFirstMate) %>% filter(genoA != "A/A") %>% top_n(70000, merit)
    MultPopMales <- MultPop %>% filter(sex == "M") %>% filter(age >= AgeFirstMate) %>% filter(genoA != "A/A")%>% top_n(10000, merit)
    MultPopPiglets <- MultPop %>% filter(sex == "F") %>% filter(age < AgeFirstMate & age >= 1) %>% filter(genoA != "A/A") %>% top_n(20000, merit) ### may be better to transfer down, let them age, and then determine whether to breed
    MultPop <- rbind.fill(MultPopFem, MultPopMales, MultPopPiglets) ### combines all piglets in or from Multpop for next generation selection
  }
  
  MultPop <- rbind.fill(Mult_Females, MultPop) %>% distinct() #### may need to make mult females table. ### remerge multpop breeders so they can be used again??
  
  MultPop <- ageing(MultPop)
  
  print(paste0("Mult Males:", sum(MultPop$sex == "M"), sep = " "))
  print(paste0 ("Mult Females:", sum(MultPop$sex =="F"), sep = " "))
  
  print(paste0 ("MultPop:", round(nrow(filter(MultPop, genoA == "a/a")) / nrow(filter(MultPop)) *100), "%"))
  
  # MultPop <- anti_join(MultPop, ProdPop, by = "ID") %>% mutate(ID = as.character(ID)) ## ensures there is no cross over between tiers
  
  if (sum(BWpop$age > 8) > 150000) {  ### may be too many!!!
    
    BWpopFem <- BWpop %>% filter(sex == "F") %>% filter(age >= AgeFirstMate) %>% filter(genoA != "A/A") %>% top_n(100000, merit)
    BWpopMales <- BWpop %>% filter(sex == "M") %>% filter(age >= AgeFirstMate) %>% filter(genoA != "A/A") %>% top_n(10000, merit)
    BWpopPiglets <- BWpop %>% filter(sex == "F") %>% filter(age < AgeFirstMate & age >= 1) %>% filter(genoA != "A/A") %>% top_n(20000, merit) ### may be better to transfer down, let them age, and then determine whether to breed
    BWpop <- rbind.fill(BWpopFem, BWpopMales, BWpopPiglets) ### combines all piglets in or from Multpop for next generation selection
  }
  
  BWpop <- ageing(BWpop)
  
  print(paste0("BW Males:", sum(BWpop$sex == "M"), sep = " "))
  print(paste0 ("BW Females:", sum(BWpop$sex =="F"), sep = " "))
  
  print(paste0 ("BWpop:", round(nrow(filter(BWpop, genoA == "a/a")) / nrow(filter(BWpop)) *100), "%"))
  
  # BWpop <- anti_join(BWpop, MultPop, by = "ID") %>% mutate(ID = as.character(ID))## ensures there is no cross over between tiers
  ### remerge BWpop breeders so they can be used again??
  
  print(paste0 ("Commercial Pigs:", nrow(CommercialPop), sep = " ")) ## need a table binding all commercial piglets to see value of IAV to farmers
  print(paste0 ("Commercial Pigs:", round(nrow(filter(CommercialPop, genoA == "a/a")) / nrow(filter(CommercialPop)) *100), "%"))
  
  
  Allpop <- rbind.fill(SPFpop, ProdPop, MultPop, BWpop) # doesn't include dying pigs and designate as dead. missing pigs! 
  ### need to assign pigs removed by joining into Allpop or another table of cast offs
  
  
  if (g == 0) { 
    Allpop2 <- data.frame(Allpop)
  } else  (Allpop2 <- rbind.fill(Allpop2, Allpop))
  
  
  ################################
  #### need data.frame compiling all commercial pigs created and distributed as piglets/breeding animals to farmers
  
} ### end of g loop

#### 

###### need to create table of all pigs that have ever been created and their genotypes/ages
###### need to create table of all pigs that have ever bred and their genotypes


