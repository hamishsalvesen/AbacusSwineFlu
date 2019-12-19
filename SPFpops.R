for (g in -50:0){
  print(paste0("Gen:", g, sep = " "))
  
### select top 25% of eligible breeding females
### select males from top percentage. Should leave other breeding age pigs for below tiers
BreedSPFA_Male <- SPFpop %>% filter (sex == "M" & herd == "A" & age >= AgeFirstMate) %>% top_frac(0.02, merit) 
BreedSPFA_Fem <- filter(SPFpop, sex == "F" & herd == "A" & age >= AgeFirstMate & age %% FarrowInt == rem)  %>% top_frac(0.25, merit) ### remove these pigs from prod tier. Will leave some pigs that can be bred from!

newGenSPFNucA <- CreatePiglets(BreedSPFA_Male, BreedSPFA_Fem, indexSD, g, paste0("SPF",g,"_"), littersize)

##### separate into B herd ###########
BreedSPFB_Male <- SPFpop %>% filter (sex == "M" & herd == "B" & age >= AgeFirstMate) %>% top_frac(0.1, merit) 
BreedSPFB_Fem <- filter(SPFpop, sex == "F" & herd == "B" & age >= AgeFirstMate & age %% FarrowInt == rem) %>% top_frac(0.25, merit) ### remove these pigs from prod tier. Will leave some pigs that can be bred from!

newGenSPFNucB <- CreatePiglets(BreedSPFB_Male, BreedSPFB_Fem, indexSD, g, paste0("SPF",g,"_"), littersize)

SPFNucA_Male <- filter(newGenSPFNucA, sex == "M") %>% top_(0.02, merit)


##### separate into T herd ###########
BreedSPFT_Male <- SPFpop %>% filter (sex == "M" & herd == "T" & age >= AgeFirstMate) %>% top_frac(0.1, merit) 
BreedSPFT_Fem <- filter(SPFpop, sex == "F" & herd == "T" & age >= AgeFirstMate & age %% FarrowInt == rem ) %>% top_frac(0.25, merit) ### remove these pigs from prod tier. Will leave some pigs that can be bred from!

newGenSPFNucT <- CreatePiglets(BreedSPFT_Male, BreedSPFT_Fem, indexSD, g, paste0("SPF",g,"_"), littersize)

######### ### new SPFpop piglets join rows at the end of the code ###


print(paste0("SPF Males:", sum(SPFpop$sex == "M"), sep = " "))
print(paste0 ("SPF Females:", sum(SPFpop$sex =="F"), sep = " "))

SPF_Breeders <- rbind.fill(BreedSPFA_Fem, BreedSPFB_Fem, BreedSPFT_Fem) ### animals for SPF breeding
### will rejoin SPFpop at end of code incase they are no longer in top % for mating and can be used in below tiers  
### males are not removed as they can inseminate multiple tiers through an AI program


#####################
############################################


SPFpop <- rbind.fill(SPFpop, SPF_Breeders, newGenSPFNucT, newGenSPFNucA, newGenSPFNucB) %>% distinct() ### includes all SPF pigs, no selection. Cull after ageing function
## All SPF pigs. SPF pigs moved to another tier are excluded by antijoin from SPF tables
### remerges all SPFpops, older pigs with lower merit will not get stuck in breeding population - all SPF pigs as a single herd!!! ###
## need to rebind as top 2% of boars from total SPF pop will not remain the same!
## put breeders back with SPFpop, if they are not used for SPF pop breeding next round they may be used for PN. Need to select breeders from same pool as before ##

### population controls ###
## SPF pop retainers need to consider herd distributions with such low numbers initially breeding

if (nrow(SPFpop) > 9000) {
  NextBreeders <- SPFpop %>% filter (sex =="F" & age > 9 & age %% FarrowInt == 2) ### filters all pigs that will breed next loop to be included for selection

SPFpopA_females <- SPFpop %>% filter(sex == "F" & herd == "A") %>% filter(age >= 8) %>% top_n(4500, merit) ### same as numbers for initial basepop
SPFpopB_females <- SPFpop %>% filter(sex == "F" & herd == "B") %>% filter(age >= 8) %>% top_n(2250, merit) ### not removing females available for breeding in the next generation??? ###
SPFpopT_females <- SPFpop %>% filter(sex == "F" & herd == "T") %>% filter(age >= 8) %>% top_n(2250, merit)

SPFpop_males <- SPFpop %>% filter(sex == "M") %>% filter(age >= 8) %>% top_n(500, merit) #### may need to ensure not all old pigs are culled so that there are enough older animals for breeding!
SPFpopB_males <- SPFpop %>% filter(sex == "M" & herd == "B") %>% filter(age >= 8) %>% top_n(500, merit) 
SPFpopT_males <- SPFpop %>% filter(sex == "M", herd == "T") %>% filter(age >= 8) %>% top_n(500, merit) 

#SPFpopA_piglets <- SPFpop %>% filter(herd == "A" & sex =="M") %>% filter(age < 4 & age >= 1) %>% top_n(1000, merit) 
#SPFpopB_piglets <- SPFpop %>% filter(herd == "B" & sex =="M") %>% filter(age < 4 & age >= 1) %>% top_n(1000, merit) ### must be 1 to cull. allows for birth and genomic/pedigree selection
#SPFpopT_piglets <- SPFpop %>% filter(herd == "T" & sex =="M") %>% filter(age < 4 & age >= 1) %>% top_n(1000, merit)

topSPFpop <-  rbind.fill(NextBreeders, SPFpopA_females, SPFpopB_females, SPFpopT_females, SPFpop_males) ## remove low performing males from SPF population 


### careful that older pigs are not killed as piglets will be higher merit but can't breed/ put in selection for if over 16 months???
#### may need to ensure not all old pigs are culled so that there are enough older animals for breeding!
#already only 1000 males selected.

SPFpop <- semi_join(topSPFpop, SPFpop, by = "ID") %>% distinct() ### joins only the selected SPFpop number and removes any pigs duplicated. 


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

}