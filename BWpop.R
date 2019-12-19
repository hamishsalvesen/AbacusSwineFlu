MultPop <-  anti_join(MultPop, MultPop_Fem_Breed, by = "ID") #### keep mult pop females that have bred out of potential BW pop

BW_Fem <- MultPop %>% filter(sex == "F", herd == "B" | herd =="A", age >= AgeFirstMate) %>% top_frac(0.7, merit) ### need to be creaming off the percentage and not returning!!! 
## no A herd as they must have been bred in Multt tier ### no T herd either as the piglets are transferred directly to commercial tier ###

#### not selecting correctly here.... 
Mult_BW_Fem_Breed <- BW_Fem %>% filter(age %% FarrowInt == rem) ## probably best to filter for age here, not realistic to redo it for each cycle ## No BW_Fem_Breed as filtered here

SPF_BWpop_Males <- SPFpop %>% filter(sex == "M", herd == "T", age >= AgeFirstMate) %>% top_frac(0.1, merit) #take away top 2% for SPF breeding, select next 10% for PN 
#SPFpop <-  anti_join(SPFpop, SPF_MultPop_Males, by = "ID") ## returns all the pigs in SPFpop that are not present in SPF_ProdNuc_Breeders_M

BWpop <- rbind.fill(BWpop, Mult_BW_Fem_Breed) ## put multipler females into BW population

BW_Fem_Breed <- BWpop %>% filter(sex == "F", age %% FarrowInt == rem)

MultPop <- anti_join(MultPop, Mult_BW_Fem_Breed, by = "ID") ### remove BW bred females, leave rest in multPop

MultPop <- rbind.fill(MultPop, MultPop_Fem_Breed) ### rejoin bred multpop females for next generation ###

CommercialPop <- CreatePiglets(SPF_BWpop_Males, BW_Fem_Breed, indexSD, g, paste0("BW",g,"_"), littersize) #piglets will be T herd, females from B herd. 
### maybe transfer BW males and retain some BW females for breeding if they are of better merit than piglets born??? ###

## maybe need culling function to remove low perfroming animals. May also not need population max as these animals are passed to the commercial breeders for growing. 

print(paste0("BW Males:", sum(BWpop$sex == "M"), sep = " "))
print(paste0 ("BW Females:", sum(BWpop$sex =="F"), sep = " "))

print(paste0 ("Commercial Pigs:", nrow(CommercialPop), sep = " "))
