Prod_MultPop_Females <- ProdPop %>% filter(sex == "F", herd == "A" | herd == "B", age >= AgeFirstMate) %>% top_frac(0.5, merit) #selects top 50% for mult tier from Prod. 
MultPop_Fem_Breed <- filter(Prod_MultPop_Females, age %% FarrowInt == rem) #already filtered on age and merit

SPF_MultPop_Males <- SPFpop %>% filter(sex == "M", herd == "B", age >= AgeFirstMate) %>% top_frac(0.1, merit) #takes top 10% available boars, SPF pop excluded.## breeders removed above
#SPFpop <-  anti_join(SPFpop, SPF_MultPop_Males, by = "ID") ## removes mult pigs from SPF pop 

ProdPop <- anti_join(ProdPop, MultPop_Fem_Breed, by = "ID") ## removes pigs transferred down to multpop from ProdPop permanently

MultPop <- rbind.fill(MultPop, MultPop_Fem_Breed, SPF_MultPop_Males) #puts new Multpop with new piglets from SPF & Prod. Are these females now permanently in this tier?

Mult_Males <- MultPop %>% filter(sex == "M") ### are these required ??? #### maybe to keep the pigs that have been put into these herds.. ###
Mult_Females <- MultPop %>% filter(sex == "F" & age %% FarrowInt == rem) ## probably best to filter for age here, not realistic to redo it for each cycle ##

NewMultPop <- CreatePiglets(Mult_Males, Mult_Females, indexSD, g, paste0("Mult",g,"_"), littersize)

MultPop <- rbind.fill(MultPop, NewMultPop, MultPop_Fem_Breed)
MultPop <- MultPop %>% filter(sex == "F")


print(paste0("Mult Males:", sum(MultPop$sex == "M"), sep = " "))
print(paste0 ("Mult Females:", sum(MultPop$sex =="F"), sep = " "))
