SPFpop <-  anti_join(SPFpop, SPF_Breeders, by = "ID") ## removes breeding SPF breeding animals from SPF pop,
# selection will exclude nucleus breeding animals ## keeps out for mult & BW tier 

SPF_ProdPop_Fem <- SPFpop %>% filter(sex == "F", herd == "A", age >= AgeFirstMate) %>% top_frac(0.25, merit) #takes away top 25% for PN 
ProdPop_Fem_Breed <- filter(SPF_ProdPop_Fem, age %% FarrowInt == rem) #already filtered on age and merit

SPF_ProdPop_Males <- SPFpop %>% filter(sex == "M", herd == "A", age >= AgeFirstMate) %>% top_frac(0.1, merit) #takes top 10% available boars
#SPFpop <-  anti_join(SPFpop, SPF_ProdPop_Males, by = "ID") ## removes prod pop males from the SPFpop #can't be reused in SPFpop ### don't need as AI is performed on SPF pops

SPFpop <-  anti_join(SPFpop, ProdPop_Fem_Breed, by = "ID") ## SPF pop loses the animals transferred to prod pop. Not reintegrated and cannot be used in SPF in future 

ProdPop <- rbind.fill(ProdPop, ProdPop_Fem_Breed) ### No males from SPF stored as used by AI. Only PN bred males will be in PN ######

Prod_Males <- rbind.fill(ProdPop, SPF_ProdPop_Males) %>% filter(sex == "M") # selects for males from ProdPop and SPF Nucleus
Prod_Females <- ProdPop %>% filter(sex == "F", age %% FarrowInt == rem) ### already over age first mate, ensure they are correct farrowing interval. 

NewProdPop <- CreatePiglets(Prod_Males, Prod_Females, indexSD, g, paste0("Prod",g,"_"), littersize)

NewProd_Females <- NewProdPop %>% filter(sex == "F") %>% top_frac(0.2, merit) ## takes the top 20% of new females created, all that is required. #### increase if I want more PN breeding ###
NewProd_Males <- NewProdPop %>% filter(sex == "M") %>% top_frac(0.05, merit) ## top 10%% of new males are added to the production population

ProdPop <- rbind.fill(ProdPop, NewProd_Females, NewProd_Males) ### combine all animals to stay in ProdPop ### prod_Females must be added in after MultPophas selected ##

print(paste0("Prod Males:", sum(ProdPop$sex == "M"), sep = " "))
print(paste0 ("Prod Females:", sum(ProdPop$sex =="F"), sep = " "))
