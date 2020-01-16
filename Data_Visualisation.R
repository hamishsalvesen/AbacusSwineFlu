library (ggplot2)

### graph the proportion of homozygote resistant animals per generation
##### AllPop GENO ########
 Allpop_GenoCount <- Allpop2 %>% filter(gen >= 0) %>% 
  group_by(gen, genoA) %>%
   summarise (geno_count = n()) %>%
   mutate(prop = prop.table(geno_count) *100) 
 
Allpop_GenoCount <- Allpop_GenoCount %>% filter(genoA %in% c("a/a"))
 
ggplot(Allpop_GenoCount) + 
  geom_line (mapping = aes(x = gen , y = prop, colour = 5)) +
labs (x = "Generation", y = "% resistant animals", title = "IAV resistant pigs in AllPop w sows, boars & piglets SPF Selection") +
  theme(legend.position = "none")

####### SPFpop GENO #########
SPFpop_GenoCount <- SPFpop %>%
  group_by(gen, genoA) %>%
  summarise (geno_count = n()) %>%
  mutate(prop = prop.table(geno_count) *100) 

SPFpop_GenoCount <- SPFpop_GenoCount %>% filter(genoA %in% c("a/a"))

ggplot(SPFpop_GenoCount) + 
  geom_line (mapping = aes(x = gen , y = prop, colour = 5)) +
  labs (x = "Generation", y = "% resistant animals", title = "IAV resistant pigs in SPF Tier w Selection") +
  theme(legend.position = "none")
######################################

####### ProdPop GENO #########
ProdPop_GenoCount <- ProdPop %>%
  group_by(gen, genoA) %>%
  summarise (geno_count = n()) %>%
  mutate(prop = prop.table(geno_count) *100) 

ProdPop_GenoCount <- ProdPop_GenoCount %>% filter(genoA %in% c("a/a"))

ggplot(ProdPop_GenoCount) + 
  geom_line (mapping = aes(x = gen , y = prop, colour = 5)) +
  labs (x = "Generation", y = "% resistant animals", title = "IAV resistant pigs in Prod Tier w sows, boars & piglets SPF Selection") +
  theme(legend.position = "none")
######################################

### want to graph the merit of each animal of their age as scatter plot. then create line
###### AllPop MERIT ###############
Allpop_Merit <- Allpop2 %>% filter(gen >= 0) %>%
  group_by(gen) %>%
  summarise (mean_merit = mean(merit)) 

ggplot(Allpop_Merit) + 
  geom_line (mapping = aes(x = gen , y = mean_merit, colour = 5)) +
  labs (x = "Generation", y = "Mean Genetic Merit", title = "Genetic Merit of All Pigs (GE_MI) SPF Selection sows, boars & piglets") +
  theme(legend.position = "none")

###### SPFpop MERIT ###############
SPFpop_Merit <- SPFpop %>%
  group_by(gen) %>%
  summarise (mean_merit = mean(merit)) 

ggplot(SPFpop_Merit) + 
  geom_line (mapping = aes(x = gen , y = mean_merit, colour = 5)) +
  labs (x = "Generation", y = "Mean Genetic Merit", title = "Genetic Merit of SPF Tier (GE_MI) SPF Selection sows, boars & piglets") +
  theme(legend.position = "none")

### allele count of a as a total proportion in population over time
##############################

SPFpop_GenoCount <- SPFpop %>%
  group_by(gen, genoA) %>%
  summarise (geno_count = n()) %>%
  mutate(prop = prop.table(geno_count) *100) 

SPFpop_GenoCount <- separate(SPFpop, genoA, into = c("Allele1", "Allele2"), sep = "/", remove = TRUE) %>% mutate(Allele1 = as.character(Allele1), Allele2 = as.character(Allele2))

SPFpop_GenoCount <- SPFpop_GenoCount %>%
  group_by(gen, Allele1, Allele2) %>%
  summarise (sum(Allele1 == "a"), sum (Allele2 == "a")) 

SPFpop_GenoCount <- SPFpop_GenoCount %>% filter(genoA %in% c("a/a"))

ggplot(SPFpop_GenoCount) + 
  geom_line (mapping = aes(x = gen , y = prop, colour = 5)) +
  labs (x = "Generation", y = "% resistant animals", title = "IAV resistant pigs in SPF Tier after 72 months") +
  theme(legend.position = "none")

