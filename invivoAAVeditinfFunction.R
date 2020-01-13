edit_genes <- function(popdata) { 
 
  A_ToEdit <- filter(popdata, genoA %in% c("A/A", "A/a", "a/A"))
  A_Others <- filter (popdata , genoA == "a/a")
  
  ## Change to piglets for function format
  if (nrow(A_ToEdit) > 1){

A_ToEdit <- separate(A_ToEdit, genoA, into = c("Allele1", "Allele2"), sep = "/", remove = TRUE) %>% mutate(Allele1 = as.character(Allele1), Allele2 = as.character(Allele2))
  
#  SELECTION OF PIGS ON MERIT ## would need to be on sires and dams the selection, not on the progeny
# A_ToEdit <- arrange(A_ToEdit, desc(merit)) #%>% top_frac(0.5) ## selection criteria for pigs to be edited. May need to be applied to males or females
# filter out top_frac if only certain animals are to be edited? filter sires/dams or piglets born?

### need to split table as selection of allele2 must draw from animals with big or little A/a
#for 70% of alleles, convert to little a resistance gene
A_ToEditA1 <- A_ToEdit  %>% filter(Allele1 == "A") %>% sample_frac(0.7, replace = FALSE)## need to select at the same time as genes are in the same pig vessel
A_ToEditA1$Allele1 <- tolower(A_ToEditA1$Allele1) ## use pipe for putting through and avoid splitting table?
A_ToEditA1$Allele2 <- NULL
A_ToEditA2 <- A_ToEdit  %>% filter(Allele2 == "A") %>% sample_frac(0.7, replace = FALSE) ## need to select at the same time as genes are in the same pig vessel
A_ToEditA2$Allele2 <- tolower(A_ToEditA2$Allele2) ## use pipe for putting through and avoid splitting table?
A_ToEditA2$Allele1 <- NULL

#A_ToEditA1 <- arrange(A_ToEditA1, desc(ID)) 
#A_ToEditA2 <- arrange(A_ToEditA2, desc(ID)) 

A_ToEditJoinA1 <- A_ToEditA1 %>% dplyr::select(ID, Allele1)
A_ToEditJoinA2 <- A_ToEditA2 %>% dplyr::select(ID, Allele2)

A_ToEdit$Allele1 <- NULL ## remove previous Allele1 & Allele2 Columns #pipe and select(-allele1, allele2) would work
A_ToEdit$Allele2 <- NULL

## Bind Alleles back together by ID match. replace NA (non-succesfully edited alleles) with A. 
A_ToEdit <- A_ToEdit %>% left_join (A_ToEditJoinA1, by = "ID") %>% left_join(A_ToEditJoinA2, by = "ID") 
A_ToEdit [c("Allele1", "Allele2")][is.na(A_ToEdit[c("Allele1", "Allele2")])] <- "A" ### 

### Both alleles put back together in a single column
A_ToEdit <- A_ToEdit %>% unite(col = genoA, c(Allele1, Allele2), sep = "/", remove = TRUE) ## remerge alleles columns ### need to merge the column of Allele2 data to allele1

### extra column present in new data for some reason

#A_ToEditJoin <- merge(A_ToEditA1, A_ToEditA2, by = "ID", no.dups = TRUE, all = TRUE) #### faulty merge function

popdata <- rbind(A_ToEdit, A_Others) ### add in edited animals and overwrite previous unedited rows. Recreate data frame with succuseful editing, unsuccesful editing

### need mortality rate from editing process. Use Ageing code function as frame

### better if it appears as a proportion %%%

} ## end of if function
  
return (popdata)  } 
  
  