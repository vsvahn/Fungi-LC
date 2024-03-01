##############################################################################
# Match and filter names using GBIF taxonomic backbone

# VJ Svahnstrom
##############################################################################

# load packages ----
library(tidyverse)
library(readr)
library(rgbif)

## read in species of interest ----

# red listed species
redlist <- inner_join(read_csv('C:/Users/vsv10kg/Documents/Files_for_fungi_repo/01_raw_data/Redlist_download_NorthAmerica_2.11.23/assessments.csv'),
                                  read_csv('C:/Users/vsv10kg/Documents/Files_for_fungi_repo/01_raw_data//Redlist_download_NorthAmerica_2.11.23/taxonomy.csv'),
                                  by = "scientificName") %>%
  mutate(NameAuthor = paste(scientificName, authority), 
         Threatened = ifelse(redlistCategory %in% c('Least Concern', 'Near Threatened'), 0, 1), #  add binary threat status column
         Source = 'RedList') %>%
  filter(redlistCategory != 'Data Deficient') # remove DD species

# NatureServe assessed species
load('C:/Users/vsv10kg/Documents/Files_for_fungi_repo/01_raw_data/Fungi_Natureserve_Assessments.RData')
# filter out infraspecifics
natserv <- natserv %>% 
  mutate(NameAuthor = paste(scientificName, scientificNameAuthor),
         Threatened = ifelse(roundedGRank %in% c('G1', 'G2', 'G3', 'GH'), 1, 0), #  add binary threat status column
         Source = 'NatureServe') %>%
  filter(!(scientificName %in% redlist$scientificName))

# combine assessed species 
training.sp <- bind_rows(redlist, natserv) %>%
  select(c(NameAuthor, Source, Threatened)) %>%
  filter(!(NameAuthor %in% c('Ramaria largentii Marr & Stuntz', 'Montagnea radiosa (Pall.) Šebek'))) # remove species with known issues)

# Unassessed North american species
checklist <- read_csv('C:/Users/vsv10kg/Documents/Files_for_fungi_repo/02_cleaned_data/checklist.csv') %>%
  filter(!(NameAuthor %in% training.sp$NameAuthor)) # remove assessed sp.


## function to find accepted names using gbif name search 
# note column of names + authority should be called NameAuthor 
find_accepted_names <- function(names_df){
  
  # run the gbif name search for the input names - get the usage key from GBIF
  batch_out <-  purrr::map_dfr(names_df$NameAuthor, name_search_gbif)
  
  # filter out names with multiple matches
  grouped_names <- batch_out %>% group_by(searchName) %>%
    count() %>%
    filter(n !=1) 
  # filter out names not matched with 100% confidence
  batch_out_unique <- batch_out %>%
    filter(!(searchName %in% grouped_names$searchName),
           confidence == 100) 
  
  # now use that key to get the higher level taxonomy and taxonomic status
  # create empty df
  gbif_usage <- data.frame()
  
  # use map_df to iterate over rows and bind rows
  gbif_usage <- purrr::map_df(batch_out_unique$usageKey, ~name_usage(key = .)$data)
  
  # left join with batch_out_unique
  gbif_usage <- left_join(gbif_usage, batch_out_unique[,c('searchName', 'scientificName')], by='scientificName') %>%
    distinct(scientificName, .keep_all = T)
  return(gbif_usage)
  
}

## find accepted names, filter, and save ----

# assessed species
accepted.training <- find_accepted_names(training.sp) 
# filter out species assessed under a heterotypic synonym of accepted name, doubtful names
accepted.training <- accepted.training %>%
  filter(!(taxonomicStatus %in% c('DOUBTFUL', 'HETEROTYPIC_SYNONYM')) &
           !(species %in% c('Hygrophorus camarophyllus', 'Turbinellus floccosus',
                            'Cortinarius sublanatus', 'Lactarius silviae', 'Russula adusta'))) %>% # remove names assessed under heterotypic syn
  mutate(accepted = ifelse(is.na(accepted), scientificName, accepted), 
         acceptedKey = ifelse(is.na(acceptedKey), key, acceptedKey)) # put accepted name and key in accepted columns
# save assessed species for training
training.sp <- training.sp %>%
  right_join(accepted.training[,c('searchName', 'accepted', 'species')], join_by('NameAuthor' == 'searchName')) %>%
  rename(acceptedNameAuthor = accepted, acceptedName = species)
write_csv(training.sp, 'training_sp_filtered.csv')


# unassessed species
accepted_names_checklist <- find_accepted_names(checklist)
# filter out doubtful names
accepted_names_checklist2 <- accepted_names_checklist %>%
  filter(!(taxonomicStatus %in% c('DOUBTFUL'))) %>%
  mutate(accepted = ifelse(is.na(accepted), scientificName, accepted),
         acceptedKey = ifelse(is.na(acceptedKey), key, acceptedKey))
# save unassessed species 
checklist <- checklist %>% 
  right_join(accepted_names_checklist[,c('searchName','species', 'accepted', 'acceptedKey')],
             join_by('NameAuthor' == 'searchName')) %>%
  rename(acceptedName = species,
         acceptedNameAuthor = accepted,
         gbifAcceptedKey = acceptedKey) %>%
  group_by(acceptedName)
write_csv(checklist, 'checklist_filtered.csv')

#######################################################################

##function to find synonyms of accepted names
# use output of find_accepted_names
find_synonyms <- function(accepted_names_output){

accepted_names_output <- distinct(accepted_names_output, acceptedKey)
  
synonyms <- purrr::map_df(accepted_names_output$acceptedKey, ~name_usage(key = ., data = 'synonyms')$data)
  
# filter out synonyms matched to multiple accepted names
grouped_names <- synonyms %>% group_by(scientificName) %>%
  count() %>%
  filter(n !=1) 

# remove names matching to >1 accepted name
synonyms <- synonyms %>%
  filter(!(scientificName %in% grouped_names$scientificName))

return(syn)
}

## find synonyms, reformat and save data ----

# assessed species
training.synonyms <- find_synonyms(accepted.training)
training.synonyms <- training.synonyms %>% 
  bind_rows(accepted.training[,c('species','accepted', 'acceptedKey')]) %>%
  dplyr::select(c(scientificName, species, accepted, acceptedKey)) %>%
  mutate(scientificName = ifelse(is.na(scientificName), accepted, scientificName)) %>%
  rename(acceptedName = species,
         acceptedNameAuthor = accepted,
         gbifAcceptedKey = acceptedKey,
         NameAuthor = scientificName) %>%
  filter(!grepl("\\?", NameAuthor)) %>% # remove mistakes in gbif
  arrange(acceptedName)

# unassessed species
checklist.synonyms <- find_synonyms(accepted_names_checklist)
checklist.synonyms2 <- checklist.synonyms %>%
  bind_rows(accepted_names_checklist[,c('species','accepted', 'acceptedKey')]) %>% # add row for accepted name
  dplyr::select(c(scientificName, species, accepted, acceptedKey)) %>%
  mutate(scientificName = ifelse(is.na(scientificName), accepted, scientificName)) %>%
  rename(acceptedName = species,
         acceptedNameAuthor = accepted,
         gbifAcceptedKey = acceptedKey,
         NameAuthor = scientificName) %>%
  filter(!grepl("\\?", NameAuthor)) %>% # remove mistakes in gbif
  distinct(.keep_all = T) %>%
  arrange(acceptedName)

## combine all synonyms and save ----
all_names <- rbind(training.synonyms, checklist.synonyms) %>%
  distinct(.keep_all = T) %>% #removes duplicates 
  mutate(mycoportalNameAuthor =  ifelse(NameAuthor == 'Cantharellus enelensis Voitk, Thorn, Lebeuf & J.I.Kim',  'Cantharellus enelensis Voitk, Thorn, Lebeuf, J.I. Kim',
                                        str_replace_all(NameAuthor, "\\.(?=[A-Z]\\p{Ll})", ". "))) %>% # for matching names in MycoPortal
  group_by(NameAuthor) %>% filter(n() == 1) %>% # remove names matching to >1 accepted name
  ungroup()
write_csv(all_names, 'all_names_synonyms.csv')
