##############################################################################
# Calculate occurrence-based predictors

# VJ Svahnstrom
##############################################################################


# read in packages
library(tidyverse)
library(rCAT)

##############################################################################

## occurrence-based predictors ----

#read in all occurrence data, filter out those with coordinate issues
occs <- read_csv('02_cleaned_data/all_occs_filtered_cleaned.csv')
occs <- occs %>% subset(georeferenced == F | (cc_equ==T & cc_zero==T & clean_tests==T & cc_err==T | is.na(cc_err)))

 
# filter duplicated georeferenced records rounded to 3 dec. places
geo.occs <- occs %>%
  filter(georeferenced==T) %>% 
  mutate(roundDecimalLatitude = round(decimalLatitude, digits = 3),
         roundDecimalLongitude = round(decimalLongitude, digits = 3)) %>% 
  distinct(acceptedNameAuthor, roundDecimalLatitude, roundDecimalLongitude, .keep_all = T)

#filter nongeoreferenced records from same locality per species 
clean.occs <- occs %>%
  filter(!is.na(locality) & georeferenced==F) %>%
  distinct(acceptedNameAuthor, locality, .keep_all = TRUE) %>%
  bind_rows(geo.occs) %>%
  bind_rows(occs %>% filter(is.na(locality) &  georeferenced==F))

# calculate specimen count 
spec.count <- clean.occs %>%
  count(acceptedNameAuthor, name = 'specimenCount')

## Calculate number of 'contemporary' (>1992) records
contemp_records <- occs %>%
  distinct(acceptedNameAuthor, locality, year, .keep_all = T) %>%
  distinct(acceptedNameAuthor, decimalLatitude, decimalLongitude, year, .keep_all=T) %>%
  group_by(acceptedNameAuthor) %>%
  summarize(propContemporary = mean(year > 1992, na.rm = T)) %>%
  mutate(propContemporary = ifelse(is.na(propContemporary), 0, propContemporary))

# calculate final variables
spec.count <- spec.count %>%
  left_join(contemp_records, by='acceptedNameAuthor') %>%
  mutate(logSpecimenCount = log(specimenCount)) 


#################################################################

## year of publication ----


# read in list of all names associated with accepted names
all_names <- read.csv('02_cleaned_data/accepted_names_synonyms/all_names_synonyms.csv') %>% dplyr::select(-c(NameAuthor, mycoportalNameAuthor, gbifAcceptedKey))

# read in Species Fungorum downloaded backbone
backbone <- read.csv('02_cleaned_data/species_fungorum_backbone.csv', encoding = "UTF-8") %>%
  dplyr::select(-c(X.1)) %>%# remove pesky row index column  
  filter(!is.na(AUTHORS)) %>% # remove names without authority
  mutate(CURRENT.NAME = ifelse(NAME.OF.FUNGUS == 'Amylocystis lapponicus', 'Amylocystis lapponica', CURRENT.NAME),
         CURRENT.NAME = ifelse(NAME.OF.FUNGUS == 'Agaricus pattersonae', 'Agaricus pattersoniae', CURRENT.NAME)) # fix known orthography errors

# find year of publication for basionyms
basionyms <- backbone %>%
  filter(RECORD.NUMBER == BASIONYM.RECORD.NUMBER) %>% # filter to basionyms
  group_by(CURRENT.NAME) %>% 
  filter(YEAR.OF.PUBLICATION == min(YEAR.OF.PUBLICATION, na.rm = T))

# find year of publication of basionyms for species in list
year <- backbone %>%
  filter(CURRENT.NAME %in% all_names$acceptedNameAuthor &
           (RECORD.NUMBER == CURRENT.NAME.RECORD.NUMBER |
              NAME.OF.FUNGUS == CURRENT.NAME)) %>% 
  left_join(basionyms, join_by(CURRENT.NAME)) %>%
  mutate(yearPublished = YEAR.OF.PUBLICATION.y) %>%
  dplyr::select(c(CURRENT.NAME, yearPublished)) %>%
  rename(acceptedNameAuthor = CURRENT.NAME)  %>% 
  distinct(.keep_all = T) 

# find species missed
missed <- all_names %>% filter(!(acceptedNameAuthor %in% year.pub$acceptedNameAuthor)) %>% distinct(acceptedNameAuthor)
# find year of publication for basionyms of unplaced species
basionyms2 <- backbone %>%
  filter(RECORD.NUMBER == BASIONYM.RECORD.NUMBER) %>%
  group_by(NAME.OF.FUNGUS) %>% 
  filter(YEAR.OF.PUBLICATION == min(YEAR.OF.PUBLICATION, na.rm = T))

# find year of publication of basionyms of unplcaed species
year2 <- backbone %>%
  filter(NAME.OF.FUNGUS %in% missed$acceptedNameAuthor) %>%
  left_join(basionyms2, join_by(NAME.OF.FUNGUS)) %>%
  mutate(yearPublished = YEAR.OF.PUBLICATION.y) %>%
  dplyr::select(c(NAME.OF.FUNGUS, yearPublished)) %>%
  rename(acceptedNameAuthor = NAME.OF.FUNGUS) %>% 
  distinct(.keep_all = T) 

# add assessed species missing from species fungorum
manual_years <- data.frame(
  acceptedNameAuthor = c('Sarcodon joeides', 'Gomphidius smithii', 'Elaphomyces decipiens'),
  yearPublished = c(1872, 1949, 1831)
)

#combine names matched to years by different methods
year.all <- rbind(year, year2, manual_years) %>%
  drop_na() 

#####################################################################

# combine predictors

# load in threat status of training species
threat <- read_csv('02_cleaned_data/training_sp_filtered.csv') %>%
  select(c(acceptedNameAuthor, Threatened)) 

# combine predictors
predictors <- left_join(spec.count, year, by='acceptedNameAuthor') %>%
  drop_na() # drop species with missing values
predictors <- left_join(sp_predictors, threat, by='acceptedNameAuthor')

# save predictor data
write_csv(predictors, '02_cleaned_data/predictors.csv')
