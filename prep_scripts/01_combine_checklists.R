##############################################################################
# Prepare raw list of North american macrofungal names

# VJ Svahnstrom
##############################################################################

# load packages ----
library(tidyverse)
library(readr)
library(stringr)

# load possible microfungi in MycoPortal checklist ----
microfungi <- read_csv('microfungi.csv')

# read in Mycoportal checklist, reformat ----
# From https://www.mycoportal.org/portal/checklists/checklist.php?clid=60&pid=2
myco.list <- read_csv('01_raw_data/mycoportal_checklist.csv') %>%
  mutate(Genus = str_extract(ScientificName, "\\w+"), # create genus column
         Source = 'Mycoportal', # specify source of name
         NameAuthor = ifelse(is.na(ScientificNameAuthorship), 
                             ScientificName,
                             paste(ScientificName, ScientificNameAuthorship))) %>% # paste name + author together
  filter(!grepl("\\.", ScientificName) & # removes infraspecific taxa
           !(ScientificName %in% microfungi$Species) &
           !(Genus %in% microfungi$Species))  #removes known microfungi
  

# load in download of SpeciesFungorum higher taxonomy, filtered to macrofungal groups ----
macro.taxa <- read_csv('01_raw_data/dbo_FundicClassification.csv') %>%
  filter(`Subphylum name` == 'Agaricomycotina' | `Class name` == 'Geoglossomycetes' |  `Class name` == 'Pezizomycetes' | `Order name` == 'Xylariales')

# read in protochecklist, filter to macrofungal groups, reformat ----
# https://doi.org/10.1080/00275514.2018.1515410
proto.list <- read_csv('01_raw_data/protochecklist.csv') %>%
  mutate(Genus = str_extract(MycoBank, "\\w+"), # create genus column
         Source = 'Protochecklist', # specify source of names
         NameAuthor = ifelse(is.na(Author),
                             MycoBank, 
                             paste(MycoBank, Author))) %>% # paste name + author together
  filter(!is.na('Current/Correct Name') & # remove unplaced names
           str_count(MycoBank, "\\w+") >= 2 & 
           Genus %in% macro.taxa$`Genus name` & # filter to macrofungi 
           !(MycoBank %in% checklist$ScientificName) & # remove duplicates in MycoPortal checklist
           !grepl("\\.", MycoBank)) %>% # remove infraspecifics
  rename(ScientificName = 'MycoBank',
         ScientificNameAuthorship = 'Author') %>% # rename columns to match MycoPortal checklist
  dplyr::select(c(ScientificName, ScientificNameAuthorship, Genus, NameAuthor, Source)) # select desired columns


# combine checklists ----
checklist_total <- bind_rows(myco.list, proto.list) 

# write combined checklist ----
write_csv(checklist_total, 'checklist.csv')



