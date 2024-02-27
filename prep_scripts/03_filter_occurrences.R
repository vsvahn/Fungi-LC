##############################################################################
# Filter records to North American macrofungi

# VJ Svahnstrom
##############################################################################

# load packages ----
library(tidyverse)
library(readr)

## read in names ---- 

# read in list of all names associated with accepted names
all_names <- read_csv('all_names_synonyms.csv') 

#####################################################################################

## read in Mycoportal files ----

# find names of all myco occs files
csv_files <- list.files(path = "01_raw_data/02_occurrence_records/mycoportal", 
                        pattern = ".csv", 
                        full.names = TRUE,
                        recursive = TRUE)

# USA dats
# remove unneccessary columns
occs.myco.usa <- read_tsv(csv_files[1]
                          #locale = readr::locale(encoding ="ISO-8859-1")
) %>%
  dplyr::select(-c(institutionCode, collectionCode, ownerInstitutionCode, catalogNumber,
                   otherCatalogNumbers, higherClassification, subgenus, verbatimTaxonRank, identifiedBy,
                   dateIdentified, identificationRemarks, taxonRemarks, identificationQualifier, typeStatus,
                   eventDate, startDayOfYear, endDayOfYear, verbatimEventDate, occurrenceRemarks, fieldNumber, eventID,
                   dataGeneralizations, dynamicProperties, associatedOccurrences, associatedSequences, 
                   reproductiveCondition, establishmentMeans, lifeStage, sex, preparations, locationID,
                   waterBody, islandGroup, island, locationRemarks, verbatimCoordinates, georeferencedBy,
                   georeferenceProtocol, georeferenceSources, georeferenceVerificationStatus, georeferenceRemarks,
                   minimumDepthInMeters, maximumDepthInMeters, verbatimDepth, verbatimElevation, disposition,
                   language, modified, references, eventDate, identificationReferences, recordEnteredBy,
                   verbatimAttributes, cultivationStatus, 'sourcePrimaryKey-dbpk', eventDate2, substrate,
                   individualCount, localitySecurity, localitySecurityReason, collID))

# Canada data
# remove unneccessary columns
occs.myco.can <- read_csv(csv_files[2],
                          #locale = readr::locale(encoding ="ISO-8859-1"), 
                          col_types = sapply(occs.myco.usa, class)) %>%
  dplyr::select(-c(institutionCode, collectionCode, ownerInstitutionCode, collectionID, catalogNumber,
                   otherCatalogNumbers, higherClassification, subgenus, verbatimTaxonRank, identifiedBy,
                   dateIdentified, identificationRemarks, taxonRemarks, identificationQualifier, typeStatus,
                   eventDate, startDayOfYear, endDayOfYear, verbatimEventDate, occurrenceRemarks, fieldNumber, eventID,
                   dataGeneralizations, dynamicProperties, associatedOccurrences, associatedSequences, 
                   reproductiveCondition, establishmentMeans, lifeStage, sex, preparations, locationID,
                   waterBody, islandGroup, island, locationRemarks, verbatimCoordinates, georeferencedBy,
                   georeferenceProtocol, georeferenceSources, georeferenceVerificationStatus, georeferenceRemarks,
                   minimumDepthInMeters, maximumDepthInMeters, verbatimDepth, verbatimElevation, disposition,
                   language, modified, rights, rightsHolder, accessRights, references, eventDate,
                   identificationReferences,  recordEnteredBy, informationWithheld))


# combine mycoportal data
# remove duplicated records (mycoportal issue)
myco.occs <- bind_rows(occs.myco.can, occs.myco.usa) %>%
  distinct(id, .keep_all = T)
rm(occs.myco.can, occs.myco.usa) # remove big files

#################################################################################################

## filter species - MycoPortal ----

# 1. Extract records for taxa with taxonID in the Mycportal checklist
# picks up some misspelled, etc. records not captured by name-matching below
# read in checklist 
checklist <- read.csv('checklist_filtered.csv') %>%
  dplyr::select(c(TaxonId,acceptedNameAuthor, acceptedName)) %>%
  drop_na(TaxonId)
myco.filtered <- myco.occs %>%
  filter(!is.na(taxonID),
         taxonID %in% checklist$TaxonId) %>% # taxonID in checklist
  inner_join(checklist, join_by(taxonID == TaxonId)) %>% # find acceptedname
  mutate(NameAuthor = paste(genus, specificEpithet, scientificNameAuthorship)) # allow binding with below occs

myco.filtered2 <- myco.occs %>%
  mutate(NameAuthor = paste(scientificName, scientificNameAuthorship)) %>% # for matching with backbone
  filter(NameAuthor %in% all_names$mycoportalNameAuthor | NameAuthor %in% all_names$NameAuthor) %>% #filter
  left_join(select(all_names, -c(NameAuthor, gbifAcceptedKey)), join_by('NameAuthor'=='mycoportalNameAuthor')) # add acceptedName

# merge filtered data back together
all.filtered.myco <- distinct(rbind(myco.filtered, myco.filtered2))

# remove duplicated records
all.filtered.myco <- distinct(all.filtered.myco, id, .keep_all = T)
rm(myco.filtered, myco.filtered2, myco.filtered3) # remove big files


#####################################################################################

## read in GBIF files ---- 

# get keys for filtering
gbifKeys <- all_names %>% distinct(gbifAcceptedKey, acceptedName, acceptedNameAuthor)

# Ascos download
# DOI 10.15468/dl.bebqbr
ascos <- read_tsv('C:/Users/vsv10kg/Downloads/0000493-240202131308920/0000493-240202131308920.csv') %>% 
  filter(speciesKey %in% gbifKeys$gbifAcceptedKey) %>% # filter to our species of interest
  left_join(gbifKeys, join_by('speciesKey'=='gbifAcceptedKey')) %>% # get accepted name
  rename(NameAuthor = scientificName, scientificName = species) %>% # reformat for rbinding 
  filter(!(countryCode %in% c('US', 'CA')), # filter to records OUTSIDE US/CA
         occurrenceStatus != 'ABSENT') %>% # remove absence records
  select(-c(datasetKey, publishingOrgKey, depth, depthAccuracy, eventDate, institutionCode,
            collectionCode, catalogNumber, identifiedBy, dateIdentified, license, rightsHolder, recordedBy,
            typeStatus, establishmentMeans, lastInterpreted, mediaType, eventDate)) # remove unneccessary cols

#basidios download
# DOI 10.15468/dl.seszkv
basidios <- read_tsv('C:/Users/vsv10kg/Downloads/0000299-240202131308920/0000299-240202131308920.csv') %>%
  filter(speciesKey %in% gbifKeys$gbifAcceptedKey) %>% # filter to our species of interest
  left_join(gbifKeys, join_by('speciesKey'=='gbifAcceptedKey')) %>% # get accepted name
  rename(NameAuthor = scientificName, scientificName = species) %>% # reformat for rbinding 
  filter(!(countryCode %in% c('US', 'CA')), # filter to records OUTSIDE US/CA
         occurrenceStatus != 'ABSENT') %>%  #remove absence records
  select(-c(datasetKey, publishingOrgKey, depth, depthAccuracy, eventDate, institutionCode,
            collectionCode, catalogNumber, identifiedBy, dateIdentified, license, rightsHolder, recordedBy,
            typeStatus, establishmentMeans, lastInterpreted, mediaType, eventDate)) # remove unneccessary cols


gbif <- rbind(ascos, basidios)
rm(basidios, ascos) # remove big files


#####################################################################################

## combining gbif and mycoportal data ---- 

# reformat for rbinding 
gbif$year <- as.numeric(gbif$year)

# bind all records together
all.occs <- bind_rows(all.filtered.myco, gbif) 
# add unique identifier for each record
all.occs <- all.occs %>% mutate(Source = ifelse(is.na(gbifID), 'MycoPortal', 'GBIF'), 
                                uniqueCombinedID = ifelse(is.na(id), paste0(gbifID,'g'),
                                                          paste0(id,'m'))) # give each occ a unique ID

# add ISO code to mycoportal records
all.occs <- all.occs %>%
  mutate(countryCode = ifelse(!is.na(countryCode), countryCode,
                              ifelse(country %in% c('Canada', 'canada', 'CANADA'), 'CA',
                                     ifelse(country %in% c('United States', 'USA', 'US', 'U.S.A.',
                                                           'United States of America', 'usa'),
                                            'US', NA))))

# remove species with no North American records 
species.remove <- setdiff(gbif$acceptedNameAuthor, all.filtered.myco$acceptedNameAuthor)
all.occs <- all.occs %>% filter(!(acceptedNameAuthor %in% species.remove))

# write all filtered data as csv
write_csv(all.occs, 'all_occs_filtered.csv')

