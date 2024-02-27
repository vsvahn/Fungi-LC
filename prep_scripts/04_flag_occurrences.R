##############################################################################
# Occurrence record processing 
# Flag erroneous records, georeference using matched localities

# VJ Svahnstrom
##############################################################################

# load packages ----
library(tidyverse)
library(CoordinateCleaner)
library(readr)

# import occurrence records for processing ----

# Data filtered to target species
all.occs <- read_csv('02_cleaned_data/occurrence_records/all_occs_filtered_matched.csv')

# remove NAs to run coordinate cleaner
occs_ll <- all.occs %>%
  drop_na("decimalLongitude", "decimalLatitude")  %>%
  mutate(georeferenced = TRUE)

# rogue single record removed due to invalid coordinates
occs_error <- occs_ll %>% filter(id == 3758416)
occs_ll <- occs_ll %>% filter(!(id %in% occs_error$id))

# 1. Flag zero lat long ----
clean_cc_zero <-
  occs_ll %>% # 1  remove zero lat long
  cc_zero(
    lon = "decimalLongitude",
    lat = "decimalLatitude",
    verbose = FALSE,
    value = "flagged"
  )
clean_cc_zero_df <- data.frame("cc_zero" = clean_cc_zero)
clean_cc_zero <- bind_cols(occs_ll, clean_cc_zero_df)

# 2. Flag records with Identical lat/lon ----
clean_cc_equ <-
  cc_equ(
    clean_cc_zero,
    lat = "decimalLatitude",
    lon = "decimalLongitude",
    verbose = FALSE,
    value = "flagged"
  )
clean_cc_equ_df <- data.frame("cc_equ" = clean_cc_equ)
clean_cc_equ <- bind_cols(clean_cc_zero, clean_cc_equ_df)

# 3. flag country centroids, identified institutions  ----
# NOTE: not checking institutions, capitals because most records for these are good
clean_cc_tests <-
  clean_coordinates(
    clean_cc_equ,
    lat = "decimalLatitude",
    lon = "decimalLongitude",
    species = "scientificName",
    tests = c("centroids"),
    centroids_detail = 'country',
    value = "flagged"
  )
clean_cc_tests_df <- data.frame("clean_tests" = clean_cc_tests)
clean_cc_tests <- bind_cols(clean_cc_equ, clean_cc_tests_df)


# 4. flag based on coordinate uncertainty field ----
clean_cc_err <-
  clean_cc_tests %>% 
  mutate(cc_err = ifelse(coordinateUncertaintyInMeters >= 100000, "FALSE", "TRUE")) 


# Results ----
# simple results table to show result of data cleaning

cc_0 = clean_cc_err %>% filter(cc_zero == "FALSE") %>% nrow()
cc_equ = clean_cc_err %>% filter(cc_equ == "FALSE") %>% nrow()
cc_clean = clean_cc_err %>% filter(clean_tests == "FALSE") %>% nrow()
cc_err = clean_cc_err %>% filter(cc_err == "FALSE") %>% nrow()
# make the final table
cleaning_step = c("cc_zero", "cc_equ", "cc_clean", "cc_err")
records_cleaned = c(cc_0, cc_equ, cc_clean, cc_err)
results = data.frame(cleaning_step, records_cleaned)


## recombine (flagged) georeferenced and nongeoreferenced records ----
flagged.occs <- all.occs %>%
  filter(is.na(decimalLatitude) | is.na(decimalLongitude)) %>%
  mutate(cc_zero = NA, cc_equ = NA, clean_tests = NA, cc_err = NA, georeferenced = FALSE) %>%
  rbind(clean_cc_err) 

## georeference records based on matched localities ----

# function that adds georeferences based on state-county-locality matches
match_locality <- function(occs) {
  # subset georeferenced records used for matching
  geo.for.match <- occs %>% 
    filter(georeferenced == T & locality != '' & !is.na(locality) & !is.na(stateProvince)) %>%  # useful records for matching
    distinct(countryCode, stateProvince, county, locality, .keep_all = T) %>%
    dplyr::select(c(countryCode, stateProvince, county, locality, decimalLatitude, decimalLongitude)) 
  # subset nongeoreferenced records with localities
  nongeo <- occs %>%
    filter(is.na(decimalLatitude) & locality != '' & !is.na(locality))
  
  # subset records which won't be touched 
  geo <- setdiff(occs, nongeo) %>%
    mutate(coordsAdded = F)
  
  # match locality 
  nongeo.matched <- left_join(nongeo, geo.for.match, by = c('countryCode', 'stateProvince','county', 'locality'), suffix = c('', '.y'))
  # reformat
  nongeo.matched <- nongeo.matched %>%
    mutate(decimalLatitude = ifelse(is.na(decimalLatitude) | decimalLatitude == "", decimalLatitude.y, decimalLatitude),
           decimalLongitude = ifelse(is.na(decimalLongitude) | decimalLongitude == "", decimalLongitude.y, decimalLongitude)) %>%
    dplyr::select(-c(decimalLatitude.y, decimalLongitude.y)) %>%
    mutate(coordsAdded = ifelse(is.na(decimalLatitude), F,T),
           cc_zero = ifelse(is.na(decimalLatitude), cc_zero, T),
           cc_equ = ifelse(is.na(decimalLatitude), cc_equ, T),
           clean_tests = ifelse(is.na(decimalLatitude), clean_tests, T),
           cc_err = ifelse(is.na(decimalLatitude), cc_err, T),
           georeferenced = ifelse(is.na(decimalLatitude), F, T))
  
  # bind matched and untouched records together
  all <- rbind(geo, nongeo.matched)# %>%
  # put occs back in original order
  reordered <- all[ order(match(all$uniqueCombinedID, occs$uniqueCombinedID)), ]
  
  return(reordered)
}


# apply function to all data 
matched.occs <- match_locality(flagged.occs)

# save data
write_csv(matched.occs, '02_cleaned_data/all_occs_cleaned.csv')


