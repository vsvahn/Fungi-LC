##############################################################################
# Build Bayesian logistic regression, make predictions 

# VJ Svahnstrom
##############################################################################

# load packages ----

library(brms)
library(tidyverse)
library(readr)
library(glue)
library(tidybayes)
library(pROC)

## load and prep data ----

# read in predictor data 
predictors <- read.csv('02_cleaned_data/predictors/all_sp_predictors_5.2.24.csv') %>% mutate(Threatened = as.factor(Threatened))

#split the species into ones that have been assessed (labelled)
# and those that haven’t (unlabelled) ----

labelled <- filter(predictors, ! is.na(Threatened)) 
unlabelled <- filter(predictors,  is.na(Threatened)) 


## fit the model ----

priors <- c(
  prior(student_t(3, -1, 1), class="Intercept")
)

fit <- brm(Threatened ~ yearPublished*logSpecimenCount + propContemporary,
            data = labelled, family = bernoulli())

## check fit ----
plot(fit)
summary(fit)

## model validation - k fold cross validation ----
kf = kfold(fit, save_fits = T)

## plot results for training species ----
labelled %>%
  mutate(Predicted = ifelse(Threatened == 0, 'Not threatened', 'Threatened')) %>%
  add_epred_draws(fit) %>%
  median_qi(.epred) %>% 
  arrange(.epred) %>% 
  mutate(index=1:nrow(.)) %>% 
  ggplot(aes(x=index, y=.epred, colour=Predicted)) +
  geom_point()+
  geom_errorbar(aes(ymin=.lower, ymax=.upper, x=index)) +
  scale_color_manual(values = c('#60C659', "#D81E05")) +
  theme_bw() +
  theme(axis.text.x = element_blank(), axis.ticks.x = element_blank()) +
  labs(x= '', y = 'Predicted probability threatened', color = 'Actual \nThreat status') +
  theme(legend.title.align=0.5) +  # Center align legend text
  geom_hline(yintercept = 0.5, linetype = "dashed") +
  geom_text(aes(x = 50, y = .6, label = "Predicted threatened"),
            color = "#D81E05", size = 4) + 
  geom_text(aes(x = 50, y = .4, label = "Predicted not threatened"),
            color = "#60C659", size = 4)

## make predictions ----
labelled  <- labelled %>%
  add_epred_draws(fit) %>%
  median_qi(.epred) 


## make predictions ----
predictions <-  unlabelled %>%
  add_epred_draws(fit) %>%
  median_qi(.epred) %>%
  mutate(pred.class = ifelse(.epred<0.5, 0, 1))


## calculate optimal threshold for confident LCs

# function to calculate threshold for a given false negative rate on assessed species
calculate_threshold_for_fnr <- function(observed, predicted, desired_fnr) {
  
  # Sort predictions in descending order
  sorted_predictions <- sort(predicted, decreasing = TRUE)
  
  # Calculate FNR for each threshold
  tp_fn <- sapply(sorted_predictions, function(threshold) {
    predicted_classes <- ifelse(predicted >= threshold, 1, 0)
    tp <- sum(predicted_classes == 1 & observed == 1)
    fn <- sum(predicted_classes == 0 & observed == 1)
    thresh <- threshold
    fns <- 1 - (tp / (tp + fn))
    c(tp, fn, thresh, fns)
  })
  
  # find optimal threshold for desired FNR
  optimal_threshold <- tp_fn[3,which.max(tp_fn[4,] <= desired_fnr)]
  
  return(optimal_threshold)
}

desired_fnr <- 0.05
threshold <- calculate_threshold_for_fnr(labelled$Threatened, labelled$.epred, 0.05)
cat("Threshold for FNR of", desired_fnr, ":", threshold, "\n")


## find confidently LC species --

# use threshold from above
# find species with 95% conf. int. of posterior dist. within threshold
num_lc <- sum(predictions$.upper < threshold)
glue("{(num_lc)} of macrofungi species confidently least concern")
lc <- predictions %>% filter(.upper < threshold)


