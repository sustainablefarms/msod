# msod (Multi-Species Occupancy Detection Modelling)

__This package is in draft form: it is in the early stage of being designed__

Package for Bayesian fitting of a small number of occupancy detection models for multiple species.
Models are closely related to a JSDM by Tobler and models seen in the Boral package.

## Code Structure
+ Computations Common to All Models
   + linear combinations of environmental occupancy predictors for many parameter draws
   + linear combination of LV
   + apply a modelsite x theta --> likelihood function to get array of likelihood of sites and theta
   + combine occupancy probability with conditional detection probability to get full detection likilihoods (full probability distribution has too many dimensions)
   + computing a median, hpd, and expectation for occupancy from given occupancy probability per theta x LV x RE
   + computing a median, hpd, and expectation for expected species richness given occupancy probability per theta x LV x RE
   + processing of environmental data, and their saving to a fitted model

+ Computations that will differ between models
   + occupancy RV mean and standard deviation (it will always be normal though)
     + mean will be sum of external predictor component, LV component and RE component depending on the model
   + detection prob when occupied (but only for some models)
   + parameters to check and initialise for running JAGS
   + preparing random effects
   
+ Model Classes
  + inheritance structure: JSODM without LV, JSODM with LV > (JSODM with LV and RE, JSODM with spatially correlated LV)
  + Use S3 classes
    + generic functions:
      + likelihood
      + (predict) compute probabilities
      + simulate

## Functionality (draft)
### Core
+ fitting JSDMs with occupancy and detection components
   + latent variables
   + spatial correlated latent variables
   + random effects
   + probit tranfer function with linear predictor for occupancy
   + logit model for detection
   + some missing data
   
+ computing likelihood and log posterior density (lpd) for each ModelSite

+ building models with artifical parameters (core because very important for code testing)

+ predicting for multiple species and multiple sites

+ simulating from given models (core because very important for code testing)

### Small Helping Functionality
+ standardising and preparing data
+ ViF removal of highly correlated predictors
+ Rmd template for assessing model quality

### Non-Core Functionality (may not be available for all models)
+ Occupancy and Detection Residuals
+ Species richness predictions and comparison to observed

