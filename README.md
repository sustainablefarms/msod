# msod (Multi-Species Occupancy Detection Modelling)

__This package is in draft form: it is in the early stage of being designed__

Package for Bayesian fitting of a small number of occupancy detection models for multiple species.
Models are closely related to a JSDM by Tobler and models seen in the Boral package.

## Code Structure
+ Computations Common to All Models
   + linear combinations of environmental occupancy predictors for many parameter draws
   + linear combination of LV
   + apply a modelsite x theta --> likelihood function to get array of likelihood of sites and theta
      + probably modelsite x theta x species arrays, with each value a probability of occupancy (detection), conditional on certain parameters, marginal on others.
   + combine occupancy probability with conditional detection probability to get full detection likilihoods (full probability distribution has too many dimensions) [doing this in multisite way - e.g. through large arrays]
   + computing a median, hpd, and expectation for occupancy from given occupancy probability per theta x LV x RE
   + computing a median, hpd, and expectation for expected species richness given occupancy probability per theta x LV x RE
   + processing of environmental data, and their saving to a fitted model
   + common arrays:
     + value: prob occupancy conditional on some, marginal on others dim: draw x modelsite x species x marginal param draw
     + value: prob detection conditional on occupancy, marginal on others. dim: draw x modelsite x species x marginal param draw
     + value: stdnormthresh for occupancy with marginal param sets. dim: draw x modelsite x species x marginal param draw

+ Computations that will differ between models
   + occupancy RV mean and standard deviation (it will always be normal though)
     + mean will be sum of external predictor component, LV component and RE component depending on the model
   + detection prob when occupied (but only for some models)
   + parameters to check and initialise for running JAGS
   + preparing random effects
   
+ Model Classes
  + inheritance structure: runjags > JSODM without LV > JSODM with LV > (JSODM with LV and RE, JSODM with spatially correlated LV). Actually this seems poor, I think starting with broadest and going more narrow (narrowest == JSODM witout LV) would be more appropriate.
  + Use S3 classes (I found this page useful https://elvinouyang.github.io/study%20notes/oop-with-r-s3-and-r6/)
    + generic functions:
      + logLik
      + (predict) compute probabilities
      + simulate
      + prepdata (missing values, LV and random effects parts)

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

+ remove the standardising of data to 1 with mean 0

### Small Helping Functionality
+ standardising and preparing data
+ ViF removal of highly correlated predictors
+ Rmd template for assessing model quality

### Non-Core Functionality (may not be available for all models)
+ Occupancy and Detection Residuals
+ Species richness predictions and comparison to observed

