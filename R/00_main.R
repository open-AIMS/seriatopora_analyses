## ---- loadLibraries
library(tidyverse)
library(INLA)
library(glmmTMB)
library(brms)
library(posterior)
library(tidybayes)
library(DHARMa)
library(emmeans)
library(INLAutils)
library(patchwork)
library(standist)
library(inlatools)
library(bayesplot)
library(DT)
library(modelsummary)
config_modelsummary(factory_default = 'tinytable')
## ----end

## ---- loadFunctions
##setwd("R")
source("../R/helper_functions.R")
## ----end

## ---- preparePaths
source("../R/02_prepare_paths.R")
## ----end

## Read in the data provided by Mike (23/02/2024 via email)
## source("05_read_data.R")

## Process the data
source("10_process_data.R")
if (1 == 2) {
  ## Question 1 - temporal (quasi-spatial) patterns

  ## Prepare the data
  source("21_prepare_data.R")

  ## Exploratory data analysis
  source("22_eda.R")

  ## Fit the models
  source("23_fit_models.R")
}

## Question 2 - Before/After analyses
source("31_prepare_data.R")

## Exploratory data analysis

source("32_eda.R")

## Fit model

source("33_fit_models.R")

