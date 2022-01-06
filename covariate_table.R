
# Clear workspace of all objects and unload all extra (non-base) packages
rm(list = ls(all = TRUE))
if (!is.null(sessionInfo()$otherPkgs)) {
  res <- suppressWarnings(
    lapply(paste('package:', names(sessionInfo()$otherPkgs), sep=""),
           detach, character.only=TRUE, unload=TRUE, force=TRUE))
}

pacman::p_load(knitr, kableExtra,
               tidyverse
               )

source("0_Functions.R")

################################################################################################################
cov_train <- readRDS(file.path("..", "..", "ACT HEI Supp", "act_hei_aim1a", "Output", "mm_cov_train_set.rda"))

################################################################################################################

cov_train %>%
  select(log_m_to_a1:last_col()) %>%
  names() %>%
  as.data.frame() %>%
  rename(., covariate=.) %>%
  covariate_table_fn(dt = ., cov_name = "covariate") %>%
  kable(caption = "Geographic covariates (geocovariates) used in UK-PLS models") %>% 
  kable_styling()
