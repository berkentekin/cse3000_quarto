library("tidyverse")
library("ggpmisc")
library("pcaPP")
library("VGAM")
library("purrr")
library("MonoPoly")
library("minpack.lm")

library(dplyr)
library(tidyr)
library(ggplot2)
library(MonoPoly)       # for monpol()
library(minpack.lm)     # for nlsLM if you want it
library(LaplacesDemon)  # for AICc()
library(cv)
library(EnvStats)
library(dplyr)
library(tidyr)
library(ggplot2)
library(MonoPoly)       # for monpol()
library(minpack.lm)     # for nlsLM if you want it
library(LaplacesDemon)  # for AICc()
library(cv)



snr_order <- c("n4", "n2", "0", "2", "4", "6", "8", "Q")

read_as_tibble <- function(filepath, sep = ",") {
  read.csv(filepath, sep = sep) %>% as_tibble()
}

noq <- function(data) {
  data %>% filter(snr != 'Q')
}

`%sp%` <- function(df, val) {
  filter(df, speaker == val)
}

`%snr%` <- function(df, val) {
  filter(df, snr == val)
}

tidier <- function(data) {
  data <- data %>% noq() %>% mutate(snr = factor(snr,snr_order)) %>% mutate(
    snr = as.numeric(gsub("^n", "-", as.character(snr))))
}

richards_start <- c(
  B = 0.04,
  M = 0.75,
  Q = 1
  
)

f_siib <- function(d, a, b) {
  (1 - exp(-1 * a * d)) ^ b
}

f_siib3 <- function(d, a, b, c) {
  ((1 - exp(-1 * a * d)) ^ b) * (1-c)
}

gomp <- function(d, b, c) {
  1-exp(-b * exp(-c * d))
}

gomp_3 <- function(d, a, b, c) {
  a - a * exp(-b * exp(-c * d))
}

AICc <- function(fit, n, k) {
  AIC(fit) +  (2*(k^2) + 2*k) / (n - k - 1)
}


logistic_1p1 <- function(d, a) {
  (1 - a) / (1 + exp(-d))
}

logistic_1p2 <- function(d, a) {
  1 / (1 + exp(a - d))
}

logistic_taal10 <- function(d, a, b) {
  1 / (1 + exp(a * d + b))
}

logistic_3_taal09 <- function(d, a, b, c) {
  1 / (1 + exp(a * log(d+c) + b))
}

logistic_4_taal09 <- function(d, a, b, c, h) {
  h / (1 + exp(a * log(d+c) + b))
}

dau_taal10 <- function(d, a, b, c) {
  1 / (1 + (a*d + b)^ c)
}


logistic_gen <- function(d, a, b, c) { #Three variables
  (1 - c) / (1 + exp(a * d + b))
}

logistic_gen_wrong <- function(d, a, b, c) { #Three variables
  c + ((1 - c) / (1 + exp(a * d + b)))
}

logistic_gen2 <- function(d, a, b, c, h) { #Four variables
  h + ( (1 - c - h) / (1 + exp(a * d + b)) )
}



prepare_allsstar <- function(metadata_csv_loc, ...) {
  # Helper function to read CSV files and convert them to tibbles
  
  
  # Gather the additional score file paths
  scores <- list(...)
  
  # Read the input files
  l1_data <- read_as_tibble(metadata_csv_loc, sep = ";")
  l1_stoi_scores <- read_as_tibble(scores[["stoi_scores_loc"]])
  l1_miknn_scores <- read_as_tibble(scores[["miknn_scores_loc"]])
  l1_siib_scores <- read_as_tibble(scores[["siib_scores_loc"]])
  l1_siib_scores_gaussian <- read_as_tibble(scores[["siib_scores_gaussian_loc"]])
  l1_siib_scores_gaussian_gapped <- read_as_tibble(scores[["siib_scores_gaussian_gapped_loc"]])
  
  # Merge all data
  allsstar_data <- l1_data %>%
    merge(l1_stoi_scores, by.x = c("audio", "snr"), by.y = c("filename", "snr"), all.x = TRUE) %>%
    merge(l1_miknn_scores, by.x = c("audio", "snr"), by.y = c("CleanFile", "Degradation"), all.x = TRUE) %>%
    merge(l1_siib_scores, by.x = c("audio", "snr"), by.y = c("filename", "snr"), all.x = TRUE) %>%
    merge(l1_siib_scores_gaussian, by.x = c("audio", "snr"), by.y = c("filename", "snr"), all.x = TRUE) %>%
    merge(l1_siib_scores_gaussian_gapped, by.x = c("audio", "snr"), by.y = c("filename", "snr"), all.x = TRUE)
  
  # Perform additional transformations
  allsstar_data <- allsstar_data %>% tidier() %>%
    rename(stoi_score = score) %>%
    #filter(trial > 3) %>%
    replace_na(list(stoi_score = 1)) %>%
    mutate(sim_wcr = autoscore / numwords) %>%
    mutate(parent_filename = sub("_S[0-9]+$", "", audio))
  
  # Return the final data
  return(allsstar_data)
}

prepare_dantale <- function(dantale_data, dantale_conditions, dantale_trials, dantale_oim) {
  
  
  dantale_oim <- dantale_oim %>% mutate(snr = as.numeric(snr), cond = as.numeric(cond), trial = as.numeric(trial))
  
  dantale_all <- dantale_data %>%
    merge(dantale_oim, by.x=c("experiment", "cond", "snr_no"), by.y=c("experiment", "cond", "snr"), all.x = TRUE, all.y = TRUE)
  
  dantale_all <- dantale_all %>% group_by(experiment,cond, snr) %>% summarise(
    experiment = first(experiment),
    across(
      everything(),
      ~ mean(.x, na.rm = TRUE),                  # The function to apply
      .names = "avg_{.col}"                      # How to name the resulting columns
    )
  ) %>% merge(
    dantale_conditions, by.x=c("experiment", "cond"), by.y=c("experiment", "cond")
  ) %>% mutate(snr = as.factor(snr)) %>%  filter(experiment != 'the_experiment_Training')
  
}

dantale_data <- read_csv("data/dantale_data.csv")
dantale_conditions <- read_csv("data/dantale_conditions.csv")
dantale_trials <- read_csv("data/dantale_trials.csv")
dantale_oim <- read_csv("data/dantale_results.csv")
dantale_other <- read_csv("data/dantale_results_other.csv")
dantale_miknn <- read_csv("data/dantale_results_miknn_2.csv")
dantale_oim <- dantale_oim %>% merge(dantale_other) %>% merge(dantale_miknn)
dantale_all <- prepare_dantale(dantale_data,dantale_conditions, dantale_trials, dantale_oim) 

allsstar_data <- prepare_allsstar(metadata_csv_loc = "data/l1res.csv",
                                  stoi_scores_loc="data/stoi_scores.csv",
                                  miknn_scores_loc="data/miknn_scores.csv",
                                  siib_scores_loc="data/siib_scores.csv",
                                  siib_scores_gaussian_loc="data/siib_scores_gaussian.csv",
                                  siib_scores_gaussian_gapped_loc="data/siib_scores_gaussian_gapped.csv")


data_group_by_snr_gender <- allsstar_data %>%
  group_by(snr) %>%
  summarise(
    avg_wcr = mean(sim_wcr),
    avg_stoi_score = mean(stoi_score),
    avg_siib_score_gaussian = mean(siib_score_gaussian_gapped),
    avg_miknn_score = mean(MIKNNScore2x),
        .groups = "drop"
  ) %>% select(where(~ !any(is.na(.)))) %>% mutate(snr = as.factor(snr))


### SIIB Patch


siib_scores <- read_csv("data/siib_scores_all.csv") %>% group_by(snr) %>% summarise(
  avg_siib_score = mean(siib_score),
  .groups = "drop"
) %>% select(where(~ !any(is.na(.)))) %>% mutate(snr = as.factor(snr))

stoi_scores_grouped <- read_csv("data/stoi_scores_grouped.csv") %>% group_by(snr) %>% summarise(
  avg_stoi_score_grouped = mean(stoi_score_grouped),
  .groups = "drop"
) %>% select(where(~ !any(is.na(.)))) %>% mutate(snr = as.factor(snr))

sim_wcr_grouped <- read_csv("data/sim_wcr_grouped.csv") %>% group_by(snr) %>% summarise(
  avg_sim_wcr_grouped = mean(sim_wcr_grouped),
  .groups = "drop"
) %>% select(where(~ !any(is.na(.)))) %>% mutate(snr = as.factor(snr))

data_group_by_snr_gender <- data_group_by_snr_gender %>% merge(siib_scores, by="snr") %>% merge(stoi_scores_grouped, by="snr")


###

dantale_conds <- dantale_all %>% group_by(cond, snr) %>% summarise(across(everything(), ~ mean(.x, na.rm = TRUE)), .groups = "drop") %>%
  select(where(~ !any(is.na(.)))) %>% mutate(snr = as.factor(snr))

dat_all <- bind_rows(dantale_conds, data_group_by_snr_gender)

dat_allsstar <- data_group_by_snr_gender


model_specs <- list(
  #linear = list(
  #  method       = "lm",
  #  formula_expr = quote(y ~ x)  # 'quote()' so we can substitute later
  #),
  f_siib = list(
    name = "f_siib",
    method       = "nlsLM",
    formula_expr = quote(y ~ f_siib(x, a, b)),
    start        = c(a = 0.01, b = 1)
  ),
  f_siib3 = list(
    name = "f_siib3",
    method       = "nlsLM",
    formula_expr = quote(y ~ f_siib3(x, a, b, c)),
    start        = c(a = 0.01, b = 1, c = 0.95)
  ),
  monpoly3 = list(
    name = "monpoly3",
    method       = "monpol",
    formula_expr = quote(y ~ x),
    degree       = 3
  ),
  logistic = list(
    name = "logistic_2",
    method       = "nlsLM",
    formula_expr = quote(y ~ logistic_taal10(x, a, b)),
    start        = c(a = 0, b = 0.5)
  ),
  logistic_gen = list(
    name = "logistic_3",
    method       = "nlsLM",
    formula_expr = quote(y ~ logistic_gen(x, a, b, c)),
    start        = c(a = 0, b = 0.5, c=0),
    lower = c(a = -Inf, b = -Inf, c = 0),
    upper = c(a = Inf, b = Inf, c = 0.05)
  ),
  logistic_gen2 = list(
    name = "logistic_4",
    method       = "nlsLM",
    formula_expr = quote(y ~ logistic_gen2(x, a, b, c, h)),
    start        = c(a = 0, b = 0.5, c=0, h = 0),
    lower = c(a = -Inf, b = -Inf, c = 0, h = 0),
    upper = c(a = Inf, b = Inf, c = 0.05, h = 0.05)
  ),
  dau = list(
    name = "dau",
    method       = "nlsLM",
    formula_expr = quote(y ~ dau_taal10(x, a, b, c)),
    start        = c(a = 0, b = 1, c = 0.5)
  ), 
  gompertz = list(
    name = "gompertz",
    method       = "nlsLM",
    formula_expr = quote(y ~ gomp(x, b, c)),
    start        = c(b=1, c=-0.02)
  ),
  gompertz_3 = list(
    name = "gompertz_3",
    method       = "nlsLM",
    formula_expr = quote(y ~ gomp_3(x,a, b, c)),
    start        = c(a = 1, b=1, c=-0.02)
  ),
  gompertz_4 = list(
    name = "gompertz_4",
    method       = "nlsLM",
    formula_expr = quote(y ~ gomp_4(x, a, b, c, h)),
    start        = c(a = 1, b = 1, c = -0.01, h = 0)
  )
  
)



model_specs_10 <- list(
  linear = list(
    name = "linear",
     method       = "lm",
    formula_expr = quote(y ~ x)  # 'quote()' so we can substitute later
  ),

  monpoly3 = list(
    name = "monpoly3",
    method       = "monpol",
    formula_expr = quote(y ~ x),
    degree       = 3
  ),
  logistic = list(
    name = "logistic_2",
    method       = "nlsLM",
    formula_expr = quote(y ~ logistic_taal10(x, a, b)),
    start        = c(a = 0, b = 0.5)
  ),
  logistic_gen = list(
    name = "logistic_3",
    method       = "nlsLM",
    formula_expr = quote(y ~ logistic_gen(x, a, b, c)),
    start        = c(a = 0, b = 0.5, c=0),
    lower = c(a = -Inf, b = -Inf, c = 0),
    upper = c(a = Inf, b = Inf, c = 0.05)
  ),
  logistic_gen2 = list(
    name = "logistic_4",
    method       = "nlsLM",
    formula_expr = quote(y ~ logistic_gen2(x, a, b, c, h)),
    start        = c(a = 0, b = 0.5, c=0, h = 0),
    lower = c(a = -Inf, b = -Inf, c = 0, h = 0),
    upper = c(a = Inf, b = Inf, c = 0.05, h = 0.05)
  ),
  dau = list(
    name = "dau",
    method       = "nlsLM",
    formula_expr = quote(y ~ dau_taal10(x, a, b, c)),
    start        = c(a = 0, b = 1, c = 0.5)
  ), 
  gompertz = list(
    name = "gompertz",
    method       = "nlsLM",
    formula_expr = quote(y ~ gomp(x, b, c)),
    start        = c(b=1, c=-0.02)
  ),
  gompertz_3 = list(
    name = "gompertz_3",
    method       = "nlsLM",
    formula_expr = quote(y ~ gomp_3(x,a, b, c)),
    start        = c(a = 1, b=-0.2, c=-0.24)
  )
  #gompertz_4 = list(
  #  name = "gompertz_4",
  #  method       = "nlsLM",
  #  formula_expr = quote(y ~ gomp_4(x, a, b, c, h)),
  #  start        = c(a = 1, b = 1, c = -0.01, h = 0)
  #)
  
)

dat_allsstar_snr_wcr_res <- fit_many_models(
  data       = dat_n,
  x_col      = "snr",
  y_col      = "avg_wcr",
  model_list = model_specs_10,
  n          = nrow(dat_allsstar)
)
dat_allsstar_snr_wcr_res$summary
dat_allsstar_snr_wcr_res$coefs
ggplot(dat_allsstar_snr_wcr_res$predictions, aes(x = snr, y = avg_wcr)) +
  geom_point(data = dat_allsstar, aes(x = snr, y = avg_wcr)) +
  geom_line(aes(y = FittedValue)) +
  facet_wrap(~Model, scales = "free_y")


dat_allsstar_stoi_wcr_res <- fit_many_models(
  data       = dat_n,
  x_col      = "avg_stoi_score",
  y_col      = "avg_wcr",
  model_list = model_specs_10,
  n          = nrow(dat_allsstar)
)

dat_allsstar_stoi_wcr_res$summary
dat_allsstar_stoi_wcr_res$coefs
ggplot(dat_allsstar_stoi_wcr_res$predictions, aes(x = avg_stoi_score, y = avg_wcr)) +
  geom_point(data = dat_allsstar, aes(x = avg_stoi_score, y = avg_wcr)) +
  geom_line(aes(y = FittedValue)) +
  facet_wrap(~Model, scales = "free_y")

model_specs_1_stoi_log <- list(
  logistic = list(
    name = "logistic",
    method       = "nlsLM",
    formula_expr = quote(y ~ logistic_taal10(x, a, b)),
    start        = c(a = 0, b = -0.5)
  ),
  logistic_3_taal09 = list(
    name = "logistic_3_taal09",
    method       = "nlsLM",
    formula_expr = quote(y ~ logistic_3_taal09(x, a, b, c)),
    start        = c(a = 1, b = -0.5, c=0.1),
    lower = c(a = -Inf, b = -Inf, c = 0),
    upper = c(a = Inf, b = Inf, c = 0.05)
  ),
  logistic_3_free = list(
    name = "logistic_3_free",
    method       = "nlsLM",
    formula_expr = quote(y ~ logistic_gen(x, a, b, c)),
    start        = c(a = -10.878974, b = 6.210622, c=0),
    lower = c(a = -Inf, b = -Inf, c = -Inf),
    upper = c(a = Inf, b = Inf, c = Inf)
  ),
  logistic_4_taal09 = list(
    name = "logistic_4_taal09",
    method       = "nlsLM",
    formula_expr = quote(y ~ logistic_4_taal09(x, a, b, c, h)),
    start        = c(a = 1, b = -0.5, c=0.05, h = 0.95),
    lower = c(a = -Inf, b = -Inf, c = 0, h = 0.90),
    upper = c(a = Inf, b = Inf, c = 0.01, h = 1)
  ),
  logistic_4_free = list(
      name = "logistic_4",
      method       = "nlsLM",
      formula_expr = quote(y ~ logistic_gen2(x, a, b, c, h)),
      start        = c(a = -10.878974, b = 6.210622, c=0, h = 0),
      lower = c(a = -Inf, b = -Inf, c = -Inf, h = -Inf),
      upper = c(a = Inf, b = Inf, c = Inf, h = Inf)
  )
)

dat_allsstar_stoi_res <- fit_many_models(
  data       = dat_allsstar,
  x_col      = "avg_stoi_score",
  y_col      = "avg_wcr",
  model_list = model_specs_1_stoi_log,
  n          = nrow(dat_allsstar)
)


dat_allsstar_stoi_res$summary
dat_allsstar_stoi_res$coefs
dat_allsstar_stoi_res$predictions
p <- ggplot(
  dat_allsstar_stoi_res$predictions, 
  aes(x = avg_stoi_score, y = FittedValue, color = Model)
) +
  # Thicker lines
  geom_line(size = 1.5) +
  # Points from the original data, without inheriting color/model
  geom_point(
    data = dat_allsstar, 
    aes(x = avg_stoi_score, y = avg_wcr), 
    shape = 1,
    inherit.aes = FALSE
  )+
  # Change the legend title and individual label names:
  scale_color_discrete(
    name = "Model Type",
    labels = c(
      "logistic_3_free_pred"    = "Logistic, 3 parameters",
      "logistic_3_taal09_pred"  = "Logistic, 3 parameters (Taal '09)",
      "logistic_4_free_pred"    = "Logistic, 4 parameters",
      "logistic_4_taal09_pred"  = "Logistic, 4 parameters (Taal '09)",
      "logistic_pred"           = "Logistic, 2 parameters (Taal)"
    )
  ) +
  # Force y to end at 1
  scale_y_continuous(labels = function(y) y * 100, limits = c(NA, 1)) +
  labs(x = "Average STOI Score",
         y = "Average Word Correct Ratio (%)")
  theme_minimal()
# Print the plot in R to check it
print(p)



model_specs_1_miknn_log <- list(
  logistic = list(
    name = "logistic",
    method       = "nlsLM",
    formula_expr = quote(y ~ logistic_taal10(x, a, b)),
    start        = c(a = 0, b = -0.5)
  ),
  logistic_3_taal09 = list(
    name = "logistic_3_taal09",
    method       = "nlsLM",
    formula_expr = quote(y ~ logistic_3_taal09(x, a, b, c)),
    start        = c(a = 0, b = -0.5, c=0.1),
    lower = c(a = -Inf, b = -Inf, c = 0),
    upper = c(a = Inf, b = Inf, c = 0.05)
  ),
  logistic_3_free = list(
    name = "logistic_3_wrong",
    method       = "nlsLM",
    formula_expr = quote(y ~ logistic_gen(x, a, b, c)),
    start        = c(a = -0.232, b = 8.992, c=-0.1),
    lower = c(a = -Inf, b = -Inf, c = -Inf),
    upper = c(a = Inf, b = Inf, c = Inf)
  ),
  logistic_4_taal09 = list(
    name = "logistic_4_taal09",
    method       = "nlsLM",
    formula_expr = quote(y ~ logistic_4_taal09(x, a, b, c, h)),
    start        = c(a = 1, b = -0.5, c=0.1, h = 1),
    lower = c(a = -Inf, b = -Inf, c = 0, h = 0.90),
    upper = c(a = Inf, b = Inf, c = 0.05, h = 1)
  ),
  logistic_4_free = list(
    name = "logistic_4",
    method       = "nlsLM",
    formula_expr = quote(y ~ logistic_gen2(x, a, b, c, h)),
    start        = c(a = -0.232, b = 8.992, c=-0.1, h = -0.1),
    lower = c(a = -Inf, b = -Inf, c = -Inf, h = -Inf),
    upper = c(a = Inf, b = Inf, c = Inf, h = Inf)
  )
)


dat_allsstar_miknn_res <- fit_many_models(
  data       = dat_allsstar,
  x_col      = "avg_miknn_score",
  y_col      = "avg_wcr",
  model_list = model_specs_1_miknn_log,
  n          = nrow(dat_allsstar)
)

dat_allsstar_miknn_res$summary
dat_allsstar_miknn_res$coefs

p <- ggplot(
  dat_allsstar_miknn_res$predictions, 
  aes(x = avg_miknn_score, y = FittedValue, color = Model)
) +
  # Thicker lines
  geom_line(size = 1.5) +
  # Points from the original data, without inheriting color/model
  geom_point(
    data = dat_allsstar, 
    aes(x = avg_miknn_score, y = avg_wcr), 
    shape = 1,
    inherit.aes = FALSE
  )+
  # Change the legend title and individual label names:
  scale_color_discrete(
    name = "Model Type",
    labels = c(
      "logistic_3_free_pred"    = "Logistic, 3 parameters",
      "logistic_3_taal09_pred"  = "Logistic, 3 parameters (Taal '09)",
      "logistic_4_free_pred"    = "Logistic, 4 parameters",
      "logistic_4_taal09_pred"  = "Logistic, 4 parameters (Taal '09)",
      "logistic_pred"           = "Logistic, 2 parameters (Taal)"
    )
  ) +
  # Force y to end at 1
  scale_y_continuous(labels = function(y) y * 100, limits = c(NA, 1)) +
  labs(x = "Average MIKNN score",
       y = "Average Word Correct Ratio")
theme_minimal()
# Print the plot in R to check it
print(p)

ggplot(dat_allsstar_miknn_res$predictions, aes(x = avg_miknn_score, y = avg_wcr)) +
  geom_point(data = dat_allsstar, aes(x = avg_miknn_score, y = avg_wcr)) +
  geom_line(aes(y = FittedValue)) +
  facet_wrap(~Model, scales = "free_y")

model_specs_1_siib <- list(
  f_siib = list(
    name = "f_siib3",
    method       = "nlsLM",
    formula_expr = quote(y ~ f_siib(x, a, b)),
    start        = c(a = 0.01, b = 1)
  ),
  f_siib3 = list(
    name = "f_siib3",
    method       = "nlsLM",
    formula_expr = quote(y ~ f_siib3(x, a, b, c)),
    start        = c(a = 0.01, b = 1, c = 0.05)
  ),
  logistic = list(
    name = "logistic",
    method       = "nlsLM",
    formula_expr = quote(y ~ logistic_taal10(x, a, b)),
    start        = c(a = 0, b = -0.5)
  ),
  logistic_3 = list(
    name = "logistic_3",
    method       = "nlsLM",
    formula_expr = quote(y ~ logistic_gen(x, a, b, c)),
    start        = c(a = -0.0196, b = 1.0052 , c=0) # Actual 0.05724

  ),
  logistic_3_taal09 = list(
    name = "logistic_3_taal09",
    method       = "nlsLM",
    formula_expr = quote(y ~ logistic_3_taal09(x, a, b, c)),
    start        = c(a = 0, b = -0.5, c=0.1),
    lower = c(a = -Inf, b = -Inf, c = 0),
    upper = c(a = Inf, b = Inf, c = 0.05)
  )
)

dat_allsstar_siib_res <- fit_many_models(
  data       = dat_allsstar,
  x_col      = "avg_siib_score.x",
  y_col      = "avg_wcr",
  model_list = model_specs_1_siib,
  n          = nrow(dat_allsstar)
)

dat_allsstar_siib_res$summary
dat_allsstar_siib_res$coefs
ggplot(dat_allsstar_siib_res$predictions, aes(x = avg_siib_score, y = avg_wcr)) +
  geom_point(data = dat_allsstar, aes(x = avg_siib_score, y = avg_wcr)) +
  geom_line(aes(y = FittedValue)) +
  facet_wrap(~Model, scales = "free_y")

p <- ggplot(
  dat_allsstar_siib_res$predictions, 
  aes(x = avg_siib_score.x, y = FittedValue, color = Model)
) +
  # Thicker lines
  geom_line(size = 1.5) +
  # Points from the original data, without inheriting color/model
  geom_point(
    data = dat_allsstar, 
    aes(x = avg_siib_score.x, y = avg_wcr), 
    shape = 1,
    inherit.aes = FALSE
  )+
  # Change the legend title and individual label names:
  scale_color_discrete(
    name = "Model Type",
    labels = c(
      "f_siib_pred"    = "SIIB, 2 parameters",
      "f_siib3_pred"  = "SIIB, 3 parameters",
      "logistic_3_pred"    = "Logistic, 3 parameters",
      "logistic_3_taal09_pred"  = "Logistic, 3 parameters (Taal '09)",
      "logistic_pred"           = "Logistic, 2 parameters (Taal)"
    ),
    limits = c("f_siib_pred", "f_siib3_pred", "logistic_pred", "logistic_3_pred", "logistic_3_taal09_pred")
  ) +
  # Force y to end at 1
  scale_y_continuous(labels = function(y) y * 100, limits = c(NA, 1)) +
  labs(x = "Average SIIB score",
       y = "Average Word Correct Ratio (%)")
# Print the plot in R to check it
print(p)

model_specs_1_siib_log <- list(
  logistic = list(
    name = "logistic",
    method       = "nlsLM",
    formula_expr = quote(y ~ logistic_taal10(x, a, b)),
    start        = c(a = 0, b = -0.5)
  ),
  logistic_3 = list(
    name = "logistic_3",
    method       = "nlsLM",
    formula_expr = quote(y ~ logistic_gen(x, a, b, c)),
    start        = c(a = -0.0196, b = 1.0052 , c=0), # Actual 0.05724
    lower = c(a = -Inf, b = -Inf, c = 0),
    upper = c(a = Inf, b = Inf, c = 0.1)
  ),
  logistic_3_taal09 = list(
    name = "logistic_3_taal09",
    method       = "nlsLM",
    formula_expr = quote(y ~ logistic_3_taal09(x, a, b, c)),
    start        = c(a = 0, b = -0.5, c=0.1),
    lower = c(a = -Inf, b = -Inf, c = 0),
    upper = c(a = Inf, b = Inf, c = 0.05)
  )
  #logistic_4 = list(
  #  name = "logistic_4",
  #  method       = "nls",
  #  formula_expr = quote(y ~ logistic_gen2(x, a, b, c, h)),
  #  start        = c(a = -0.0196, b = 1.0052, c=0, h = 0),
  #  lower = c(a = -Inf, b = -Inf, c = 0, h = 0),
  #  upper = c(a = Inf, b = Inf, c = Inf, h = Inf)
  #),
  #logistic_4_free = list(
  #  name = "logistic_4",
  #  method       = "nls",
  #  formula_expr = quote(y ~ logistic_gen2(x, a, b, c, h)),
  #  start        = c(a = -0.0196, b = 1.0052, c=0, h = 0),
  #  lower = c(a = -Inf, b = -Inf, c = 0, h = 0),
  #  upper = c(a = Inf, b = Inf, c = Inf, h = Inf)
  #)
)


dat_allsstar_siib_res <- fit_many_models(
  data       = dat_allsstar,
  x_col      = "avg_siib_score",
  y_col      = "avg_wcr",
  model_list = model_specs_1_siib_log,
  n          = nrow(dat_allsstar)
)

dat_allsstar_siib_res$summary
dat_allsstar_siib_res$coefs
ggplot(dat_allsstar_siib_res$predictions, aes(x = avg_siib_score, y = avg_wcr)) +
  geom_point(data = dat_allsstar, aes(x = avg_siib_score, y = avg_wcr)) +
  geom_line(aes(y = FittedValue)) +
  facet_wrap(~Model, scales = "free_y")


model_specs_1_siib_gauss_log <- list(
  f_siib = list(
    name = "f_siib",
    method       = "nlsLM",
    formula_expr = quote(y ~ f_siib(x, a, b)),
    start        = c(a = 0.01, b = 1)
  ),
  f_siib3 = list(
    name = "f_siib3",
    method       = "nlsLM",
    formula_expr = quote(y ~ f_siib3(x, a, b, c)),
    start        = c(a = 0.01, b = 1, c = 0)
  ),
  logistic = list(
    name = "logistic",
    method       = "nlsLM",
    formula_expr = quote(y ~ logistic_taal10(x, a, b)),
    start        = c(a = 0, b = -0.5)
  ),
  logistic_3 = list(
    name = "logistic_3",
    method       = "nlsLM",
    formula_expr = quote(y ~ logistic_gen(x, a, b, c)),
    start        = c(a = -0.0722, b = 1.5338 , c=0), # Actual 0.05724
    lower = c(a = -Inf, b = -Inf, c = 0),
    upper = c(a = Inf, b = Inf, c = 1)
  ),
  logistic_3_taal09 = list(
    name = "logistic_3_taal09",
    method       = "nlsLM",
    formula_expr = quote(y ~ logistic_3_taal09(x, a, b, c)),
    start        = c(a = 0, b = -0.5, c=0.1),
    lower = c(a = -Inf, b = -Inf, c = 0),
    upper = c(a = Inf, b = Inf, c = 0.05)
  )
)

dat_allsstar_siib_gauss_log_res <- fit_many_models(
  data       = dat_allsstar,
  x_col      = "avg_siib_score_gaussian",
  y_col      = "avg_wcr",
  model_list = model_specs_1_siib_gauss_log,
  n          = nrow(dat_allsstar)
)

dat_allsstar_siib_gauss_log_res$summary
dat_allsstar_siib_gauss_log_res$coefs
ggplot(dat_allsstar_siib_gauss_log_res$predictions, aes(x = avg_siib_score_gaussian, y = avg_wcr)) +
  geom_point(data = dat_allsstar, aes(x = avg_siib_score_gaussian, y = avg_wcr)) +
  geom_line(aes(y = FittedValue)) +
  facet_wrap(~Model, scales = "free_y")


p <- ggplot(
  dat_allsstar_siib_gauss_log_res$predictions, 
  aes(x = avg_siib_score_gaussian, y = FittedValue, color = Model)
) +
  # Thicker lines
  geom_line(size = 1.5) +
  # Points from the original data, without inheriting color/model
  geom_point(
    data = dat_allsstar, 
    aes(x = avg_siib_score_gaussian, y = avg_wcr), 
    shape = 1,
    inherit.aes = FALSE
  )+
  # Change the legend title and individual label names:
  scale_color_discrete(
    name = "Model Type",
    labels = c(
      "f_siib_pred"    = "SIIB, 2 parameters",
      "f_siib3_pred"  = "SIIB, 3 parameters",
      "logistic_3_pred"    = "Logistic, 3 parameters",
      "logistic_3_taal09_pred"  = "Logistic, 3 parameters (Taal '09)",
      "logistic_pred"           = "Logistic, 2 parameters (Taal)"
    ),
    limits = c("f_siib_pred", "f_siib3_pred", "logistic_pred", "logistic_3_pred", "logistic_3_taal09_pred")
    
  ) +
  # Force y to end at 1
  scale_y_continuous(labels = function(y) y * 100, limits = c(NA, 1)) +
  labs(x = "Average SIIB_gauss score",
       y = "Average Word Correct Ratio (%)")
print(p)





 model_specs_siib <- list(
  #linear = list(
  #  method       = "lm",
  #  formula_expr = quote(y ~ x)  # 'quote()' so we can substitute later
  #),
  f_siib = list(
    name = "f_siib",
    method       = "nlsLM",
    formula_expr = quote(y ~ f_siib(x, a, b)),
    start        = c(a = 0.01, b = 1)
  ),
  f_siib3 = list(
    name = "f_siib3",
    method       = "nlsLM",
    formula_expr = quote(y ~ f_siib3(x, a, b, c)),
    start        = c(a = 0.01, b = 1, c = 1)
  ),
  logistic = list(
    name = "logistic",
    method       = "nlsLM",
    formula_expr = quote(y ~ logistic_taal10(x, a, b)),
    start        = c(a = 0, b = -0.5)
  ),
  logistic_gen = list(
    name = "logistic_gen",
    method       = "nlsLM",
    formula_expr = quote(y ~ logistic_gen(x, a, b, c)),
    start        = c(a = 0, b = -0.5, c=0),
    lower = c(a = -Inf, b = -Inf, c = 0),
    upper = c(a = Inf, b = Inf, c = 0.05)
  ),
  logistic_gen2 = list(
    name = "logistic_gen2",
    method       = "nlsLM",
    formula_expr = quote(y ~ logistic_gen2(x, a, b, c, h)),
    start        = c(a = 0, b = -0.5, c=0, h = 0),
    lower = c(a = -Inf, b = -Inf, c = 0),
    upper = c(a = Inf, b = Inf, c = 0.05)
  ),
  dau = list(
    name = "dau",
    method       = "nlsLM",
    formula_expr = quote(y ~ dau_taal10(x, a, b, c)),
    start        = c(a = 0, b = 1, c = 0.5)
  ), 
  #gompertz = list(
  #  name = "gompertz",
  #  method       = "nlsLM",
  #  formula_expr = quote(y ~ gomp(x, b, c)),
  #  start        = c(b=1, c=-0.02)
  #),
  gompertz_4 = list(
    name = "gompertz_4",
    method       = "nlsLM",
    formula_expr = quote(y ~ gomp_4(x, a, b, c, d0)),
    start        = c(a = 1, b = 1, c = -0.01, d0 = 0)
  )
  
)


model_specs_1 <- list(
  #linear = list(
  #  method       = "lm",
  #  formula_expr = quote(y ~ x)  # 'quote()' so we can substitute later
  #),
  #monpoly3 = list(
  #  name = "monpoly3",
  #  method       = "monpol",
  #  formula_expr = quote(y ~ x),
  #  degree       = 3
  #),
  logistic = list(
    name = "logistic",
    method       = "nlsLM",
    formula_expr = quote(y ~ logistic_taal10(x, a, b)),
    start        = c(a = 0, b = -0.5)
  ),
  logistic_3 = list(
    name = "logistic_gen",
    method       = "nlsLM",
    formula_expr = quote(y ~ logistic_gen(x, a, b, c)),
    start        = c(a = 0, b = -0.5, c=0),
    lower = c(a = -Inf, b = -Inf, c = 0),
    upper = c(a = Inf, b = Inf, c = 0.05)
  ),
  logistic_3_wrong = list(
    name = "logistic_gen_wrong",
    method       = "nlsLM",
    formula_expr = quote(y ~ logistic_gen_wrong(x, a, b, c)),
    start        = c(a = 0, b = -0.5, c=0),
    lower = c(a = -Inf, b = -Inf, c = 0),
    upper = c(a = Inf, b = Inf, c = 0.05)
  ),
  logistic_4 = list(
    name = "logistic_gen2",
    method       = "nlsLM",
    formula_expr = quote(y ~ logistic_gen2(x, a, b, c, h)),
    start        = c(a = 0, b = -0.5, c=0, h = 0),
    lower = c(a = -Inf, b = -Inf, c = 0),
    upper = c(a = Inf, b = Inf, c = 0.05)
  )
  #dau = list(
  #  name = "dau",
  #  method       = "nlsLM",
  #  formula_expr = quote(y ~ dau_taal10(x, a, b, c)),
  #  start        = c(a = 0, b = 1, c = 0.5)
  #), 
  #gompertz = list(
  #  name = "gompertz",
  #  method       = "nlsLM",
  #  formula_expr = quote(y ~ gomp(x, b, c)),
  #  start        = c(b=1, c=-0.02)
  #),
  #gompertz_4 = list(
  #  name = "gompertz_4",
  #  method       = "nlsLM",
  #  formula_expr = quote(y ~ gomp_4(x, a, b, c, d0)),
  #  start        = c(a = 1, b = 1, c = -0.01, d0 = 0)
  #)
  
)



model_specs_1_siib_gaus_log <- list(
  #linear = list(
  #  method       = "lm",
  #  formula_expr = quote(y ~ x)  # 'quote()' so we can substitute later
  #),
  #monpoly3 = list(
  #  name = "monpoly3",
  #  method       = "monpol",
  #  formula_expr = quote(y ~ x),
  #  degree       = 3
  #),
  logistic = list(
    name = "logistic",
    method       = "nlsLM",
    formula_expr = quote(y ~ logistic_taal10(x, a, b)),
    start        = c(a = 0, b = -0.5)
  ),
  logistic_3 = list(
    name = "logistic_gen",
    method       = "nlsLM",
    formula_expr = quote(y ~ logistic_gen(x, a, b, c)),
    start        = c(a = -0.05, b = -0.5, c=0),
    lower = c(a = -Inf, b = -Inf, c = 0),
    upper = c(a = Inf, b = Inf, c = 0.05)
  ),
  logistic_3_wrong = list(
    name = "logistic_gen_wrong",
    method       = "nlsLM",
    formula_expr = quote(y ~ logistic_gen_wrong(x, a, b, c)),
    start        = c(a = -0.05, b = -0.5, c=0),
    lower = c(a = -Inf, b = -Inf, c = 0),
    upper = c(a = Inf, b = Inf, c = 0.05)
  ),
  logistic_4 = list(
    name = "logistic_gen2",
    method       = "nlsLM",
    formula_expr = quote(y ~ logistic_gen2(x, a, b, c, h)),
    start        = c(a = 0, b = -0.5, c=0, h = 0),
    lower = c(a = -Inf, b = -Inf, c = 0),
    upper = c(a = Inf, b = Inf, c = 0.05)
  )
  #dau = list(
  #  name = "dau",
  #  method       = "nlsLM",
  #  formula_expr = quote(y ~ dau_taal10(x, a, b, c)),
  #  start        = c(a = 0, b = 1, c = 0.5)
  #), 
  #gompertz = list(
  #  name = "gompertz",
  #  method       = "nlsLM",
  #  formula_expr = quote(y ~ gomp(x, b, c)),
  #  start        = c(b=1, c=-0.02)
  #),
  #gompertz_4 = list(
  #  name = "gompertz_4",
  #  method       = "nlsLM",
  #  formula_expr = quote(y ~ gomp_4(x, a, b, c, d0)),
  #  start        = c(a = 1, b = 1, c = -0.01, d0 = 0)
  #)
  
)



model_specs_2_stoi <- list(
  #linear = list(
  #  method       = "lm",
  #  formula_expr = quote(y ~ x)  # 'quote()' so we can substitute later
  #),
  f_siib = list(
    name = "f_siib",
    method       = "nlsLM",
    formula_expr = quote(y ~ f_siib(x, a, b)),
    start        = c(a = 0.01, b = 1)
  ),
  f_siib3 = list(
    name = "f_siib3",
    method       = "nlsLM",
    formula_expr = quote(y ~ f_siib3(x, a, b, c)),
    start        = c(a = 9.490175, b = 164.232140, c = 1)
  ),
  monpoly3 = list(
    name = "monpoly3",
    method       = "monpol",
    formula_expr = quote(y ~ x),
    degree       = 3
  ),
  logistic = list(
    name = "logistic",
    method       = "nlsLM",
    formula_expr = quote(y ~ logistic_taal10(x, a, b)),
    start        = c(a = 0, b = -0.5)
  ),
  logistic_3 = list(
    name = "logistic_3",
    method       = "nlsLM",
    formula_expr = quote(y ~ logistic_gen(x, a, b, c)),
    start        = c(a = 0, b = -0.5 , c=0), # Actual 0.05724
    lower = c(a = -Inf, b = -Inf, c = 0),
    upper = c(a = Inf, b = Inf, c = 1)
  ),
  dau = list(
    name = "dau",
    method       = "nlsLM",
    formula_expr = quote(y ~ dau_taal10(x, a, b, c)),
    start        = c(a = 2, b = 0, c = -3)
  ), 
  gompertz = list(
    name = "gompertz",
    method       = "nlsLM",
    formula_expr = quote(y ~ gomp(x, b, c)),
    start        = c(b=0.1, c=-3.7)
  ),
  gompertz_3 = list(
    name = "gompertz_3",
    method       = "nlsLM",
    formula_expr = quote(y ~ gomp_3(x,a, b, c)),
    start        = c(a = 1, b=0.1, c=--3.7)
  )
  #gompertz_4 = list(
  #  name = "gompertz_4",
  #  method       = "nls",
  #  formula_expr = quote(y ~ gomp_4(x, a, b, c, h)),
  #  start        = c(a = 1, b =10.98, c = -1.44, h = 0)
  #)
  
)

model_specs_2_stoi_res <- fit_many_models(
  data       = dat_allsstar,
  x_col      = "avg_stoi_score",
  y_col      = "avg_wcr",
  model_list = model_specs_2_stoi,
  n          = nrow(dat_allsstar)
)

model_specs_2_stoi_res$summary
model_specs_2_stoi_res$coefs
ggplot(model_specs_2_stoi_res$predictions, aes(x = avg_stoi_score, y = avg_wcr)) +
  geom_point(data = dat_allsstar, aes(x = avg_stoi_score, y = avg_wcr)) +
  geom_line(aes(y = FittedValue)) +
  facet_wrap(~Model, scales = "free_y")

model_specs_2_miknn <- list(
  #linear = list(
  #  method       = "lm",
  #  formula_expr = quote(y ~ x)  # 'quote()' so we can substitute later
  #),
  f_siib = list(
    name = "f_siib",
    method       = "nlsLM",
    formula_expr = quote(y ~ f_siib(x, a, b)),
    start        = c(a = 0.01, b = 1)
  ),
  f_siib3 = list(
    name = "f_siib3",
    method       = "nlsLM",
    formula_expr = quote(y ~ f_siib3(x, a, b, c)),
    start        = c(a = 0.202, b = 1827.545, c = 0.95)
  ),
  monpoly3 = list(
    name = "monpoly3",
    method       = "monpol",
    formula_expr = quote(y ~ x),
    degree       = 3
  ),
  logistic = list(
    name = "logistic",
    method       = "nlsLM",
    formula_expr = quote(y ~ logistic_taal10(x, a, b)),
    start        = c(a = 0, b = -0.5)
  ),
  logistic_3 = list(
    name = "logistic_3",
    method       = "nlsLM",
    formula_expr = quote(y ~ logistic_gen(x, a, b, c)),
    start        = c(a = -0.0722, b = 1.5338 , c=0), # Actual 0.05724
    lower = c(a = -Inf, b = -Inf, c = 0),
    upper = c(a = Inf, b = Inf, c = 0.1)
  ),
  dau = list(
    name = "dau",
    method       = "nlsLM",
    formula_expr = quote(y ~ dau_taal10(x, a, b, c)),
    start        = c(a = 0.03, b = -1, c = -0.5)
  ), 
  gompertz = list(
    name = "gompertz",
    method       = "nlsLM",
    formula_expr = quote(y ~ gomp(x, b, c)),
    start        = c(b=0.05, c=-0.1)
  ),
  gompertz_3 = list(
    name = "gompertz_3",
    method       = "nlsLM",
    formula_expr = quote(y ~ gomp_3(x,a, b, c)),
    start        = c(a = 1, b=0.00909, c=--0.05)
  )
  #gompertz_4 = list(
  #  name = "gompertz_4",
  #  method       = "nls",
  #  formula_expr = quote(y ~ gomp_4(x, a, b, c, h)),
  #  start        = c(a = 1, b =10.98, c = -1.44, h = 0)
  #)
)

model_specs_2_miknn_res <- fit_many_models(
  data       = dat_allsstar,
  x_col      = "avg_miknn_score",
  y_col      = "avg_wcr",
  model_list = model_specs_2_miknn,
  n          = nrow(dat_allsstar)
)

model_specs_2_miknn_res$summary
model_specs_2_miknn_res$coefs
ggplot(model_specs_2_miknn_res$predictions, aes(x = avg_miknn_score, y = avg_wcr)) +
  geom_point(data = dat_allsstar, aes(x = avg_miknn_score, y = avg_wcr)) +
  geom_line(aes(y = FittedValue)) +
  facet_wrap(~Model, scales = "free_y")


model_specs_2_siib <- list(
  #linear = list(
  #  method       = "lm",
  #  formula_expr = quote(y ~ x)  # 'quote()' so we can substitute later
  #),
  f_siib = list(
    name = "f_siib",
    method       = "nlsLM",
    formula_expr = quote(y ~ f_siib(x, a, b)),
    start        = c(a = 0.01, b = 1)
  ),
  f_siib3 = list(
    name = "f_siib3",
    method       = "nlsLM",
    formula_expr = quote(y ~ f_siib3(x, a, b, c)),
    start        = c(a = 0.0152, b = 1.2799, c = 1)
  ),
  monpoly3 = list(
    name = "monpoly3",
    method       = "monpol",
    formula_expr = quote(y ~ x),
    degree       = 3
  ),
  logistic = list(
    name = "logistic",
    method       = "nlsLM",
    formula_expr = quote(y ~ logistic_taal10(x, a, b)),
    start        = c(a = 0, b = -0.5)
  ),
  logistic_3 = list(
    name = "logistic_3",
    method       = "nlsLM",
    formula_expr = quote(y ~ logistic_gen(x, a, b, c)),
    start        = c(a = -0.0722, b = 1.5338 , c=0), # Actual 0.05724
    lower = c(a = -Inf, b = -Inf, c = 0),
    upper = c(a = Inf, b = Inf, c = 0.1)
  ),
  dau = list(
    name = "dau",
    method       = "nlsLM",
    formula_expr = quote(y ~ dau_taal10(x, a, b, c)),
    start        = c(a = 0.0001, b = 1, c = -1)
  ) ,
  gompertz = list(
    name = "gompertz",
    method       = "nlsLM",
    formula_expr = quote(y ~ gomp(x, b, c)),
    start        = c(b=0.05, c=-0.01)
  ),
  gompertz_3 = list(
  
      name = "gompertz_3",
    method       = "nlsLM",
    formula_expr = quote(y ~ gomp_3(x,a, b, c)),
    start        = c(a = 1, b=0.466771252, c=-0.009882308)
  )
  #gompertz_4 = list(
  #  name = "gompertz_4",
  #  method       = "nls",
  #  formula_expr = quote(y ~ gomp_4(x, a, b, c, h)),
  #  start        = c(a = 1, b =0.00909, c = -0.1149, h = 0)
  #)
)

model_specs_2_siib_res <- fit_many_models(
  data       = dat_allsstar,
  x_col      = "avg_siib_score",
  y_col      = "avg_wcr",
  model_list = model_specs_2_siib,
  n          = nrow(dat_allsstar)
)

model_specs_2_siib_res$summary
model_specs_2_siib_res$coefs
ggplot(model_specs_2_siib_res$predictions, aes(x = avg_siib_score, y = avg_wcr)) +
  geom_point(data = dat_allsstar, aes(x = avg_siib_score, y = avg_wcr)) +
  geom_line(aes(y = FittedValue)) +
  facet_wrap(~Model, scales = "free_y")


model_specs_2_snr <- list(
  linear = list(
    method       = "lm",
    formula_expr = quote(y ~ x)  # 'quote()' so we can substitute later
  ),
  #f_siib = list(
  #  name = "f_siib",
  #  method       = "nlsLM",
  #  formula_expr = quote(y ~ f_siib(x, a, b)),
  #  start        = c(a = 0.0152, b = 1.2799)
  #),
  #f_siib3 = list(
  #  name = "f_siib3",
  #  method       = "nlsLM",
  #  formula_expr = quote(y ~ f_siib3(x, a, b, c)),
  #  start        = c(a = 0.0152, b = 1.2799, c = 1)
  #),
  monpoly3 = list(
    name = "monpoly3",
    method       = "monpol",
    formula_expr = quote(y ~ x),
    degree       = 3
  ),
  logistic = list(
    name = "logistic",
    method       = "nlsLM",
    formula_expr = quote(y ~ logistic_taal10(x, a, b)),
    start        = c(a = 0, b = -0.5)
  ),
  logistic_3 = list(
    name = "logistic_3",
    method       = "nlsLM",
    formula_expr = quote(y ~ logistic_gen(x, a, b, c)),
    start        = c(a = -0.282, b = -1.428 , c=0), # Actual 0.05724
    lower = c(a = -Inf, b = -Inf, c = 0),
    upper = c(a = Inf, b = Inf, c = 0.1)
  ),  
  logistic_3_taal09 = list(
    name = "logistic_3_taal09",
    method       = "nlsLM",
    formula_expr = quote(y ~ logistic_3_taal09(x, a, b, c)),
    start        = c(a = 0, b = -0.5, c=0.1),
    lower = c(a = -Inf, b = -Inf, c = 0),
    upper = c(a = Inf, b = Inf, c = 0.05)
  ),
  dau = list(
    name = "dau",
    method       = "nlsLM",
    formula_expr = quote(y ~ dau_taal10(x, a, b, c)),
    start        = c(a = 0.0001, b = 1, c = -1)
  ) ,
  gompertz = list(
    name = "gompertz",
    method       = "nlsLM",
    formula_expr = quote(y ~ gomp(x, b, c)),
    start        = c(b=0.05, c=-0.01)
  ),
  gompertz_3 = list(
    
    name = "gompertz_3",
    method       = "nlsLM",
    formula_expr = quote(y ~ gomp_3(x,a, b, c)),
    start        = c(a = 1, b=0.466771252, c=-0.009882308)
  )
  #gompertz_4 = list(
  #  name = "gompertz_4",
  #  method       = "nls",
  #  formula_expr = quote(y ~ gomp_4(x, a, b, c, h)),
  #  start        = c(a = 1, b =0.00909, c = -0.1149, h = 0)
  #)
)

model_specs_2_snr_res <- fit_many_models(
  data       = dat_allsstar,
  x_col      = "snr",
  y_col      = "avg_wcr",
  model_list = model_specs_2_snr,
  n          = nrow(dat_allsstar)
)

model_specs_2_snr_res$summary
model_specs_2_snr_res$coefs
ggplot(model_specs_2_snr_res$predictions, aes(x = snr, y = avg_wcr)) +
  geom_point(data = dat_allsstar, aes(x = snr, y = avg_wcr)) +
  geom_line(aes(y = FittedValue)) +
  facet_wrap(~Model, scales = "free_y")

model_specs_2_snr_stoi_res <- fit_many_models(
  data       = dat_allsstar,
  x_col      = "snr",
  y_col      = "avg_stoi_score",
  model_list = model_specs_2_snr,
  n          = nrow(dat_allsstar)
)

model_specs_2_snr_stoi_res$summary
model_specs_2_snr_stoi_res$coefs
ggplot(model_specs_2_snr_stoi_res$predictions, aes(x = snr, y = avg_stoi_score)) +
  geom_point(data = dat_allsstar, aes(x = snr, y = avg_stoi_score)) +
  geom_line(aes(y = FittedValue)) +
  facet_wrap(~Model, scales = "free_y")



model_specs_2_snr_siib <- list(
  linear = list(
    method       = "lm",
    formula_expr = quote(y ~ x)  # 'quote()' so we can substitute later
  ),
  #f_siib = list(
  #  name = "f_siib",
  #  method       = "nlsLM",
  #  formula_expr = quote(y ~ f_siib(x, a, b)),
  #  start        = c(a = 0.0152, b = 1.2799)
  #),
  #f_siib3 = list(
  #  name = "f_siib3",
  #  method       = "nlsLM",
  #  formula_expr = quote(y ~ f_siib3(x, a, b, c)),
  #  start        = c(a = 0.0152, b = 1.2799, c = 1)
  #),
  monpoly3 = list(
    name = "monpoly3",
    method       = "monpol",
    formula_expr = quote(y ~ x),
    degree       = 3
  ),
  logistic = list(
    name = "logistic",
    method       = "nlsLM",
    formula_expr = quote(y ~ logistic_taal10(x, a, b)),
    start        = c(a = 1, b = 0)
  ),
  logistic_3 = list(
    name = "logistic_3",
    method       = "nlsLM",
    formula_expr = quote(y ~ logistic_gen(x, a, b, c)),
    start        = c(a = -0.282, b = -1.428 , c=0), # Actual 0.05724
    lower = c(a = -Inf, b = -Inf, c = 0),
    upper = c(a = Inf, b = Inf, c = 0.1)
  ),  
  dau = list(
    name = "dau",
    method       = "nlsLM",
    formula_expr = quote(y ~ dau_taal10(x, a, b, c)),
    start        = c(a = 0.0001, b = 1, c = -1)
  ) ,
  gompertz = list(
    name = "gompertz",
    method       = "nlsLM",
    formula_expr = quote(y ~ gomp(x, b, c)),
    start        = c(b=0.05, c=-0.01)
  ),
  gompertz_3 = list(
    
    name = "gompertz_3",
    method       = "nlsLM",
    formula_expr = quote(y ~ gomp_3(x,a, b, c)),
    start        = c(a = 1, b=0.466771252, c=-0.009882308)
  )
  #gompertz_4 = list(
  #  name = "gompertz_4",
  #  method       = "nls",
  #  formula_expr = quote(y ~ gomp_4(x, a, b, c, h)),
  #  start        = c(a = 1, b =0.00909, c = -0.1149, h = 0)
  #)
)

model_specs_2_snr_siib_res <- fit_many_models(
  data       = dat_allsstar,
  x_col      = "snr",
  y_col      = "avg_siib_score",
  model_list = model_specs_2_snr_siib,
  n          = nrow(dat_allsstar)
)

model_specs_2_snr_siib_res$summary
model_specs_2_snr_siib_res$coefs
ggplot(dat_allsstar, aes(x = snr, y = avg_siib_score)) +
  geom_point(data = dat_allsstar, aes(x = snr, y = avg_siib_score))



model_specs_2_siib_gauss <- list(
  #linear = list(
  #  method       = "lm",
  #  formula_expr = quote(y ~ x)  # 'quote()' so we can substitute later
  #),
  f_siib = list(
    name = "f_siib",
    method       = "nlsLM",
    formula_expr = quote(y ~ f_siib(x, a, b)),
    start        = c(a = 0.01, b = 1)
  ),
  f_siib3 = list(
    name = "f_siib3",
    method       = "nlsLM",
    formula_expr = quote(y ~ f_siib3(x, a, b, c)),
    start        = c(a = 0.0152, b = 1.2799, c = 1)
  ),
  monpoly3 = list(
    name = "monpoly3",
    method       = "monpol",
    formula_expr = quote(y ~ x),
    degree       = 3
  ),
  logistic = list(
    name = "logistic",
    method       = "nlsLM",
    formula_expr = quote(y ~ logistic_taal10(x, a, b)),
    start        = c(a = 0, b = -0.5)
  ),
  logistic_3 = list(
    name = "logistic_3",
    method       = "nlsLM",
    formula_expr = quote(y ~ logistic_gen(x, a, b, c)),
    start        = c(a = -0.0722, b = 1.5338 , c=0), # Actual 0.05724
    lower = c(a = -Inf, b = -Inf, c = 0),
    upper = c(a = Inf, b = Inf, c = 0.1)
  ),
  dau = list(
    name = "dau",
    method       = "nlsLM",
    formula_expr = quote(y ~ dau_taal10(x, a, b, c)),
    start        = c(a = 0.0001, b = 1, c = -1)
  ) ,
  gompertz = list(
    name = "gompertz",
    method       = "nlsLM",
    formula_expr = quote(y ~ gomp(x, b, c)),
    start        = c(b=0.05, c=-0.01)
  ),
  gompertz_3 = list(
    
    name = "gompertz_3",
    method       = "nlsLM",
    formula_expr = quote(y ~ gomp_3(x,a, b, c)),
    start        = c(a = 1, b=0.466771252, c=-0.009882308)
  )
  #gompertz_4 = list(
  #  name = "gompertz_4",
  #  method       = "nls",
  #  formula_expr = quote(y ~ gomp_4(x, a, b, c, h)),
  #  start        = c(a = 1, b =0.00909, c = -0.1149, h = 0)
  #)
)

model_specs_2_siib_gauss_res <- fit_many_models(
  data       = dat_allsstar,
  x_col      = "avg_siib_score_gaussian",
  y_col      = "avg_wcr",
  model_list = model_specs_2_siib_gauss,
  n          = nrow(dat_allsstar)
)

model_specs_2_siib_gauss_res$summary
model_specs_2_siib_gauss_res$predictions

model_specs_2_siib_gauss_res$coefs
ggplot(model_specs_2_siib_gauss_res$predictions, aes(x = avg_siib_score_gaussian, y = avg_wcr)) +
  geom_point(data = dat_allsstar, aes(x = avg_siib_score_gaussian, y = avg_wcr)) +
  geom_line(aes(y = FittedValue)) +
  facet_wrap(~Model, scales = "free_y")

# SNR to STOI graph
ggplot(dat_allsstar, aes(x = snr, y = avg_siib_score)) +
  geom_point(aes(color = avg_wcr)) +
  geom_smooth(method = "lm") +
  scale_color_viridis_c()




fit_many_models <- function(data, x_col, y_col, model_list, n = nrow(data)) {
  # data: your data frame
  # x_col, y_col: strings for the x and y columns in `data`
  # model_list: the named list of model specs
  # n: number of data points (used in AICc)
  
  # We’ll store:
  # 1) A data frame of model performance metrics
  # 2) A data frame of predictions (in "long" format for ggplot)
  
  # 1. Prepare an empty container for each model’s summary
  model_summaries <- list()
  
  # 2. Prepare a container for predicted curves
  #    We'll build it as we go, or initialize after we have predictions.
  
  # Build a formula-friendly expression from user’s x_col, y_col
  # We'll do something like:  y_col ~ x_col  inside each specification
  # but the user might specify something more advanced (poly, etc.)
  # So we’ll do it via substituting x->(x_col) and y->(y_col) in the quoted expression.
  
  # This function: Substitutes the symbol 'x' with x_col, 'y' with y_col, etc.
  make_formula <- function(expr, x_sym, y_sym) {
    # e.g., expr = quote(y ~ log(x)), x_sym=as.symbol("avg_stoi_score"), etc.
    # Then we evaluate it into an actual formula object.
    call_mod <- eval(substitute(
      substitute(EXPR, list(x = x_sym, y = y_sym)), 
      list(EXPR = expr)
    ))
    # call_mod is now something like `y_col ~ log(x_col)` as a call
    # Turn it into a formula
    return(call_mod)
  }
  
  # 3. We will define a grid of x-values for plotting predictions
  x_vals <- seq(min(data[[x_col]], na.rm = TRUE),
                max(data[[x_col]], na.rm = TRUE),
                length.out = 200)
  
  # Build a “prediction data frame” with the same column name as `x_col`
  pred_df <- data.frame(x_vals)
  names(pred_df) <- x_col
  
  coefs <- list()
  
  # We’ll build up columns in pred_df for each model’s predicted Y.
  
  for (model_name in names(model_list)) {
    print(paste("Model:", model_name))
    spec <- model_list[[model_name]]
    
    
    # We expect:
    #   spec$method in c("lm","nls","nlsLM","monpol", etc.)
    #   spec$formula_expr is a quoted expression like y ~ x, y ~ logistic_taal10(x,a,b), etc.
    #   Maybe spec$degree for monpol
    #   Maybe spec$start for nls
    
    # 1. Construct the formula from spec$formula_expr
    fmla <- make_formula(spec$formula_expr,
                         x_sym = as.symbol(x_col),
                         y_sym = as.symbol(y_col))
    
    # 2. Fit the model
    convert_monpol <- function(monpol_model) {
      l <- list(
        coefficients = coef(monpol_model),  # If available
        fitted.values = monpol_model$fitted,  # Adjust if necessary
        residuals = monpol_model$residuals  # Adjust if necessary
      )
      attr(l, "class") <- "monpol"
   #   attr(l, "names") <- attr(monpol_model, "names")
      l.function.name <- "monpol"
      return(l)
    }
    
    fit_obj <- NULL
    if (spec$method == "lm") {
      fit_obj <- lm(fmla, data = data)
    } else if (spec$method == "monpol") {
      # monpol uses "fmla", "degree", and possibly other arguments
      fit_obj <- monpol(fmla, data = data, degree = spec$degree)
      #fit_obj <- convert_monpol(fit_obj)
    } else if (spec$method == "nls") {
      
        fit_obj <- nls(fmla, data = data, start = spec$start, control = nls.control(maxiter = 200))
      
    } else if (spec$method == "nlsLM") {
      if (spec$name == "gompertz" || spec$name == "logistic_3" || spec$name == "logistic_3_wrong" || spec$name == "logistic_4") {
        fit_obj <- nlsLM(fmla, data = data, control = nls.lm.control(maxiter = 200), start = spec$start, lower = spec$lower, upper = spec$upper) 
      } else {
      fit_obj <- nlsLM(fmla, data = data, start = spec$start)
      }
    } else {
      stop("Unsupported fitting method: ", spec$method)
    }
    
    
    if (!is.null(fit_obj$converged) && !fit_obj$converged) {
      stop("Error: The model did not converge. Check starting values and constraints.")
    }
    
    print(paste("!!!!", spec$method, class(fit_obj)))
    #print(paste("Cv:", spec$method, str(cv::cv(fit_obj))))
    print(paste("Coefs:", model_name, spec$method, str(coef(fit_obj))))
    coefs[[model_name]] <- coef(fit_obj)
    
    
    
    # 3. Compute metrics: RMSE, AIC, AICc, PCC, etc.
    #    We'll define a small helper for predictions
    pred_obj <- predict(fit_obj, newdata = data)
    # monpol returns a data frame with columns x and y
    # others return a numeric vector
    if (spec$method == "monpol") {
      preds <- c(pred_obj)
      obs <- data[[y_col]]  # monpol returns a data frame
    } else {
      preds <- pred_obj
      obs <- data[[y_col]]
    } 
    
    
    
    rmse <- sqrt(mean((obs - preds)^2, na.rm = TRUE))
    aic  <- AIC(fit_obj)
    bic <- BIC(fit_obj)
    #anova <- car::Anova(fit_obj)
    #print(paste("Anova:", spec$method, str(anova)))
    
    # qq-plot
    r <- residuals(fit_obj)
    kt <- ks.test(r, "pnorm", mean=mean(r), sd=sd(r))
    s <- shapiro.test(r)
    print(paste("Shapiro:", spec$method, str(s)))
    print(paste("KS:", spec$method, str(kt)))
    

    qqnorm(residuals(fit_obj))
    qqline(residuals(fit_obj))
    
    
    
    # For number of parameters k, we can approximate with length(coef(fit_obj)).
    k    <- length(coef(fit_obj))
    aicc <- AICc(fit_obj, n = n, k = k)
    
    pcc <- cor(obs, preds, use = "complete.obs")  # Pearson corr
    # Kendall’s tau (optional)
    ktau <- cor(obs, preds, method = "kendall", use = "complete.obs")
    roundd <- function(x, digits = 3) {
      return(round(x, digits))
    }
    # Collect results
    model_summaries[[model_name]] <- data.frame(
      Model       = model_name,
      Params      = k,
      Method      = spec$method,
      RMSE        = roundd(rmse),
      PCC         = roundd(pcc),
      AICc        = roundd(aicc),
      stringsAsFactors = FALSE
    )
    
    
    # 4. Add predictions for the *new grid* of x-values
    pred_colname <- paste0(model_name, "_pred")
    pred_df[[pred_colname]] <- predict(fit_obj, newdata = pred_df)
  }
  
  # Combine the model summaries
  summary_df <- do.call(rbind, model_summaries)
  
  # We can also compute the minimum AICc and get ΔAICc, etc.
  min_aicc <- min(summary_df$AICc)
  min_aic <- min(summary_df$AIC)
  summary_df <- summary_df %>%
    mutate(
      
      dAICc = AICc - min_aicc,
      relLikAICc_raw = exp(-0.5 * dAICc),
      relLikAICc = round(exp(-0.5 * dAICc), 3),
      weightAICc = round(relLikAICc_raw / sum(relLikAICc_raw), 3)
    )
  
  # Now pivot the pred_df from wide to long for plotting
  # pred_df has columns: x_col, model1_pred, model2_pred, ...
  long_pred_df <- pred_df %>%
    tidyr::pivot_longer(
      cols = -all_of(x_col),
      names_to = "Model",
      values_to = "FittedValue"
    )
  
  # Return a list with everything we might want
  list(
    summary    = summary_df,
    predictions = long_pred_df,
    coefs = coefs
  )
}

# Example usage:
dantale_stoi_res <- fit_many_models(
  data       = dantale_all,
  x_col      = "avg_stoi_score",
  y_col      = "avg_wcr",
  model_list = model_specs,
  n          = nrow(dantale_all)
)

all_conds_stoi_res <- fit_many_models(
  data       = dat_all,
  x_col      = "avg_stoi_score",
  y_col      = "avg_wcr",
  model_list = model_specs,
  n          = nrow(dat_all)
)

dantale_conds_stoi_res <- fit_many_models(
  data       = dantale_conds,
  x_col      = "avg_stoi_score",
  y_col      = "avg_wcr",
  model_list = model_specs,
  n          = nrow(dantale_conds)
)

#allsstar_stoi_res <- fit_many_models(
#  data       = allsstar_data,
#  x_col      = "stoi_score",
#  y_col      = "sim_wcr",
#  model_list = model_specs,
#  n          = nrow(allsstar_data)
#)

##

dat_n <- dat_allsstar
dat_n$snr <- as.numeric(as.character(dat_n$snr))

dat_allsstar_snr_wcr_res <- fit_many_models(
  data       = dat_n,
  x_col      = "snr",
  y_col      = "avg_wcr",
  model_list = model_specs_1,
  n          = nrow(dat_allsstar)
)

dat_allsstar_snr_wcr_res$summary
dat_allsstar_snr_wcr_res$coefs
ggplot(dat_allsstar_snr_wcr_res$predictions, aes(x = snr, y = avg_wcr)) +
  geom_point(data = dat_allsstar, aes(x = snr, y = avg_wcr)) +
  geom_line(aes(y = FittedValue)) +
  facet_wrap(~Model, scales = "free_y")

dat_allsstar_stoi_res <- fit_many_models(
  data       = dat_allsstar,
  x_col      = "avg_stoi_score",
  y_col      = "avg_wcr",
  model_list = model_specs_1_stoi_log,
  n          = nrow(dat_allsstar)
)

dat_allsstar_stoi_res$summary
dat_allsstar_stoi_res$coefs
ggplot(dat_allsstar_stoi_res$predictions, aes(x = avg_stoi_score, y = avg_wcr)) +
  geom_point(data = dat_allsstar, aes(x = avg_stoi_score, y = avg_wcr)) +
  geom_line(aes(y = FittedValue)) +
  facet_wrap(~Model, scales = "free_y")


dat_allsstar_miknn_res <- fit_many_models(
  data       = dat_allsstar,
  x_col      = "avg_miknn_score",
  y_col      = "avg_wcr",
  model_list = model_specs_1_miknn_log,
  n          = nrow(dat_allsstar)
)

dat_allsstar_miknn_res$summary
dat_allsstar_miknn_res$coefs
ggplot(dat_allsstar_miknn_res$predictions, aes(x = avg_miknn_score, y = avg_wcr)) +
  geom_point(data = dat_allsstar, aes(x = avg_miknn_score, y = avg_wcr)) +
  geom_line(aes(y = FittedValue)) +
  facet_wrap(~Model, scales = "free_y")


dat_allsstar_siib_res  <- fit_many_models(
  data       = dat_allsstar,
  x_col      = "avg_siib_score",
  y_col      = "avg_wcr",
  model_list = model_specs,
  n          = nrow(dat_allsstar)
)

dat_allsstar_siib_res$summary
dat_allsstar_siib_res$coefs
ggplot(dat_allsstar_siib_res$predictions, aes(x = avg_siib_score, y = avg_wcr)) +
  geom_point(data = dat_allsstar, aes(x = avg_siib_score, y = avg_wcr)) +
  geom_line(aes(y = FittedValue)) +
  facet_wrap(~Model, scales = "free_y")



dat_allsstar_stoi_log_res <- fit_many_models(
  data       = dat_allsstar,
  x_col      = "avg_stoi_score",
  y_col      = "avg_wcr",
  model_list = model_specs_1_siib_gaus_log,
  n          = nrow(dat_allsstar)
)

dat_allsstar_miknn_log_res <- fit_many_models(
  data       = dat_allsstar,
  x_col      = "avg_miknn_score",
  y_col      = "avg_wcr",
  model_list = model_specs_1,
  n          = nrow(dat_allsstar)
)

dat_allsstar_miknn_res <- fit_many_models(
  data       = dat_allsstar,
  x_col      = "avg_miknn_score",
  y_col      = "avg_wcr",
  model_list = model_specs_miknn,
  n          = nrow(dat_allsstar)
)

dat_allsstar_siib_gaus_res  <- fit_many_models(
  data       = dat_allsstar,
  x_col      = "avg_siib_score_gaussian",
  y_col      = "avg_wcr",
  model_list = model_specs_1_siib_gaus_log,
  n          = nrow(dat_allsstar)
)

dat_allsstar_siib_res  <- fit_many_models(
  data       = dat_allsstar,
  x_col      = "avg_siib_score",
  y_col      = "avg_wcr",
  model_list = model_specs_1,
  n          = nrow(dat_allsstar)
)

dat_allsstar_siib_res_all <- fit_many_models(
  data       = dat_allsstar,
  x_col      = "avg_siib_score",
  y_col      = "avg_wcr",
  model_list = model_specs,
  n          = nrow(dat_allsstar)
)

dat_allsstar_siib_res_all$summary
dat_allsstar_siib_res_all$coefs
ggplot(dat_allsstar_siib_res_all$predictions, aes(x = avg_siib_score, y = avg_wcr)) +
  geom_point(data = dat_allsstar, aes(x = avg_siib_score, y = avg_wcr)) +
  geom_line(aes(y = FittedValue)) +
  facet_wrap(~Model, scales = "free_y")

dat_allsstar_siib_res_2  <- fit_many_models(
  data       = dat_allsstar,
  x_col      = "avg_siib_score",
  y_col      = "avg_wcr",
  model_list = model_specs_siib,
  n          = nrow(dat_allsstar)
)

dantale_miknn_res <- fit_many_models(
  data       = dantale_all,
  x_col      = "avg_miknn_score",
  y_col      = "avg_wcr",
  model_list = model_specs,
  n          = nrow(dantale_all)
)

dantale_conds_miknn_res <- fit_many_models(
  data       = dantale_conds,
  x_col      = "avg_miknn_score",
  y_col      = "avg_wcr",
  model_list = model_specs,
  n          = nrow(dantale_conds)
)


# Check the summary:

#        Model   Method     RMSE      AIC     AICc       PCC KendallTau    dAICc     relLik     weight
# 1     linear       lm   0.0852  10.1234  10.5811  ...   ...  
# 2   monpoly3   monpol   0.0821   9.9981  10.2234  ...   ...
# 3   logistic      nls   0.0760   6.1235   7.0012  ...   ...
# 4  gompertz      nls    0.0778   7.5233   7.9900  ...   ...
# etc.

# For plotting:
long_pred <- allsstar_stoi_res$predictions


data_sum  = allsstar_data %>% group_by(sim_wcr) %>% summarize(stoi_score = mean(stoi_score), snr = mean(snr))



all_conds_stoi_res$summary
ggplot(all_conds_stoi_res$predictions, aes(x = avg_stoi_score, y = avg_wcr)) +
  geom_point(data = dat_all, aes(x = avg_stoi_score, y = avg_wcr)) +
  geom_line(aes(y = FittedValue)) +
  facet_wrap(~Model, scales = "free_y")

# plot res
allsstar_stoi_res$summary
ggplot(long_pred, aes(x = stoi_score, y = sim_wcr)) +
  geom_point(data = allsstar_data, aes(x = stoi_score, y = sim_wcr)) +
  geom_line(aes(y = FittedValue)) +
  facet_wrap(~Model, scales = "free_y")

###
dat_allsstar_stoi_res$summary
dat_allsstar_stoi_res$coefs
ggplot(dat_allsstar_stoi_res$predictions, aes(x = avg_stoi_score, y = avg_wcr)) +
  geom_point(data = dat_allsstar, aes(x = avg_stoi_score, y = avg_wcr)) +
  geom_line(aes(y = FittedValue)) +
  facet_wrap(~Model, scales = "free_y")

# 1
dat_allsstar_stoi_log_res$summary
dat_allsstar_stoi_log_res$coefs
ggplot(dat_allsstar_stoi_log_res$predictions, aes(x = avg_stoi_score, y = avg_wcr)) +
  geom_point(data = dat_allsstar, aes(x = avg_stoi_score, y = avg_wcr)) +
  geom_line(aes(y = FittedValue)) +
  facet_wrap(~Model, scales = "free_y")
# 1

#
dat_allsstar_miknn_log_res$summary
dat_allsstar_miknn_log_res$coefs
ggplot(dat_allsstar_miknn_log_res$predictions, aes(x = avg_miknn_score, y = avg_wcr)) +
  geom_point(data = dat_allsstar, aes(x = avg_miknn_score, y = avg_wcr)) +
  geom_line(aes(y = FittedValue)) +
  facet_wrap(~Model, scales = "free_y")


dat_allsstar_siib_gaus_res$summary
dat_allsstar_siib_gaus_res$coefs
ggplot(dat_allsstar_siib_gaus_res$predictions, aes(x = avg_siib_score_gaussian, y = avg_wcr)) +
  geom_point(data = dat_allsstar, aes(x = avg_siib_score_gaussian, y = avg_wcr)) +
  geom_line(aes(y = FittedValue)) +
  facet_wrap(~Model, scales = "free_y")

dat_allsstar_miknn_res$summary
dat_allsstar_miknn_res$coefs
ggplot(dat_allsstar_miknn_res$predictions, aes(x = avg_miknn_score, y = avg_wcr)) +
  geom_point(data = dat_allsstar, aes(x = avg_miknn_score, y = avg_wcr)) +
  geom_line(aes(y = FittedValue)) +
  facet_wrap(~Model, scales = "free_y")

dat_allsstar_siib_res$summary
dat_allsstar_siib_res$coefs
ggplot(dat_allsstar_siib_res$predictions, aes(x = avg_siib_score, y = avg_wcr)) +
  geom_point(data = dat_allsstar, aes(x = avg_siib_score, y = avg_wcr)) +
  geom_line(aes(y = FittedValue)) +
  facet_wrap(~Model, scales = "free_y")

dat_allsstar_siib_res_2$summary
dat_allsstar_siib_res_2$coefs
ggplot(dat_allsstar_siib_res_2$predictions, aes(x = avg_siib_score, y = avg_wcr)) +
  geom_point(data = dat_allsstar, aes(x = avg_siib_score, y = avg_wcr)) +
  geom_line(aes(y = FittedValue)) +
  facet_wrap(~Model, scales = "free_y")

# plot res

dantale_stoi_res$summary
ggplot(dantale_stoi_res$predictions, aes(x = avg_stoi_score, y = avg_wcr)) +
  geom_point(data = dantale_all, aes(x = avg_stoi_score, y = avg_wcr)) +
  geom_line(aes(y = FittedValue)) +
  facet_wrap(~Model, scales = "free_y")


dantale_conds_stoi_res$summary
ggplot(dantale_conds_stoi_res$predictions, aes(x = avg_stoi_score, y = avg_wcr)) +
  geom_point(data = dantale_conds, aes(x = avg_stoi_score, y = avg_wcr)) +
  geom_line(aes(y = FittedValue)) +
  facet_wrap(~Model, scales = "free_y")


dat_allsstar_stoi_res$summary
ggplot(dat_allsstar_stoi_res$predictions, aes(x = avg_stoi_score, y = avg_wcr)) +
  geom_point(data = dat_allsstar, aes(x = avg_stoi_score, y = avg_wcr)) +
  geom_line(aes(y = FittedValue)) +
  facet_wrap(~Model, scales = "free_y")

dantale_miknn_res$summary
ggplot(dantale_miknn_res$predictions, aes(x = avg_miknn_score, y = avg_wcr)) +
  geom_point(data = dantale_all, aes(x = avg_miknn_score, y = avg_wcr)) +
  geom_line(aes(y = FittedValue)) +
  facet_wrap(~Model, scales = "free_y")

dantale_conds_miknn_res$summary
ggplot(dantale_conds_miknn_res$predictions, aes(x = avg_miknn_score, y = avg_wcr)) +
  geom_point(data = dantale_conds, aes(x = avg_miknn_score, y = avg_wcr)) +
  geom_line(aes(y = FittedValue)) +
  facet_wrap(~Model, scales = "free_y")

##!!!!!
ggplot() +
  geom_point(
    data  = data_sum,
    aes(x = sim_wcr, y = stoi_score, color = snr)
  ) 


ggplot() +
  geom_point(
    data  = data_sum,
    aes(x = snr, y = sim_wcr, color = stoi_score)
  ) 


#



