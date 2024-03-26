## Retrieve the data
load(file = "../data/processed/data_q2.R")


## Strip out all unnecessary variables

full_data <- data
data <- full_data

## Prepare data

data <- full_data
## ---- q2_7_1_glmmTMB_prepare_data
## Focus on only the necessary variables
load(file = "../data/processed/data_q2.R")
data <- data |>
  mutate(SecShelf = paste(A_SECTOR, SHELF)) |> 
  dplyr::select(
    n.points, total.points, Dist.time, s, c, d, b, u,
    AIMS_REEF_NAME, Site, Transect, A_SECTOR, SHELF, SecShelf, Dist.number,
    SecShelf, Dist, SSDist
  ) |>
  filter(!SecShelf %in% c("CG M", "PC M", "PC O")) |>
  mutate(across(c(s, c, d, b, u), \(x) if_else(Dist.time == "Before", "0", x))) |> 
  mutate(across(c(s, c, d, b, u), factor)) 
 
save(data, file = "../data/modelled/q2_data_raw_1.R") 
## ----end

## ---- q2_7_1_glmmTMB_fit_model
load(file = "../data/modelled/q2_data_raw_1.R") 
if (rerun) {
  mod_glmmTMB <- glmmTMB(
    cbind(n.points, total.points - n.points) ~ 1 + (s + c + d + b + u ) +
      ## (1 | Dist.time:SecShelf:(s + c + d + b + u)) +
      ## (1 | SecShelf:(s + c + d + b + u)) +
      (1 | AIMS_REEF_NAME) +
      (1 | Site) +
      (1 | Transect),
    ziformula = ~ 1 + (1 | AIMS_REEF_NAME) +
      (1 | Site) +
      (1 | Transect),
    data = data,
    family = "binomial",
    REML = TRUE
  )

  summary(mod_glmmTMB)
  save(mod_glmmTMB, file = "../data/modelled/mod_glmmTMB_7.1.RData")
}

## ----end

## ---- g2_7_1_glmmTMB_cellmeans
mod_glmmTMB |> emmeans(~s + c + d + b + u, type = "response")
## ----end

load(file = "../data/modelled/q2_data_raw_1.R") 
if (rerun) {
  mod_glmmTMB <- glmmTMB(
    cbind(n.points, total.points - n.points) ~ 1 + SSDist +
      ## (1 | Dist.time:SecShelf:(s + c + d + b + u)) +
      ## (1 | SecShelf:(s + c + d + b + u)) +
      (1 | AIMS_REEF_NAME) +
      (1 | Site) +
      (1 | Transect),
    ziformula = ~ 1 + (1 | AIMS_REEF_NAME) +
      (1 | Site) +
      (1 | Transect),
    data = data,
    family = "binomial",
    REML = TRUE
  )

  summary(mod_glmmTMB)
  save(mod_glmmTMB, file = "../data/modelled/mod_glmmTMB_7.1.RData")
}
mod_glmmTMB |> emmeans(~s + c + d + b + u, type = "response")
