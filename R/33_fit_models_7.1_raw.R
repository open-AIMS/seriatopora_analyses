## Retrieve the data
load(file = "../data/processed/data_q2.R")


## Strip out all unnecessary variables

full_data <- data
data <- full_data

## Prepare data

data <- full_data
## ---- q2_7_1_raw_prepare_data
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
  mutate(across(c(s, c, d, b, u), \(x) ifelse(Dist.time == "Before", 0, x))) |> 
  mutate(across(c(s, c, d, b, u), factor)) 
 
save(data, file = "../data/modelled/q2_data_raw_1.R") 
## ----end


## ---- 7_1_raw_cellmeans
load(file = "../data/modelled/q2_data_raw_1.R") 
cellmeans_summ_raw <-
  data |>
  mutate(Dist.time2 = ifelse(Dist == "Before", "Before", "After")) |> 
  group_by(A_SECTOR, SHELF, Dist, Dist.time) |>
  reframe(
    type = c("mean", "median"),
    Mean = c(mean(n.points / total.points), median(n.points / total.points)),
    SD = c(sd(n.points / total.points), MAD = mad(n.points / total.points)),
    N = c(n(), n())
  ) |>
  mutate(
    lower = Mean - 2 * (SD / sqrt(N)),
    upper = Mean + 2 * (SD / sqrt(N))
  ) 

save(cellmeans_summ_raw, file = "../data/modelled/cellmeans_summ_raw_7.1.RData")
cellmeans_summ_raw |>
  ggplot(aes(y = Mean, x = Dist, colour = Dist.time, shape = type)) +
  geom_pointrange(aes(ymin = lower, ymax = upper)) +
  scale_shape_manual(values = c(16, 21)) +
  facet_wrap(A_SECTOR ~ SHELF,
    scales = "free",
    labeller = label_wrap_gen(multi_line = FALSE)
  )
## ----end

## ---- 7_1_raw_effects
load(file = "../data/modelled/cellmeans_summ_raw_7.1.RData")
eff_summ_raw <- cellmeans_summ_raw |>
  mutate(Values = Mean) |>
  nest(.by = c(A_SECTOR, SHELF, type)) |>
  mutate(eff = map(
    .x = data,
    .f = ~ {
      before_vs_afters(.x)
    }
  )) |>
  dplyr::select(-data) |> 
  unnest(c(eff))
save(eff_summ_raw, file = "../data/modelled/eff_summ_raw_7.1.RData")
eff_summ_raw |> ggplot(aes(y = Dist, x = Values, colour = type)) +
  geom_vline(xintercept = 1, linetype = "dashed") +
  geom_point() +
  facet_grid(A_SECTOR ~ SHELF, scales = "free")
## ----end
