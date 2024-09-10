
## Retrieve the data
## ---- q3_retrieve_data
load(file = "../data/processed/data_q3.RData")
## ----end


## Abundance

## Explore the spatio-temporal domain using simple means and CI's
## Note in particular, the CI's are completely invalid
## ---- q3_eda_1
data |>
  group_by(A_SECTOR, SHELF, Site, cREPORT_YEAR) |>
  summarise(across(juv.dens.aa, list(
    mu = mean, med = median,
    sd = sd, n = length
  ))) |>
  mutate(
    lower = juv.dens.aa_mu - 2 * juv.dens.aa_sd / sqrt(juv.dens.aa_n),
    upper = juv.dens.aa_mu + 2 * juv.dens.aa_sd / sqrt(juv.dens.aa_n)
  ) |>
  mutate(A_SECTOR = factor(A_SECTOR, levels = c(
    "CG", "PC", "CL", "CA", "IN", "TO", "WH", "PO", "SW", "CB"
  ))) |>
  ggplot(aes(y = juv.dens.aa_mu, colour = Site, x = as.numeric(as.character(cREPORT_YEAR)))) +
  geom_line(show.legend = FALSE) +
  geom_pointrange(aes(ymin = lower, ymax = upper), show.legend = FALSE) +
  facet_grid(A_SECTOR ~ SHELF, scales = "free") +
  scale_y_continuous("Juvenile density") +
  scale_x_continuous("")
## ----end

## Conclusions:
## - substantial combinations with all 0's
## - this could be highly problematic

## Explore which Sector/Shelf/Year combinations have only 0 counts
## ---- q3_eda_2
data |>
  group_by(A_SECTOR, SHELF, cREPORT_YEAR) |>
  summarise(Abund = sum(abund)) |>
  filter(Abund == 0) |>
  as.data.frame()
## ----end

## 61 Sector/Shelf/Years with all 0's It is likely we will need to
## exclude these from the analyses and put them back in at the end -
## they will destabalise the model without very specific and strong
## priors

## these are likely to be very problematic

## Lets explore distributions over the space
## ---- q3_eda_3
data |>
  group_by(A_SECTOR, SHELF) |>
  ggplot(aes(x = abund)) +
  geom_histogram() +
  facet_grid(A_SECTOR ~ SHELF)
## ----end

## Could be zero-inflated - more so in some sector/shelf
## ---- q3_eda_4
data |>
  group_by(cREPORT_YEAR) |>
  ggplot(aes(x = abund)) +
  geom_histogram() +
  facet_wrap(~cREPORT_YEAR)
## ----end

## More zero-inflation


## Relative abundance

## ---- q3_eda_5
data |>
  group_by(A_SECTOR, SHELF, Site, cREPORT_YEAR) |>
  summarise(across(propn.juv, list(mu = mean, med = median,
    sd = sd, n = length))) |>
  mutate(
    lower = propn.juv_mu - 2 * propn.juv_sd / sqrt(propn.juv_n),
    upper = propn.juv_mu + 2 * propn.juv_sd / sqrt(propn.juv_n)
  ) |> 
  mutate(A_SECTOR = factor(A_SECTOR, levels = c(
    "CG", "PC", "CL", "CA", "IN", "TO", "WH", "PO", "SW", "CB"
  ))) |>
  ggplot(aes(y = propn.juv_mu, colour = Site, x = as.numeric(as.character(cREPORT_YEAR)))) +
  geom_line(show.legend = FALSE) +
  geom_pointrange(aes(ymin = lower, ymax = upper), show.legend = FALSE) +
  facet_grid(A_SECTOR ~ SHELF, scales = "free")
## ----end

## Conclusions:
## - substantial combinations with all 0's
## - this could be highly problematic

## Explore which Sector/Shelf/Year combinations have only 0 counts
## ---- q3_eda_2
data |>
  group_by(A_SECTOR, SHELF, cREPORT_YEAR) |>
  summarise(Propn.Juv = sum(propn.juv)) |>
  filter(Propn.Juv == 0) |>
  as.data.frame()
## ----end

## 57 Sector/Shelf/Years with all 0's It is likely we will need to
## exclude these from the analyses and put them back in at the end -
## they will destabalise the model without very specific and strong
## priors

## these are likely to be very problematic

## Lets explore distributions over the space
## ---- q3_eda_3
data |>
  group_by(A_SECTOR, SHELF) |>
  ggplot(aes(x = propn.juv)) +
  geom_histogram() +
  facet_grid(A_SECTOR ~ SHELF)
## ----end

## Could be zero-inflated - more so in some sector/shelf
## ---- q3_eda_4
data |>
  group_by(cREPORT_YEAR) |>
  ggplot(aes(x = propn.juv)) +
  geom_histogram() +
  facet_wrap(~cREPORT_YEAR)
## ----end

## More zero-inflation
