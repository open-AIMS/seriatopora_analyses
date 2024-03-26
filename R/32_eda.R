## Retrieve the data
## ---- q2_retrieve_data
load(file = "../data/processed/data_q2.R")
## ----end

## This is similar to above (for lookup data), except that was
## performed on a reef level, whereas the current is performed on the
## transect levels

## tests <- data |>
##   mutate(COVER = n.points / total.points) |> 
##   nest(.by = c("FULLREEF_ID", "REEF", "Dist.number")) |>
##     mutate(BeforeAfter = map(
##       .x = data,
##       .f = ~ before_after_tests(.x)
##     )) 


## tests |>
##   dplyr::select(-data) |>
##   unnest(c(BeforeAfter)) |>
##   filter(before_flag == "no before") |>
##   mutate(across(where(is.numeric), round, 3)) |>
##   as.data.frame() |>
##   filter(REEF_NAME == "NORTH DIRECTION REEF", TRANSECT_NO == 1, Dist.number == 1)


## data |> filter(Transect == "North Direction Island 1 1", Dist.number == 1) |> as.data.frame() |> head()
## data |> filter(Transect == "North Direction Island 1 1") |> as.data.frame() |> head()
## data |> filter(AIMS_REEF_NAME == "North Direction Island", SITE_NO == 1, TRANSECT_NO == 1) |> as.data.frame() |> head()

## lookup |> filter(REEF == "North Direction Island") |> as.data.frame() |> head()
## ## the issue is that for this reef, although the lookup has before
## ## listed for 1994, the benthic data only starts in 1995.




## ---- eda_basic
data |> dim()
data |>
        group_by(Dist.number) |>
        count()
## ----end


## How many combinations have only 0's
data |>
  group_by(Dist.time, A_SECTOR, SHELF, Site, Transect) |>
  summarise(Points = sum(n.points)) |>
  filter(Points == 0)

## Raw means to get a sense of the data
data |>
  group_by(Dist.time) |>
  summarise(cover = mean(n.points / total.points),
    cover_sd = sd(n.points / total.points),
    lower = cover - cover_sd,
    upper = cover + cover_sd) |>
  ggplot(aes(y = cover, x = Dist.time, colour = Dist.time)) +
  geom_pointrange(aes(ymin = lower, ymax = upper), position = position_dodge(width = 0.4))

data |>
  group_by(Dist.time, A_SECTOR, SHELF) |>
  summarise(cover = mean(n.points / total.points),
    cover_sd = sd(n.points / total.points),
    lower = cover - cover_sd,
    upper = cover + cover_sd) |>
  ggplot(aes(y = cover, x = A_SECTOR, colour = Dist.time)) +
  geom_pointrange(aes(ymin = lower, ymax = upper), position = position_dodge(width = 0.4)) +
  facet_grid(~SHELF)

data |>
  group_by(Dist.time, A_SECTOR, SHELF) |>
  summarise(cover = mean(n.points / total.points),
    cover_sd = sd(n.points / total.points),
    lower = cover - cover_sd,
    upper = cover + cover_sd) |>
  ggplot(aes(y = cover, x = Dist.time, colour = Dist.time)) +
  geom_hline(yintercept = 0) +
  geom_pointrange(aes(ymin = lower, ymax = upper),
    position = position_dodge(width = 0.4)) +
  facet_grid(A_SECTOR ~ SHELF, scales = "free")

## ---- eda_before_after_sector_shelf
data |>
  group_by(Dist.time, A_SECTOR, SHELF) |>
  summarise(cover = mean(n.points / total.points),
    cover_sd = sd(n.points / total.points),
    lower = cover - cover_sd,
    upper = cover + cover_sd) |>
  ggplot(aes(y = cover, x = Dist.time, colour = Dist.time)) +
  geom_hline(yintercept = 0) +
  geom_pointrange(aes(ymin = lower, ymax = upper),
    position = position_dodge(width = 0.4)) +
  facet_wrap(A_SECTOR ~ SHELF,
    scales = "free",
    labeller = label_wrap_gen(multi_line = FALSE)
  )
## ----end


data |>
  mutate(Dist.time2 = ifelse(Dist == "Before", "Before", "After"),
    Dist.time2 = forcats::fct_relevel(Dist.time2, "Before")) |>
  ## pull(Dist.time) |> levels()
 filter(Dist.time != Dist.time2) |> as.data.frame() |> head()


## ---- eda_before_after_disturbances
data |>
  mutate(Dist.time2 = ifelse(Dist == "Before", "Before", "After")) |> 
  group_by(A_SECTOR, SHELF, Dist, Dist.time) |>
  reframe(
    Mean = c(mean(n.points / total.points)), 
    SD = c(sd(n.points / total.points)), 
    N = c(n(), n())
  ) |>
  mutate(
    lower = Mean - 2 * (SD / sqrt(N)),
    upper = Mean + 2 * (SD / sqrt(N))
  ) |>
  ggplot(aes(y = Mean, x = Dist, colour = Dist.time)) +
  geom_pointrange(aes(ymin = lower, ymax = upper)) +
  facet_wrap(A_SECTOR ~ SHELF,
    scales = "free",
    labeller = label_wrap_gen(multi_line = FALSE)
  )
## ----end



data |>
  ggplot(aes(x = n.points)) +
  geom_histogram()

## ---- eda_zero_inflation
data |>
  ggplot(aes(x = n.points)) +
  geom_histogram() +
  facet_grid(SecShelf ~ Dist.time)
## ----end

