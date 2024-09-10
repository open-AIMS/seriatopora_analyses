
## Read in the primary data
## ---- read_data_1
## data <- get(load("../data/primary/seriatopora.RData"))
data1 <- get(load("../data/primary/seriatopora1.RData"))
## ----end

## ---- make_data_consistent_with_previous_version
data <-
  data1 |>
  dplyr::select(-LATITUDE, -LONGITUDE, -cover) |>
  group_by(
    P_CODE, A_SECTOR, SHELF, REEF, REEF_ZONE, DEPTH, REPORT_YEAR,
    SITE_NO, TRANSECT_NO, Decade, GROUP_CODE, COMP_2021,
    COMP_2021_DESCRIPTION
  ) |>
  summarise(n.points = sum(n.points), total.points = sum(total.points)) |>
  ungroup() |>
  dplyr::rename(AIMS_REEF_NAME = REEF) 
## ----end


## ---- glimpse_data_1
data |> glimpse()
## ----end


## We restrict to just the Seriatopora (G_SER)
## ---- ser_only
data <- data |>
        filter(COMP_2021 == "G_SER") |>
        droplevels()
## ----end

## ---- save_data_1
save(data, file = "../data/primary/data.RData")
## ----end


## read in the disturbance lookup
## ---- read_lookup
lookup <- read_csv("../data/primary/dist.lookup 10.csv")
lookup <- lookup |> dplyr::select(-DEPTH) # until the most recent dist.lookup data, DEPTH was not included
## ----end

## ---- save_data_2
save(lookup, file = "../data/primary/lookup.RData")
## ----end

## ---- read_lookup_mmp
lookup_mmp <- read_csv("../data/primary/dist_lookup_MMP.csv")
lookup_mmp <- lookup_mmp |>
  dplyr::rename(year = report_year) |>
  dplyr::rename(REEF = AIMS_REEF_NAME) |>
  mutate(REEF_ZONE = ifelse(is.na(REEF_ZONE), "", REEF_ZONE),
    REEF_ZONE = factor(REEF_ZONE))
## ----end

## ---- save_data_2_mmp
save(lookup_mmp, file = "../data/primary/lookup_mmp.RData")
## ----end
