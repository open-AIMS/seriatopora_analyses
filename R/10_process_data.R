
## Read in the primary data
## ---- read_data_1
data <- get(load("../data/primary/seriatopora.RData"))
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
lookup <- read_csv("../data/primary/dist.lookup 5.csv")
## ----end

## ---- save_data_2
save(lookup, file = "../data/primary/lookup.RData")
## ----end

