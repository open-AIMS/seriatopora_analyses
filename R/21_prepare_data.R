## load the data
## ---- q1_load_data
load(file = "../data/primary/data.RData")
## ----end

## ---- q1_process_data
data <-
  data |>
  mutate(
    Site = paste(AIMS_REEF_NAME, SITE_NO),
    Transect = paste(Site, TRANSECT_NO),
    SecShelf = factor(paste(A_SECTOR, SHELF)),
    zone_depth = factor(paste(REEF_ZONE, DEPTH)),
    cREPORT_YEAR = factor(REPORT_YEAR)
  )
## ----end

## Save the data
## ---- q1_save_data
save(data, file = "../data/processed/data_q1.RData")
## ----end
