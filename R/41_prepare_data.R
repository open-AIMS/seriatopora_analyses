## load the data
## ---- q3_load_data
## data <- get(load(file = "../data/primary/juv.ser.ltmp.mmp.RData"))
data <- get(load(file = "../data/primary/ser.juv.density.Rdata"))
## ----end

## ---- q3_process_data
data <-
  data |>
  mutate(
    Site = paste(AIMS_REEF_NAME, SITE_NO),
    SecShelf = factor(paste(A_SECTOR, SHELF)),
    zone_depth = factor(paste(REEF_ZONE, DEPTH)),
    cREPORT_YEAR = factor(REPORT_YEAR),
    abund = ABUNDANCE
  )
## ----end

## Save the data
## ---- q3_save_data
save(data, file = "../data/processed/data_q3.RData")
## ----end
