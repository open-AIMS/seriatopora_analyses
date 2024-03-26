
## load the data
## ---- load_data_3
load(file = "../data/primary/data.RData")
load(file = "../data/primary/lookup.RData")
## ----end

old_lookup <- lookup

## join the lookup to the data
## ---- common_fields
data |> colnames()
lookup |> colnames()

## What fields do they have in common
colnames(data)[colnames(data) %in% colnames(lookup)]
## ----end


## But some fields have changed names!!
## ---- join_lookup
data |>
  left_join(lookup,
    by = c("A_SECTOR", "FULLREEF_ID", "REPORT_YEAR" = "year")
  )
## ----end


## Pay attention to that error message

## ---- explore_join_warning_1
data |>
  slice(4716) |>
  dplyr::select(P_CODE, A_SECTOR, SHELF, FULLREEF_ID, REPORT_YEAR)
## ----end


## Ah - multiple disturbances in a single year
## ---- explore_join_warning_2
lookup |>
  filter(A_SECTOR == "CA", FULLREEF_ID == "16064S", year == 2012) |>
  as.data.frame()
## ----end


## lets collapse these together and dummy code the disturbance types
lookup |> dim()

## We just want Before and After, and each needs to have the same DISTURBANCE_TYPE

## 1. Remove cases where:
##    - DISTURBANCE_TYPE is NA
##    - Dist.number is 0 (sequence without a known Before)
##    - Dist.type is "Pre"
##    - Dist.time is "Recovery"
## 2. Remove sequences (Dist.number) within REEF that do not have a
##    Before and an After
## 3. Create a DIST_TYPE that is the collapsed concatenation of the
##    non "n" and non NA DISTURBANCE_TYPE values within each sequence
## 4. Use pivot wider to create dummy codes for each disturbance type
## lookup1 <- lookup
## lookup <- lookup1

## ---- explore_join_warning_3
lookup |>
        filter(REEF == "Pompey Reef No.1") |>
        as.data.frame()
## ----end

## - Dist.number 1, does not have a Before
## - Dist.time == "Before", Dist.number == 2 has DISTURBANCE_TYPE of "n"
## - Before (n), During (s), After (u)
## ---- explore_join_warning_4
lookup |>
        filter(REEF == "Arlington Reef") |>
        as.data.frame()
## ----end

## - Dist.time == "Before" for Dist.number == 1 is a "n"
## - Dist.time == "After" for Dist.number == 1 has two rows
## - Dist.number == 2, Before, During and After all have different DISTURBANCE_TYPE (n, d, b)

## lookup <- old_lookup 
## ---- before_after_tests
tests <- lookup |>
  filter(!is.na(DISTURBANCE_TYPE),
    Dist.number > 0,
    Dist.time != "Recovery",
    Dist.type != "Pre") |>
  nest(.by = c("FULLREEF_ID", "REEF", "Dist.number")) |>
    mutate(BeforeAfter = map(
      .x = data,
      .f = ~ before_after_tests(.x)
    )) 
save(tests, file = "../data/processed/q2_tests_1.RData")
## ----end
## .x <- tests[14, "data"][[1]][[1]]

## Instances where there are no before cases
## ---- before_after_tests_no_before
load(file = "../data/processed/q2_tests_1.RData")
tests |>
  dplyr::select(-data) |>
  unnest(c(BeforeAfter)) |>
  filter(before_flag == "no before") |>
  mutate(across(where(is.numeric), round, 3)) |> 
  as.data.frame() |>
  datatable(class = "compact") 
  ## formatRound(columns = 1:5, digits = 3)
  ## formatStyle(columns = "REEF", fontSize = "8pt")
## ----END

## Instances where there are multiple before cases
## ---- before_after_tests_multiple_before
load(file = "../data/processed/q2_tests_1.RData")
tests |>
  dplyr::select(-data) |>
  unnest(c(BeforeAfter)) |>
  filter(before_flag == "multiple before") |>
  mutate(across(where(is.numeric), round, 3)) |> 
  as.data.frame() |> 
  datatable(class = "compact")  
## ----end

## Instances where there After Disturbance Type is "n"
## ---- before_after_tests_n_after
load(file = "../data/processed/q2_tests_1.RData")
tests |>
  dplyr::select(-data) |>
  unnest(c(BeforeAfter)) |>
  filter(n_flag == "after disturb n") |>
  mutate(across(where(is.numeric), round, 3)) |> 
  as.data.frame() |> 
  datatable(class = "compact")
## ----end

## Instances where there are no after cases
## ---- before_after_tests_no_after
load(file = "../data/processed/q2_tests_1.RData")
tests |>
  dplyr::select(-data) |>
  unnest(c(BeforeAfter)) |>
  filter(after_flag == "no after") |>
  mutate(across(where(is.numeric), round, 3)) |> 
  as.data.frame() |> 
  datatable(class = "compact")
## ----end

## Instances where cover was greater after than before
## ---- before_after_tests_cover_greater_after
load(file = "../data/processed/q2_tests_1.RData")
tests |>
  dplyr::select(-data) |>
  unnest(c(BeforeAfter)) |>
  filter(cover_flag == "cover increased") |>
  mutate(across(where(is.numeric), round, 3)) |> 
  as.data.frame() |> 
  datatable(class = "compact")
## ----end

## ---- make_excludes
load(file = "../data/processed/q2_tests_1.RData")
excludes <- tests |>
  dplyr::select(-data) |>
  unnest(c(BeforeAfter)) |>
  filter(!is.na(before_flag) |
           !is.na(after_flag) |
            !is.na(n_flag) |
            !is.na(cover_flag)) |>
  dplyr::select(FULLREEF_ID, REEF, Dist.number) |>
  distinct() |>
  mutate(EXCLUDE = TRUE)
save(excludes, file = "../data/processed/q2_excludes_1.RData")
## ----end



## lookup <- old_lookup 
## ---- process_lookup
load(file = "../data/processed/q2_excludes_1.RData")
lookup <- 
  lookup |>
  ## remove cases without disturbances
  filter(
    !is.na(DISTURBANCE_TYPE),
    Dist.number > 0,
    Dist.time != "Recovery",
    Dist.type != "Pre"
  ) |>
  left_join(excludes) |>
  filter(is.na(EXCLUDE)) |> 
  pivot_wider(id_cols = everything(),
    values_fn = \(x) sum(!is.na(x)),
    values_fill = 0,
    names_from = DISTURBANCE_TYPE, values_from = DISTURBANCE_TYPE) |>
  dplyr::select(-EXCLUDE)
save(lookup, file = "../data/processed/q2_lookup_1.R") 
## ----end

## old
if (1 == 2) {
  lookup <-
    lookup |>
    ## remove cases without disturbances
    filter(
      !is.na(DISTURBANCE_TYPE),
      Dist.number > 0,
      Dist.time != "Recovery",
      Dist.type != "Pre"
    ) |>
    ## within each REEF and Dist.number, assign a DISTURBANCE_TYPE for Before to be the same as After
    ## collapse disturbances
    group_by(FULLREEF_ID, REEF, Dist.number) |>
    filter(str_detect(Dist.time, "Before|After")) |>
    filter(sum(str_detect(unique(Dist.time), "Before|After")) == 2) |>
    ungroup() |>
    nest(.by = c(FULLREEF_ID, REEF, Dist.number)) |>
    ## slice(1:20) |>
    # filter(REEF == "Agincourt Reef No.1") |>
    mutate(data1 = purrr::map(
      .x = data,
      .f = ~ {
        print(.x)
        Before <- .x |>
          filter(Dist.time == "Before") |>
          dplyr::select(-DISTURBANCE_TYPE)
        After <- .x |>
          filter(Dist.time == "After") |>
          dplyr::select(-DISTURBANCE_TYPE)
        Dists <- .x |>
          pull(DISTURBANCE_TYPE) |>
          unique()
        Dists <- Dists[Dists != "n"]
        Before <- Before |> reframe(Before, DISTURBANCE_TYPE = Dists)
        After <- After |> reframe(After, DISTURBANCE_TYPE = Dists)
        bind_rows(Before, After)
        # print("hi")
        # print(as.data.frame(a))
        ## asdf
      }
    )) |>
    dplyr::select(-data) |>
    unnest(c(data1)) |>
    pivot_wider(
      id_cols = everything(),
      values_fn = \(x) sum(!is.na(x)),
      values_fill = 0,
      names_from = DISTURBANCE_TYPE, values_from = DISTURBANCE_TYPE
    )
}

## .x <- lookup[36, "data"][[1]][[1]]

## Try the join again
## ---- join_lookup_2
data <-
  data |>
  left_join(lookup,
    by = c("A_SECTOR", "FULLREEF_ID", "REPORT_YEAR" = "year")
  )
save(data, file = "../data/processed/q2_data_1.R") 
## ----end


## Final prep is to add some derived variables
## ---- add_derived
load(file = "../data/processed/q2_data_1.R") 
data <-
  data |>
  mutate(
    P_CODE2 = ifelse(P_CODE == "IN", "IN", "LTMP"),
    Reef = paste(P_CODE2, AIMS_REEF_NAME),
    Site = paste(Reef, DEPTH, SITE_NO),
    Transect = paste(Site, TRANSECT_NO)
  ) |>
  filter(Dist.time %in% c("Before", "After")) |>
  filter(!is.na(Dist.number)) |>
  droplevels() |> 
  mutate(Dist.time = factor(Dist.time, levels = c("Before", "After"))) |>
  mutate(SecShelf = factor(paste(A_SECTOR, SHELF)))

save(data, file = "../data/processed/q2_data_2.R") 
## ----end




## Having now joined the disturbance data to the benthic data, we need
## to perform the tests again.

## ---- before_after_tests_2
load(file = "../data/processed/q2_data_2.R") 
tests <- data |>
  mutate(COVER = n.points / total.points) |> 
  nest(.by = c("Reef", "FULLREEF_ID", "REEF", "Site", "Transect", "Dist.number")) |>
    mutate(BeforeAfter = map(
      .x = data,
      .f = ~ before_after_tests(.x)
    )) 
save(tests, file = "../data/processed/q2_tests_2.RData")
## ----end

## Instances where there are no before cases
## ---- before_after_tests_no_before_2
load(file = "../data/processed/q2_tests_2.RData")
tests |>
  dplyr::select(-data) |>
  unnest(c(BeforeAfter)) |>
  filter(before_flag == "no before") |>
  mutate(across(where(is.numeric), round, 3)) |> 
  as.data.frame() |>
  datatable(class = "compact") 
## ----END

## Instances where there are multiple before cases
## ---- before_after_tests_multiple_before_2
load(file = "../data/processed/q2_tests_2.RData")
tests |>
  dplyr::select(-data) |>
  unnest(c(BeforeAfter)) |>
  filter(before_flag == "multiple before") |>
  mutate(across(where(is.numeric), round, 3)) |> 
  as.data.frame() |> 
  datatable(class = "compact")  
## ----end

## Instances where there After Disturbance Type is "n"
## ---- before_after_tests_n_after_2
load(file = "../data/processed/q2_tests_2.RData")
tests |>
  dplyr::select(-data) |>
  unnest(c(BeforeAfter)) |>
  filter(n_flag == "after disturb n") |>
  mutate(across(where(is.numeric), round, 3)) |> 
  as.data.frame() |> 
  datatable(class = "compact")
## ----end

## Instances where there are no after cases
## ---- before_after_tests_no_after_2
load(file = "../data/processed/q2_tests_2.RData")
tests |>
  dplyr::select(-data) |>
  unnest(c(BeforeAfter)) |>
  filter(after_flag == "no after") |>
  mutate(across(where(is.numeric), round, 3)) |> 
  as.data.frame() |> 
  datatable(class = "compact")
## ----end

## Instances where cover was greater after than before
## ---- before_after_tests_cover_greater_after_2
load(file = "../data/processed/q2_tests_2.RData")
tests |>
  dplyr::select(-data) |>
  unnest(c(BeforeAfter)) |>
  filter(cover_flag == "cover increased") |>
  mutate(across(where(is.numeric), round, 3)) |> 
  as.data.frame() |> 
  datatable(class = "compact")
## ----end


## ---- make_excludes_2
load(file = "../data/processed/q2_tests_2.RData")
excludes <- tests |>
  dplyr::select(-data) |>
  unnest(c(BeforeAfter)) |>
  filter(!is.na(before_flag) |
           !is.na(after_flag) |
            !is.na(n_flag)) |>
  dplyr::select(Transect, Dist.number) |>
  distinct() |>
  mutate(EXCLUDE = TRUE)
save(excludes, file = "../data/processed/q2_excludes_2.RData")
## ----end

## ---- exclude_data
data <-
  data |>
  left_join(excludes) |>
  filter(is.na(EXCLUDE))
## ----end

## ---- encode_disturbances
load(file = "../data/processed/q2_data_2.R") 
data <- data |>
  mutate(across(c(s, c, d, b, u), \(x) ifelse(Dist.time == "Before", 0, x))) |> 
  mutate(Dist = case_when(
    s > 0 & (c + d + b + u == 0) ~ "s",
    s > 0 & (c + d + b + u > 0) ~ "sm",
    c > 0 & (s + d + b + u == 0) ~ "c",
    c > 0 & (s + d + b + u > 0) ~ "cm",
    d > 0 & (s + c + b + u == 0) ~ "d",
    d > 0 & (s + c + b + u > 0) ~ "dm",
    b > 0 & (s + c + d + u == 0) ~ "b",
    b > 0 & (s + c + d + u > 0) ~ "bm",
    u > 0 & (s + c + d + b == 0) ~ "u",
    u > 0 & (s + c + d + b > 0) ~ "um",
    .default = "Before"
  )) |> 
  mutate(Dist = forcats::fct_relevel(Dist, "Before")) |>
  mutate(Dist.type = ifelse(s + c + d + b + u > 1, "Cumulative", "Single")) |>
  mutate(across(c(s, c, d, b, u), \(x) ifelse(x > 1, 1, x))) |>
  mutate(across(c(s, c, d, b, u), as.factor)) |>
  mutate(event = paste(Transect, Dist.number)) |>
  mutate(SSDist = paste(SecShelf, Dist))
save(data, file = "../data/processed/data_q2.R") 
## ----end
