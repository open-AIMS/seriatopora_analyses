## load the data
## ---- load_data_3
load(file = "../data/primary/data.RData")
load(file = "../data/primary/lookup.RData")
load(file = "../data/primary/lookup_mmp.RData")
## ----end


## The structural differences between the LTMP and MMP disturbance
## lookups and data are sufficiently different that the joining of
## benthic data and disturbance data will need to take place separately
## for each. For example, whilst MMP (`P_CODE == "IN"`) reefs have a
## REEF_ZONE that distinguishes between multiple reefs (e.g. `Pandora`,
## `Havannah Reef` and `Fitzroy Island Reef`), this is not the case for
## LTMP.

## ---- split_data
data_mmp <- data |>
  filter(P_CODE == "IN") |>
  droplevels()
data <- data |>
  filter(P_CODE != "IN") |>
  droplevels()
## ----end
 
## LTMP data =================================

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
    ## by = c("A_SECTOR", "FULLREEF_ID", "REPORT_YEAR" = "year")
    by = c("A_SECTOR", "AIMS_REEF_NAME" = "REEF", "REPORT_YEAR" = "year")
  )
## ----end


## Pay attention to that error message

## ---- explore_join_warning_1
data |>
  slice(46) |>
  dplyr::select(P_CODE, A_SECTOR, SHELF, AIMS_REEF_NAME, REPORT_YEAR)
## ----end


## Ah - multiple disturbances in a single year
## ---- explore_join_warning_2
lookup |>
  filter(A_SECTOR == "CA", REEF == "Arlington Reef", year == 2012) |>
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


## for each reef/Dist.number, if Dist.type == "Cumulative", look through the Dist.time == "During" and get unique DISTURBANCE_TYPE values

## ---- correct_after_cumulative
replace_disturbance <- function(DISTURBANCE_TYPE, Dist.time) {
  val <- paste(DISTURBANCE_TYPE[Dist.time %in% c("During") & DISTURBANCE_TYPE != "n"], collapse = "")
  val <- ifelse(str_length(val) == 1,
    val,
    str_replace_all(val, "n", "")
  )
  val <- ifelse(val == "", "n", val)
  val
}
lookup <- lookup |>
  nest_by(REEF, Dist.number) |>
  mutate(data = list({
    data |>
      ## ## Handle the multiple (cumulative) disturbance cases
      ## mutate(DISTURBANCE_TYPE = ifelse(DISTURBANCE_TYPE == "n" & Dist.time == "After",
      ##   paste0(unique(DISTURBANCE_TYPE[Dist.type == "Cumulative" | Dist.time == "During"])),
      ##   DISTURBANCE_TYPE
      ## )) |>
      ## Handle the single disturbance cases
      mutate(DISTURBANCE_TYPE = ifelse((is.na(DISTURBANCE_TYPE) | DISTURBANCE_TYPE == "n") & Dist.time == "After",
        replace_disturbance(DISTURBANCE_TYPE, Dist.time),
        DISTURBANCE_TYPE
      ))
  })) |>
  unnest(data) |>
  arrange(FULLREEF_ID, year) |>
  ungroup()
save(lookup, file = "../data/processed/q2_lookup_0.RData")
## ----end


## - Dist.time == "Before" for Dist.number == 1 is a "n"
## - Dist.time == "After" for Dist.number == 1 has two rows
## - Dist.number == 2, Before, During and After all have different DISTURBANCE_TYPE (n, d, b)

## lookup <- old_lookup
## ---- before_after_tests
tests <- lookup |>
  filter(
    !is.na(DISTURBANCE_TYPE),
    Dist.number > 0,
    Dist.time != "Recovery",
    Dist.type != "Pre"
  ) |>
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
  datatable(class = c("compact", "white-space: nowrap")) 
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
  datatable(class = c("compact", "white-space: nowrap"))
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
  datatable(class = c("compact", "white-space: nowrap"))
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
  datatable(class = c("compact", "white-space: nowrap"))
## ----end

## Instances where cover was greater after than before
## ---- before_after_tests_cover_greater_after
load(file = "../data/processed/q2_tests_1.RData")
tests |>
  dplyr::select(-data) |>
  unnest(c(BeforeAfter)) |>
  filter(Dist.time %in% c("After")) |> 
  filter(cover_flag == "cover increased") |>
  mutate(across(where(is.numeric), round, 3)) |>
  as.data.frame() |>
  datatable(class = c("compact", "white-space: nowrap"))
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

lookup |> filter(str_length(DISTURBANCE_TYPE)>1) |> as.data.frame() |> head()

## lookup <- old_lookup
## ---- process_lookup
load(file = "../data/processed/q2_excludes_1.RData")
load(file = "../data/processed/q2_lookup_0.RData")
lookup <-
  lookup |>
  ## Duplicate rows with multiple disturbances
  mutate(DISTURBANCE_TYPE = str_split(DISTURBANCE_TYPE, "")) |> 
  unnest(DISTURBANCE_TYPE) |>                                    # Separate each character into a new row
  ## remove cases without disturbances
  filter(
    !is.na(DISTURBANCE_TYPE),
    Dist.number > 0,
    Dist.time != "Recovery",
    Dist.type != "Pre"
  ) |>
  left_join(excludes) |>
  filter(is.na(EXCLUDE)) |>
  pivot_wider(
    id_cols = everything(),
    values_fn = \(x) sum(!is.na(x)),
    values_fill = 0,
    names_from = DISTURBANCE_TYPE, values_from = DISTURBANCE_TYPE
  ) |>
  dplyr::select(-EXCLUDE)
save(lookup, file = "../data/processed/q2_lookup_1.RData")
## ----end

lookup |> filter(REEF == "Havannah Island", Dist.number == 1, year == 2004)|> as.data.frame() |> head()

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
    by = c("A_SECTOR", "AIMS_REEF_NAME" = "REEF", "REPORT_YEAR" = "year")
  )
save(data, file = "../data/processed/q2_data_1.RData")
## ----end

data |> filter(AIMS_REEF_NAME == "Havannah Island", Dist.number == 1, REPORT_YEAR == 2004)|> as.data.frame() |> head()

## Final prep is to add some derived variables
## ---- add_derived
load(file = "../data/processed/q2_data_1.RData")
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

save(data, file = "../data/processed/q2_data_2.RData")
## ----end

## ---- q2_data_view_2
load(file = "../data/processed/q2_data_2.RData")
data |>
  dplyr::relocate(A_SECTOR, SHELF, AIMS_REEF_NAME, REPORT_YEAR) |> 
  datatable(
    class = c("compact", "white-space: nowrap"),
    filter = "top",
    extensions = c("Buttons","FixedColumns", "KeyTable"),
    options = list(
      dom = "Bfrtip",
      scrollX = TRUE,
      fixedColumns = list(leftColumns = 5),
      keys = TRUE,       ## move via arrow keys
      autoWidth = TRUE,
      buttons = list(
        list(
          extend = "collection",
          buttons = c("copy", "csv", "excel", "pdf", "print"),
          text = "Download"
        )
      ),
      columnDefs = list(list(width = "10px", targets = 0))
    )
  )
## ----end


## Having now joined the disturbance data to the benthic data, we need
## to perform the tests again.

## ---- before_after_tests_2
load(file = "../data/processed/q2_data_2.RData")
tests <- data |>
  mutate(COVER = n.points / total.points) |>
  nest(.by = c("Reef", "FULLREEF_ID", "AIMS_REEF_NAME", "Site", "Transect", "Dist.number")) |>
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
  datatable(class = c("compact", "white-space: nowrap"))
## ----end


## Instances where there are multiple before cases
## ---- before_after_tests_multiple_before_2
load(file = "../data/processed/q2_tests_2.RData")
tests |>
  dplyr::select(-data) |>
  unnest(c(BeforeAfter)) |>
  filter(before_flag == "multiple before") |>
  mutate(across(where(is.numeric), round, 3)) |>
  as.data.frame() |>
  datatable(class = c("compact", "white-space: nowrap"))
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
  datatable(class = c("compact", "white-space: nowrap"))
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
  datatable(class = c("compact", "white-space: nowrap"))
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
load(file = "../data/processed/q2_data_2.RData")
if (nrow(excludes) > 0) {
  data <- data |>
    left_join(excludes) |>
    filter(is.na(EXCLUDE))
}
## ----end

## ---- encode_disturbances
data <- data |>
  mutate(across(c(s, c, d, b, u), \(x) ifelse(Dist.time == "Before", 0, x))) |>
  mutate(Dist = case_when(
    s == 1 & (c + d + b + u == 0) ~ "s",
    s > 1 & (c + d + b + u == 0) ~ "ss",
    s > 0 & (c + d + b + u > 0) ~ "sm",
    c == 1 & (s + d + b + u == 0) ~ "c",
    c > 1 & (s + d + b + u == 0) ~ "cc",
    c > 0 & (s + d + b + u > 0) ~ "cm",
    d == 1 & (s + c + b + u == 0) ~ "d",
    d > 1 & (s + c + b + u == 0) ~ "df",
    d > 0 & (s + c + b + u > 0) ~ "dm",
    b == 1 & (s + c + d + u == 0) ~ "b",
    b > 1 & (s + c + d + u == 0) ~ "bb",
    b > 0 & (s + c + d + u > 0) ~ "bm",
    u == 1 & (s + c + d + b == 0) ~ "u",
    u > 1 & (s + c + d + b == 0) ~ "uu",
    u > 0 & (s + c + d + b > 0) ~ "um",
    .default = "Before"
  )) |>
  mutate(Dist = forcats::fct_relevel(Dist, "Before")) |>
  mutate(Dist.type = ifelse(s + c + d + b + u > 1, "Cumulative", "Single")) |>
  mutate(across(c(s, c, d, b, u), \(x) ifelse(x > 1, 1, x))) |>
  mutate(across(c(s, c, d, b, u), as.factor)) |>
  mutate(event = paste(Transect, Dist.number)) |>
  mutate(SSDist = paste(SecShelf, Dist))
save(data, file = "../data/processed/data_q2.RData")
## ----end

## ---- q2_view_data
load(file = "../data/processed/data_q2.RData")
data |>
  dplyr::relocate(A_SECTOR, SHELF, AIMS_REEF_NAME, REPORT_YEAR) |> 
  datatable(
    class = c("compact", "white-space: nowrap"),
    filter = "top",
    extensions = c("Buttons","FixedColumns", "KeyTable"),
    options = list(
      dom = "Bfrtip",
      scrollX = TRUE,
      fixedColumns = list(leftColumns = 5),
      keys = TRUE,       ## move via arrow keys
      autoWidth = TRUE,
      buttons = list(
        list(
          extend = "collection",
          buttons = c("copy", "csv", "excel", "pdf", "print"),
          text = "Download"
        )
      ),
      columnDefs = list(list(width = "10px", targets = 0))
    )
  )
## ----end


## MMP data =================================

## join the lookup to the data
## ---- common_fields_mmp
data_mmp |> colnames()
lookup_mmp |> colnames()

## What fields do they have in common
colnames(data_mmp)[colnames(data_mmp) %in% colnames(lookup_mmp)]
## ----end

## But some fields have changed names!!
## ---- join_lookup_mmp
lookup_mmp <- lookup_mmp |> mutate(DEPTH = factor(DEPTH))
data_mmp |>
  left_join(lookup_mmp,
    ## by = c("A_SECTOR", "FULLREEF_ID", "REPORT_YEAR" = "year")
    by = c("A_SECTOR", "AIMS_REEF_NAME" = "REEF",
      "REPORT_YEAR" = "year", "DEPTH", "REEF_ZONE")
  ) 
## ----end

## Pay attention to that error message

## ---- explore_join_warning_1_mmp
data_mmp |>
  slice(2659) |>
  dplyr::select(
    P_CODE, A_SECTOR, SHELF, AIMS_REEF_NAME, DEPTH,
    REEF_ZONE, REPORT_YEAR
  )  
## ----end


## Ah - multiple disturbances in a single year
## ---- explore_join_warning_2_mmp
lookup_mmp |>
  filter(A_SECTOR == "IN", REEF == "Frankland Islands", REEF_ZONE == "East",
    year == 2017, DEPTH == 5) |>
  as.data.frame()  
## ----end

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

## ---- explore_join_warning_3_mmp
lookup_mmp |>
  filter(REEF == "Fitzroy Island", REEF_ZONE == "East", DEPTH == 5) |>
  as.data.frame()
## ----end

## ---- explore_join_warning_4_mmp
lookup_mmp |>
  filter(REEF == "High Island", REEF_ZONE == "East", DEPTH == 5) |>
  as.data.frame()
## ----end


## ---- correct_after_cumulative_mmp
replace_disturbance <- function(DISTURBANCE_TYPE, Dist.time) {
  val <- paste(DISTURBANCE_TYPE[Dist.time %in% c("During") & DISTURBANCE_TYPE != "n"], collapse = "")
  val <- ifelse(str_length(val) == 1,
    val,
    str_replace_all(val, "n", "")
  )
  val <- ifelse(val == "", "n", val)
  val
}
lookup_mmp <- lookup_mmp |>
  nest_by(REEF, REEF_ZONE, DEPTH, Dist.number) |>
  mutate(data = list({
    data |>
      ## ## Handle the multiple (cumulative) disturbance cases
      ## mutate(DISTURBANCE_TYPE = ifelse(DISTURBANCE_TYPE == "n" & Dist.time == "After",
      ##   paste0(unique(DISTURBANCE_TYPE[Dist.type == "Cumulative" | Dist.time == "During"])),
      ##   DISTURBANCE_TYPE
      ## )) |>
      ## Handle the single disturbance cases
      mutate(DISTURBANCE_TYPE = ifelse((is.na(DISTURBANCE_TYPE) | DISTURBANCE_TYPE == "n") & Dist.time == "After",
        replace_disturbance(DISTURBANCE_TYPE, Dist.time),
        DISTURBANCE_TYPE
      ))
  })) |>
  unnest(data) |>
  arrange(REEF, REEF_ZONE, DEPTH, year) |>
  ungroup()
save(lookup_mmp, file = "../data/processed/q2_lookup_mmp_0.RData")
## ----end

## - Dist.time == "Before" for Dist.number == 1 is a "n"
## - Dist.time == "After" for Dist.number == 1 has two rows
## - Dist.number == 2, Before, During and After all have different DISTURBANCE_TYPE (n, d, b)

## lookup <- old_lookup
## ---- before_after_tests_mmp
tests_mmp <- lookup_mmp |>
  filter(
    !is.na(DISTURBANCE_TYPE),
    Dist.number > 0,
    Dist.time != "Recovery",
    Dist.type != "Pre"
  ) |>
  nest(.by = c("REEF", "REEF_ZONE", "DEPTH", "Dist.number")) |>
  mutate(BeforeAfter = map(
    .x = data,
    .f = ~ before_after_tests(.x)
  ))
save(tests_mmp, file = "../data/processed/q2_tests_mmp_1.RData")
## ----end
## .x <- tests_mmp[14, "data"][[1]][[1]]


## Instances where there are no before cases
## ---- before_after_tests_no_before_mmp
load(file = "../data/processed/q2_tests_mmp_1.RData")
tests_mmp |>
  dplyr::select(-data) |>
  unnest(c(BeforeAfter)) |>
  filter(before_flag == "no before") |>
  mutate(across(where(is.numeric), round, 3)) |>
  as.data.frame() |>
  datatable(class = c("compact", "white-space: nowrap")) 
## ----end

## I looked into this for Barren Island, This one is triggered because
## the before `DISTURBANCE_TYPE` is `NA` - it needs to be something
## else.

## ---- before_after_tests_no_before_mmp_2
options(width = 250)
lookup_mmp |>
  filter(REEF == "Barren Island", DEPTH == 2, Dist.number == 1) |>
  as.data.frame() |>
  head()
options(width = 80)
## ----end

## I have looked into `Fitzroy Island` `East`. This one is triggered
## because the before case is `Dist.number` is 0 rather than 1. That
## is, the before and after need to be listed as the same
## `Dist.number`.

## ---- before_after_tests_no_before_mmp_3
options(width = 250)
lookup_mmp |>
  filter(REEF == "Fitzroy Island", REEF_ZONE == "East", DEPTH == 5) |>
  as.data.frame() |>
  head()
options(width = 80)
## ----end

## ---- before_after_tests_no_before_mmp_4
options(width = 250)
lookup_mmp |>
  filter(REEF == "Dunk Island", REEF_ZONE == "North", DEPTH == 2, Dist.number == 1) |>
  as.data.frame() |>
  head()
options(width = 80)
## ----end

## ---- before_after_tests_no_before_mmp_4
options(width = 250)
lookup_mmp |>
  filter(REEF == "Frankland Islands", REEF_ZONE == "East", DEPTH == 2, Dist.number == 1) |>
  as.data.frame() |>
  head()
options(width = 80)
## ----end

## ---- before_after_tests_no_before_mmp_4
options(width = 250)
lookup_mmp |>
  filter(REEF == "Frankland Islands", REEF_ZONE == "West", DEPTH == 5, Dist.number == 1) |>
  as.data.frame() |>
  head()
options(width = 80)
## ----end

## ---- before_after_tests_no_before_mmp_4
options(width = 250)
lookup_mmp |>
  filter(REEF == "Double Cone Island", DEPTH == 2, Dist.number %in% c(0, 1)) |>
  as.data.frame() |>
  head(10)
options(width = 80)
## ----end
## Instances where there are multiple before cases
## ---- before_after_tests_multiple_before_mmp
load(file = "../data/processed/q2_tests_mmp_1.RData")
tests_mmp |>
  dplyr::select(-data) |>
  unnest(c(BeforeAfter)) |>
  filter(before_flag == "multiple before") |>
  mutate(across(where(is.numeric), round, 3)) |>
  as.data.frame() |>
  datatable(class = c("compact", "white-space: nowrap"))
## ----end

## I have looked into `Double Cone Island` `Dist.number` 2. This one
## is triggered because for `Dist.number` 2, there seems to be two
## sequences of before and after (one starting in 2016, the other in
## 2019).

## The same seems to be the case for Snapper Island North

## ---- before_after_tests_multiple_before_mmp_2
options(width = 250)
lookup_mmp |>
  filter(REEF == "Snapper Island", REEF_ZONE == "South",
    DEPTH == 5, Dist.number == 1) |>
  as.data.frame() |>
  head()
options(width = 80)
## ----end


## Instances where there After Disturbance Type is "n"
## ---- before_after_tests_n_after_mmp
load(file = "../data/processed/q2_tests_mmp_1.RData")
tests_mmp |>
  dplyr::select(-data) |>
  unnest(c(BeforeAfter)) |>
  filter(n_flag == "after disturb n") |>
  mutate(across(where(is.numeric), round, 3)) |>
  as.data.frame() |>
  datatable(class = c("compact", "white-space: nowrap"))
## ----end

## Instances where there are no after cases
## ---- before_after_tests_no_after_mmp
load(file = "../data/processed/q2_tests_mmp_1.RData")
tests_mmp |>
  dplyr::select(-data) |>
  unnest(c(BeforeAfter)) |>
  filter(after_flag == "no after") |>
  mutate(across(where(is.numeric), round, 3)) |>
  as.data.frame() |>
  datatable(class = c("compact", "white-space: nowrap"))
## ----end


## ---- before_after_tests_multiple_before_mmp_4
options(width = 250)
lookup_mmp |>
  filter(REEF == "Fitzroy Island", REEF_ZONE == "West",
    DEPTH == 2, Dist.number == 2) |>
  as.data.frame() |>
  head()
options(width = 80)
## ----end


## Instances where cover was greater after than before
## ---- before_after_tests_cover_greater_after_mmp
load(file = "../data/processed/q2_tests_mmp_1.RData")
tests_mmp |>
  dplyr::select(-data) |>
  unnest(c(BeforeAfter)) |>
  filter(Dist.time %in% c("After")) |> 
  filter(cover_flag == "cover increased") |>
  mutate(across(where(is.numeric), round, 3)) |>
  as.data.frame() |>
  datatable(class = c("compact", "white-space: nowrap"))
## ----end


## ---- make_excludes_mmp
load(file = "../data/processed/q2_tests_mmp_1.RData")
excludes_mmp <- tests_mmp |>
  dplyr::select(-data) |>
  unnest(c(BeforeAfter)) |>
  filter(!is.na(before_flag) |
    !is.na(after_flag) |
    !is.na(n_flag) |
    !is.na(cover_flag)) |>
  dplyr::select(REEF, REEF_ZONE, DEPTH, Dist.number) |>
  distinct() |>
  mutate(EXCLUDE = TRUE)
save(excludes_mmp, file = "../data/processed/q2_excludes_mmp_1.RData")
## ----end

## ---- process_lookup_mmp
load(file = "../data/processed/q2_excludes_mmp_1.RData")
load(file = "../data/processed/q2_lookup_mmp_0.RData")
lookup_mmp <-
  lookup_mmp |>
  ## Duplicate rows with multiple disturbances
  mutate(DISTURBANCE_TYPE = str_split(DISTURBANCE_TYPE, "")) |> 
  unnest(DISTURBANCE_TYPE) |>                                    # Separate each character into a new row
  ## remove cases without disturbances
  filter(
    !is.na(DISTURBANCE_TYPE),
    Dist.number > 0,
    Dist.time != "Recovery",
    Dist.type != "Pre"
  ) |>
  left_join(excludes_mmp) |>
  filter(is.na(EXCLUDE)) |>
  pivot_wider(
    id_cols = everything(),
    values_fn = \(x) sum(!is.na(x)),
    values_fill = 0,
    names_from = DISTURBANCE_TYPE, values_from = DISTURBANCE_TYPE
  ) |>
  dplyr::select(-EXCLUDE)
save(lookup_mmp, file = "../data/processed/q2_lookup_mmp_1.RData")
## ----end

## Try the join again
## ---- join_lookup_2_mmp
data_mmp <-
  data_mmp |>
  left_join(lookup_mmp,
    by = c("A_SECTOR", "AIMS_REEF_NAME" = "REEF",
      "REPORT_YEAR" = "year", "DEPTH", "REEF_ZONE")
  )
save(data_mmp, file = "../data/processed/q2_data_mmp_1.RData")
## ----end


data_mmp |> slice(8219:8220) |> as.data.frame()
lookup_mmp |> filter(REEF == "Barren Island", DEPTH == 2, Dist.number == 3, year == 2010) |> as.data.frame()

lookup_mmp |> slice(103) |> as.data.frame()
data_mmp |>
  filter(AIMS_REEF_NAME == "Fitzroy Island", REEF_ZONE == "East", DEPTH == 2, REPORT_YEAR == 2008) |>
  as.data.frame()

data_mmp |>
  filter(b > 0, c > 0) |>
  as.data.frame() |> tail()
## # A tibble: 30 × 30
##    P_CODE A_SECTOR SHELF AIMS_REEF_NAME    REEF_ZONE DEPTH REPORT_YEAR SITE_NO
##    <fct>  <chr>    <fct> <chr>             <fct>     <fct>       <dbl> <fct>
##  1 IN     IN       I     Frankland Islands East      5            2017 1
##  2 IN     IN       I     Frankland Islands East      5            2017 1
##  3 IN     IN       I     Frankland Islands East      5            2017 1
##  4 IN     IN       I     Frankland Islands East      5            2017 1
##  5 IN     IN       I     Frankland Islands East      5            2017 1
##  6 IN     IN       I     Frankland Islands East      5            2017 2
##  7 IN     IN       I     Frankland Islands East      5            2017 2
##  8 IN     IN       I     Frankland Islands East      5            2017 2
##  9 IN     IN       I     Frankland Islands East      5            2017 2
## 10 IN     IN       I     Frankland Islands East      5            2017 2

lookup_mmp |>
  filter(REEF == "Frankland Islands", REEF_ZONE == "East", DEPTH == 5, Dist.number == 3) |>
  as.data.frame()
lookup_mmp |>
  filter(REEF == "High Island", REEF_ZONE == "East", DEPTH == 5, Dist.number == 2) |>
  as.data.frame()

## Final prep is to add some derived variables
## ---- add_derived_mmp
load(file = "../data/processed/q2_data_mmp_1.RData")
data_mmp <-
  data_mmp |>
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

save(data_mmp, file = "../data/processed/q2_data_mmp_2.RData")
## ----end

## Having now joined the disturbance data to the benthic data, we need
## to perform the tests again.

## ---- before_after_tests_2_mmp
load(file = "../data/processed/q2_data_mmp_2.RData")
tests_mmp <- data_mmp |>
  mutate(COVER = n.points / total.points) |>
  nest(.by = c("Reef", "REEF_ZONE", "DEPTH", "AIMS_REEF_NAME", "Site", "Transect", "Dist.number")) |>
  mutate(BeforeAfter = map(
    .x = data,
    .f = ~ before_after_tests(.x)
  ))
save(tests_mmp, file = "../data/processed/q2_tests_mmp_2.RData")
## ----end



## Instances where there are no before cases
## ---- before_after_tests_no_before_2_mmp
load(file = "../data/processed/q2_tests_mmp_2.RData")
tests_mmp |>
  dplyr::select(-data) |>
  unnest(c(BeforeAfter)) |>
  filter(before_flag == "no before") |>
  mutate(across(where(is.numeric), round, 3)) |>
  as.data.frame() |>
  datatable(class = c("compact", "white-space: nowrap"))
## ----end


## Instances where there are multiple before cases
## ---- before_after_tests_multiple_before_2_mmp
load(file = "../data/processed/q2_tests_mmp_2.RData")
tests_mmp |>
  dplyr::select(-data) |>
  unnest(c(BeforeAfter)) |>
  filter(before_flag == "multiple before") |>
  mutate(across(where(is.numeric), round, 3)) |>
  as.data.frame() |>
  datatable(class = c("compact", "white-space: nowrap"))
## ----end

## Instances where there After Disturbance Type is "n"
## ---- before_after_tests_n_after_2_mmp
load(file = "../data/processed/q2_tests_mmp_2.RData")
tests_mmp |>
  dplyr::select(-data) |>
  unnest(c(BeforeAfter)) |>
  filter(n_flag == "after disturb n") |>
  mutate(across(where(is.numeric), round, 3)) |>
  as.data.frame() |>
  datatable(class = "compact")
## ----end

## Instances where there are no after cases
## ---- before_after_tests_no_after_2_mmp
load(file = "../data/processed/q2_tests_mmp_2.RData")
tests_mmp |>
  dplyr::select(-data) |>
  unnest(c(BeforeAfter)) |>
  filter(after_flag == "no after") |>
  mutate(across(where(is.numeric), round, 3)) |>
  as.data.frame() |>
  datatable(class = c("compact", "white-space: nowrap"))
## ----end


data_mmp |>
  filter(AIMS_REEF_NAME == "Snapper Island", REEF_ZONE == "North",
    DEPTH == 5, Dist.number == 1, SITE_NO == 2, TRANSECT_NO == 5) |>
  as.data.frame()



## Instances where cover was greater after than before
## ---- before_after_tests_cover_greater_after_2_mmp
load(file = "../data/processed/q2_tests_mmp_2.RData")
tests_mmp |>
  dplyr::select(-data) |>
  unnest(c(BeforeAfter)) |>
  filter(cover_flag == "cover increased") |>
  mutate(across(where(is.numeric), round, 3)) |>
  as.data.frame() |>
  datatable(class = c("compact", "white-space: nowrap"))
## ----end

## ---- make_excludes_2_mmp
load(file = "../data/processed/q2_tests_mmp_2.RData")
excludes_mmp <- tests_mmp |>
  dplyr::select(-data) |>
  unnest(c(BeforeAfter)) |>
  filter(!is.na(before_flag) |
    !is.na(after_flag) |
    !is.na(n_flag)) |>
  dplyr::select(Transect, REEF_ZONE, DEPTH, Dist.number) |>
  distinct() |>
  mutate(EXCLUDE = TRUE)
save(excludes_mmp, file = "../data/processed/q2_excludes_mmp_2.RData")
## ----end

## ---- exclude_data_mmp
load(file = "../data/processed/q2_data_mmp_2.RData")
if (nrow(excludes_mmp) > 0) {
  data_mmp <-
    data_mmp |>
    left_join(excludes_mmp) |>
    filter(is.na(EXCLUDE))
}
## ----end

## ---- encode_disturbances_mmp
data_mmp <- data_mmp |>
  mutate(across(c(s, c, d, b, u), \(x) ifelse(Dist.time == "Before", 0, x))) |>
  mutate(Dist = case_when(
    s == 1 & (c + d + b + u == 0) ~ "s",
    s > 1 & (c + d + b + u == 0) ~ "ss",
    s > 0 & (c + d + b + u > 0) ~ "sm",
    c == 1 & (s + d + b + u == 0) ~ "c",
    c > 1 & (s + d + b + u == 0) ~ "cc",
    c > 0 & (s + d + b + u > 0) ~ "cm",
    d == 1 & (s + c + b + u == 0) ~ "d",
    d > 1 & (s + c + b + u == 0) ~ "df",
    d > 0 & (s + c + b + u > 0) ~ "dm",
    b == 1 & (s + c + d + u == 0) ~ "b",
    b > 1 & (s + c + d + u == 0) ~ "bb",
    b > 0 & (s + c + d + u > 0) ~ "bm",
    u == 1 & (s + c + d + b == 0) ~ "u",
    u > 1 & (s + c + d + b == 0) ~ "uu",
    u > 0 & (s + c + d + b > 0) ~ "um",
    .default = "Before"
  )) |>
  mutate(Dist = forcats::fct_relevel(Dist, "Before")) |>
  mutate(Dist.type = ifelse(s + c + d + b + u > 1, "Cumulative", "Single")) |>
  mutate(across(c(s, c, d, b, u), \(x) ifelse(x > 1, 1, x))) |>
  mutate(across(c(s, c, d, b, u), as.factor)) |>
  mutate(event = paste(Transect, Dist.number)) |>
  mutate(SSDist = paste(SecShelf, Dist))
save(data_mmp, file = "../data/processed/data_mmp_q2.RData")
## ----end

## ---- q2_view_data_mmp
load(file = "../data/processed/data_mmp_q2.RData")
data_mmp |>
  dplyr::relocate(A_SECTOR, SHELF, AIMS_REEF_NAME, REPORT_YEAR) |> 
  datatable(
    class = c("compact", "white-space: nowrap"),
    filter = "top",
    extensions = c("Buttons","FixedColumns", "KeyTable"),
    options = list(
      dom = "Bfrtip",
      scrollX = TRUE,
      fixedColumns = list(leftColumns = 5),
      keys = TRUE,       ## move via arrow keys
      autoWidth = TRUE,
      buttons = list(
        list(
          extend = "collection",
          buttons = c("copy", "csv", "excel", "pdf", "print"),
          text = "Download"
        )
      ),
      columnDefs = list(list(width = "10px", targets = 0))
    )
  )
## ----end


### Finally, put LTMP and MMP together

## ---- combine_ltmp_mmp
data <-
  data |>
  mutate(Dist.number = as.numeric(as.character(Dist.number))) |>
  bind_rows(data_mmp)
save(data, file = "../data/processed/data_q2_combined.RData")
## ----end
