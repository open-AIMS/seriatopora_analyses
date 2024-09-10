## Retrieve the data
## ---- q2_data_1.5
load(file = "../data/processed/data_q2_combined.RData")
## ----end



## Strip out all unnecessary variables

full_data <- data
data <- full_data
rerun_models <- FALSE

## This model will explore the broad difference in cover between
## before and after (Dist.time)


## Prepare data
{
  ## ---- q2_prepare_data_1.5
  ## Focus on only the necessary variables
  data <- data |>
    mutate(SecShelf = paste(A_SECTOR, SHELF)) |> 
    mutate(zone_depth = factor(paste0(REEF_ZONE, DEPTH))) |> 
    dplyr::select(
      n.points, total.points, Dist.time, s, c, d, b, u,
      AIMS_REEF_NAME, Site, Transect, SecShelf, Dist.number, Dist,
      SecShelf, SSDist, event, zone_depth
    ) |>
    filter(!SecShelf %in% c("CG M", "PC M", "PC O")) |>
    droplevels()
  ## mutate(across(c(s, c, d, b, u), \(x) ifelse(Dist.time == "Before", 0, x))) 
  ## ----end
}

## Fit the models
{
  ## Raw means
  {
    ## ---- q2_raw_fitted_1.5
    cellmeans_summ_raw <-
      data |>
      group_by(Dist, SecShelf) |> 
      reframe(
        type = c("mean", "median"),
        Mean = c(mean(n.points / total.points), median(n.points / total.points)),
        SD = c(sd(n.points / total.points), MAD = mad(n.points / total.points)),
        N = c(n(), n())
      ) |>
      mutate(
        lower = Mean - 2 * (SD / sqrt(N)),
        upper = Mean + 2 * (SD / sqrt(N))
      ) |>
      separate(SecShelf, into = c("A_SECTOR", "SHELF")) |>
        mutate(Dist.time = ifelse(Dist == "Before", "Before", "After")) |>
        mutate(Dist.time = factor(Dist.time, levels = c("Before", "After")))
    save(cellmeans_summ_raw, file = "../data/modelled/cellmeans_summ_raw_1.5.RData")
    p <-
    cellmeans_summ_raw |>
      ggplot(aes(y = Mean, x = Dist, colour = Dist.time, shape = type)) +
      geom_pointrange(aes(ymin = lower, ymax = upper)) +
      scale_shape_manual(values = c(16, 21)) +
      facet_wrap(A_SECTOR ~ SHELF,
        scales = "free",
        labeller = label_wrap_gen(multi_line = FALSE)
      )
      ggsave(
        filename = "../docs/analysis_files/figure-html/partial_raw_1.5.png",
        p,
        width = 15, height = 10
      )
    ## ----end
  }
  ## contrasts
  {
##     a <-
##     data |>
##       dplyr::select(SecShelf, event, Dist, s, c, d, b, u, n.points, total.points) |>
##       rowwise() |>
##       mutate(Dist2 = paste(c("s", "c", "d", "b", "u")[which(c(s, c, d, b, u) == 1)], collapse = "")) |>
##       ungroup() |>
##       separate_longer_position(Dist2, width = 1, keep_empty = TRUE) |>
##       mutate(cummulative = ifelse(Dist != "Before" & str_length(Dist) > 1, TRUE, FALSE)) |>
##       group_by(SecShelf, event, Dist, Dist2) |>
##       summarise(
##         Mean = mean(n.points / total.points),
##         Median = median(n.points / total.points),
##         SD = sd(n.points / total.points),
##         N = n(),
##         Dist2 = paste(unique(Dist2), collapse = ",")
##       ) |>
##       ungroup() |>
##         mutate(Dist2 = ifelse(Dist2 == "NA", as.character(Dist), Dist2))  

##     b1 <-
##       a |>
##       mutate(Values = Mean) |>
##       mutate(Dist = Dist2) |>
##       group_by(SecShelf, event) |>
##       mutate(NN = n()) |>
##       ungroup() |>
##       filter(NN > 1) |> 
##         nest(.by = c(SecShelf, event)) |>
##         mutate(eff = map(
##           .x = data,
##           .f = ~ before_vs_afters(.x)
##         )) |>
##         dplyr::select(-data) |>
##           unnest(c(eff)) |>
##           mutate(Pred = 100 * (Values - 1))

## b <- a |> 
##     group_by(SecShelf, Dist2) |>
##       summarise(
##         Mean = mean(Mean),
##         Median = median(Median),
##         SD = mean(SD),
##         N = sum(N)
##       ) |>
##       ungroup() 

##     bb <- b |>
##       mutate(Dist = Dist2, Values = Mean) |>
##         nest(.by = c(SecShelf)) |>
##         mutate(eff = map(
##           .x = data,
##           .f = ~ before_vs_afters(.x)
##         )) |>
##         dplyr::select(-data) |>
##           unnest(c(eff)) |>
##           mutate(Pred = 100 * (Values - 1))
    
        
  }
  ## glmmTMB
  {
    ## fit model
    {
      if (rerun_models) {
        ## ---- q2_glmmTMB_fitted_1.5
        mod_glmmTMB <- glmmTMB(
          cbind(n.points, total.points - n.points) ~ 1 + (s + c + d + b + u) +
            (1 | SecShelf:(s + c + d + b + u)) +
            (1 | event) +
            (1 | AIMS_REEF_NAME) +
            (1 | zone_depth) +
            (1 | Site) +
            (1 | Transect),
          data = data,
          family = "betabinomial",
          REML = TRUE
        )
        save(mod_glmmTMB, file = "../data/modelled/mod_glmmTMB_1.5.RData")
        ## ----end
      }
    }
    ## Explore residuals (model validation)
    {
      ## ---- q2_glmmTMB_validation_1.5
      load(file = "../data/modelled/mod_glmmTMB_1.5.RData")
      resids <- simulateResiduals(mod_glmmTMB, plot = FALSE)
      wrap_plots(
        wrap_elements(~testUniformity(resids)),
        wrap_elements(~ plotResiduals(resids)),
        wrap_elements(~testDispersion(resids)),
        nrow = 1
      ) 
      ## ----end
    }
    ## Model summary
    {
      ## ---- q2_glmmTMB_summary_1.5
      load(file = "../data/modelled/mod_glmmTMB_1.5.RData")
      summ <- 
        bind_rows(
          model_summary(mod_glmmTMB) |> mutate(model = "Raw"),
          model_summary(mod_glmmTMB, exponentiate = TRUE) |> mutate(model = "Exponentiated")
        ) |>
        mutate(across(c(estimate, conf.low, conf.high), \(x) sprintf("%.3f", x))) |> 
        mutate(estimate = paste0(estimate, "<br>[", conf.low, ",", conf.high, "]")) |>
        mutate(estimate = str_replace(estimate, "<br>\\[NA,NA\\]","")) |> 
        dplyr::select(-conf.low, -conf.high) |> 
        pivot_longer(cols = c("estimate"),
          names_to = "names", values_to = "values") |> 
        pivot_wider(id_cols = everything(), names_from = "model", values_from = "values") |>
        dplyr::select(-effect, -names) 

      summ_tt <- summ |> 
        tinytable::tt() 
      save(summ, file = "../data/summ_glmmTMB_1.5.RData")
      save(summ_tt, file = "../data/summ_tt_glmmTMB_1.5.RData")
      ## ----end
      ## ---- q2_glmmTMB_summary_1.5_print
      load(file = "../data/summ_tt_glmmTMB_1.5.RData")
      summ_tt
      ## ----end
    }
    ## Partial effects
    {
      ## ---- q2_glmmTMB_parial_effects_1.5
      load(file = "../data/modelled/mod_glmmTMB_1.5.RData")

      ## 1. focus on just necessary fields
      ## 2. create a field (Dist2) that concatenates all Dists based on the 0's and 1's
      newdata <-
        data |>
        mutate(zone_depth = NA) |> 
        dplyr::select(zone_depth,SecShelf, Dist, s, c, d, b, u) |>
        ## filter(SecShelf == "CA M") |> 
        distinct() |>
        rowwise() |>
        mutate(Dist2 = paste(c("s", "c", "d", "b", "u")[which(c(s, c, d, b, u) == 1)], collapse = "")) |>
        ungroup()

      newdata <- bind_rows(
        ## 1. based on the Dist2 field, duplicate the rows into their components
        ## 2. create a field to indicate whether the row represents a cummulative (multiple disturbances)
        newdata |>
          mutate(Dist3 = Dist2) |> 
          separate_longer_position(Dist2, width = 1, keep_empty = TRUE) |>
          mutate(cummulative = ifelse(Dist != "Before" & str_length(Dist)>1, TRUE, FALSE)),
        newdata |>
          mutate(Dist3 = Dist2) |> 
          separate_longer_position(Dist2, width = 1, keep_empty = TRUE) |>
          mutate(s = factor(ifelse(!str_detect(Dist2, "s") | is.na(Dist2), 0, 1))) |>
          mutate(c = factor(ifelse(!str_detect(Dist2, "c") | is.na(Dist2), 0, 1))) |>
          mutate(d = factor(ifelse(!str_detect(Dist2, "d") | is.na(Dist2), 0, 1))) |>
          mutate(b = factor(ifelse(!str_detect(Dist2, "b") | is.na(Dist2), 0, 1))) |>
          mutate(u = factor(ifelse(!str_detect(Dist2, "u") | is.na(Dist2), 0, 1))) |>
          mutate(cummulative = FALSE)) |>
        mutate(Dist2 = ifelse(is.na(Dist2), "Before", Dist2)) |>
          dplyr::select(-Dist) |>
          distinct() |>
          arrange(Dist2) |>
          mutate(AIMS_REEF_NAME = NA, Site = NA, Transect = NA, event = NA) |>
          ## mutate(cummulative = ifelse(s + c + d + b + u > 1), TRUE, FALSE) |>
          ## mutate(cummulative = ifelse(str_detect(Dist, "m") | str_length(Dist) > 1, TRUE, FALSE)) |>
          filter(!is.na(s)) |>
          mutate(Dist3 = ifelse(cummulative, Dist3, Dist2))
      
        p <- predict(mod_glmmTMB, newdata = newdata, se.fit = TRUE)
        newdata <- newdata |>
          bind_cols(fit = p$fit, se = p$se.fit) |>
          mutate(Pred = plogis(fit), lower = plogis(fit - 2 * se), upper = plogis(fit + 2 * se)) |>
          mutate(
            Dist = ifelse(Dist2 == "" | is.na(Dist2), "Before", Dist2),
            Dist = forcats::fct_relevel(Dist, "Before"),
            Dist3 = forcats::fct_relevel(Dist3, "Before")
          )
      

      cellmeans_summ_glmmTMB <- newdata |>
        dplyr::select(SecShelf, Dist, Dist3, cummulative, Pred, lower, upper) |>
        separate(SecShelf, into = c("A_SECTOR", "SHELF"), sep = " ") 
      save(cellmeans_summ_glmmTMB, file = "../data/modelled/cellmeans_summ_glmmTMB_1.5.RData")

      p <-
      cellmeans_summ_glmmTMB |>
        mutate(Dist.time = ifelse(Dist == "Before", "Before", "After")) |>
        mutate(Dist.time = factor(Dist.time, levels = c("Before", "After"))) |>
        mutate(Dist = forcats::fct_relevel(Dist, "Before")) |> 
        ggplot(aes(y = Pred, x = Dist, colour = Dist.time, shape = cummulative)) +
        geom_pointrange(aes(ymin = lower, ymax = upper), position = position_dodge(width = 0.5)) +
        scale_shape_manual(values = c(16, 21)) +
        facet_wrap(A_SECTOR ~ SHELF,
          scales = "free",
          labeller = label_wrap_gen(multi_line = FALSE)
        )
      p
      ggsave(
        filename = "../docs/analysis_files/figure-html/partial_glmmTMB_1.5a.png",
        p,
        width = 10, height = 8
      )

      p <-
      cellmeans_summ_glmmTMB |>
        mutate(Dist.time = ifelse(Dist == "Before", "Before", "After")) |>
        mutate(Dist.time = factor(Dist.time, levels = c("Before", "After"))) |>
        mutate(Dist = forcats::fct_relevel(Dist, "Before")) |> 
        ggplot(aes(y = Pred, x = Dist3, colour = Dist.time, shape = cummulative)) +
        geom_pointrange(aes(ymin = lower, ymax = upper), position = position_dodge(width = 0.5)) +
        scale_shape_manual(values = c(16, 21)) +
        facet_wrap(A_SECTOR ~ SHELF,
          scales = "free",
          labeller = label_wrap_gen(multi_line = FALSE)
        )
      ggsave(
        filename = "../docs/analysis_files/figure-html/partial_glmmTMB_1.5b.png",
        p,
        width = 10, height = 8
      )

      ## ----end
    }
    ## Contrasts
    {
      ## ---- q2_glmmTMB_contrasts_1.5
      eff_summ_glmmTMB <-
        cellmeans_summ_glmmTMB |>
        mutate(Values = Pred) |>
        nest(.by = c(A_SECTOR, SHELF)) |>
        mutate(eff = map(
          .x = data,
          .f = ~ before_vs_afters(.x)
        )) |>
        dplyr::select(-data) |>
        unnest(c(eff))
      save(eff_summ_glmmTMB, file = "../data/modelled/eff_summ_glmmTMB_7.1.RData")
      eff_summ_glmmTMB |> ggplot(aes(y = Dist, x = Values)) +
        geom_vline(xintercept = 1, linetype = "dashed") +
        geom_point() +
        facet_wrap(A_SECTOR ~ SHELF, scales = "free")

      eff_summ_glmmTMB <-
        cellmeans_summ_glmmTMB |>
        mutate(Dist = Dist3) |> 
        mutate(Values = Pred) |>
        nest(.by = c(A_SECTOR, SHELF)) |>
        mutate(eff = map(
          .x = data,
          .f = ~ before_vs_afters(.x)
        )) |>
        dplyr::select(-data) |>
        unnest(c(eff))
      eff_summ_glmmTMB |> ggplot(aes(y = Dist, x = Values)) +
        geom_vline(xintercept = 1, linetype = "dashed") +
        geom_point() +
        facet_wrap(A_SECTOR ~ SHELF, scales = "free")

      eff_summ_glmmTMB <-
        cellmeans_summ_glmmTMB |>
        mutate(Dist = Dist3) |> 
        mutate(Values = Pred) |>
        nest(.by = c(A_SECTOR, SHELF)) |>
        mutate(eff = map(
          .x = data,
          .f = ~ before_vs_afters(.x)
        )) |>
        dplyr::select(-data) |>
          unnest(c(eff)) |>
          mutate(cummulative = ifelse(str_length(Dist) > 1, TRUE, FALSE)) |>
          mutate(Dist3 = Dist) |> 
          separate_longer_position(Dist, width = 1, keep_empty = TRUE) 
        
      eff_summ_glmmTMB |> ggplot(aes(y = Dist, x = Values,
        color = Dist3, #shape = cummulative,
        size = factor(str_length(Dist3)))) +
        geom_vline(xintercept = 1, linetype = "dashed") +
          geom_point(position = position_dodge(width = 0.5)) +
          scale_shape_manual(values = c(16, 21)) +
          scale_size_discrete(str_wrap("Number of disturbances types", 15)) +
          scale_colour_discrete("Disurbances") +
        scale_x_continuous("Effect size (fold scale decline after disturbance(x))",
          trans = scales::log2_trans(),
          labels = function(x) 1/x) +
        ## facet_wrap(A_SECTOR ~ SHELF,
        ##   scales = "free",
        ##   labeller = label_wrap_gen(multi_line = FALSE)
        ## )
      facet_grid(A_SECTOR ~ SHELF,
          scales = "free",
          labeller = label_wrap_gen(multi_line = FALSE)
        )
      ## ----end
    }
  }
  ## brms
  {
    ## fit model
    {
      if (rerun_models) {
        ## ---- q2_brm_fitted_1.5
        form <- bf(
          n.points | trials(total.points) ~ 1 + (s + c + d + b + u) +
            (1 | SecShelf:(s + c + d + b + u)) +
            (1 | event) +
            (1 | AIMS_REEF_NAME) +
            (1 | zone_depth) +
            (1 | Site) +
            (1 | Transect),
          family = "beta_binomial"
        )
        data |>
          group_by(Dist) |> 
          summarise(qlogis(mean(n.points / total.points)))
        priors <- prior(normal(-5, 1), class = "Intercept") +
          prior(normal(0, 2), class = "b") +
          prior(student_t(3, 0, 1), class = "sd") +
          prior(gamma(0.01, 0.01), class = "phi")
        mod_brm <- brm(
          formula = form,
          data = data,
          prior = priors,
          iter = 5000, warmup = 1000,
          ##iter = 50, warmup = 10,
          chains = 3, cores = 3,
          sample_prior = "yes",
          thin = 10,
          backend = "cmdstanr",
          control = list(adapt_delta = 0.99)
        )
        summary(mod_brm)
        save(mod_brm, file = "../data/modelled/mod_brm_1.5.RData")
        ## ----end
      }
    }
    ## MCMC diagnostics
    {
      ## ---- q2_brm_trace_1.5a
      load(file = "../data/modelled/mod_brm_1.5.RData")
      rstan::stan_trace(mod_brm$fit)
      rstan::stan_ac(mod_brm$fit)
      rstan::stan_ess(mod_brm$fit) + rstan::stan_rhat(mod_brm$fit)
      ## ----end
      
      ## ---- q2_brm_trace_1.5
      load(file = "../data/modelled/mod_brm_1.5.RData")
      ggsave(
        filename = "../docs/analysis_files/figure-html/stan_trace_brm_1.5.png",
        rstan::stan_trace(mod_brm$fit),
        width = 15, height = 8
      )
      ggsave(
        filename = "../docs/analysis_files/figure-html/stan_ac_brm_1.5.png",
        rstan::stan_ac(mod_brm$fit),
        width = 15, height = 8
      )
      ggsave(
        filename = "../docs/analysis_files/figure-html/stan_ess_brm_1.5.png",
        rstan::stan_ess(mod_brm$fit) + rstan::stan_rhat(mod_brm$fit),
        width = 15, height = 8
      )
      ## ----end
    }
    ## Model validation
    {
      ## ---- q2_brm_validation_1.5a
      load(file = "../data/modelled/mod_brm_1.5.RData")
      mod_brm |> pp_check(type = "dens_overlay", ndraws = 100) + scale_x_log10() +
        mod_brm |> pp_check(type = "loo_pit_overlay", ndraws = 100)
      ## ----end

      ## ---- q2_brm_validation_1.5
      load(file = "../data/modelled/mod_brm_1.5.RData")
      ggsave(
        filename = "../docs/analysis_files/figure-html/pcc_brm_1.5.png",
        mod_brm |> pp_check(type = "dens_overlay", ndraws = 100) + scale_x_log10() +
          mod_brm |> pp_check(type = "loo_pit_overlay", ndraws = 100),
        width = 10,
        height = 5)
      ## ----end
    }
    ## DHARMa residuals
    {
      ## ---- q2_brm_dharma_1.5a
      load(file = "../data/modelled/mod_brm_1.5.RData")
      resids <- make_brms_dharma_res(
        mod_brm,
        integerResponse = FALSE
      )
      wrap_elements(~ testUniformity(resids)) +
        wrap_elements(~ plotResiduals(resids, form = factor(rep(1, nrow(mod_brm$data))))) +
        wrap_elements(~ plotResiduals(resids)) +
        wrap_elements(~ testDispersion(resids))
      ## ----end

      ## ---- q2_brm_dharma_1.5
      load(file = "../data/modelled/mod_brm_1.5.RData")
      resids <- make_brms_dharma_res(
        mod_brm,
        integerResponse = FALSE
      )
      ggsave(
        filename = "../docs/analysis_files/figure-html/dharma_brm_1.5.png",
        wrap_elements(~ testUniformity(resids)) +
          wrap_elements(~ plotResiduals(resids, form = factor(rep(1, nrow(mod_brm$data))))) +
          wrap_elements(~ plotResiduals(resids)) +
          wrap_elements(~ testDispersion(resids)),
        width = 12,
        height = 8,
        dpi = 100
      )
      ## ----end
      
    }
    ## Model summary
    {
      ## ---- q2_brm_summary_1.5
      load(file = "../data/modelled/mod_brm_1.5.RData")

      summ <- 
      bind_rows(
        model_summary(mod_brm) |> mutate(model = "Raw"),
        model_summary(mod_brm, exponentiate = TRUE) |> mutate(model = "Exponentiated")
      ) |>
        mutate(across(c(estimate, conf.low, conf.high), \(x) sprintf("%.3f", x))) |> 
        mutate(estimate = paste0(estimate, "<br>[", conf.low, ",", conf.high, "]")) |>
        mutate(estimate = str_replace(estimate, "<br>\\[NA,NA\\]","")) |> 
        dplyr::select(-conf.low, -conf.high) |> 
        pivot_longer(cols = c("estimate"),
          names_to = "names", values_to = "values") |> 
        pivot_wider(id_cols = everything(), names_from = "model", values_from = "values") |>
        dplyr::select(-effect, -names) 

      summ_tt <- summ |> 
        tinytable::tt() 
      save(summ, file = "../data/summ_brm_1.5.RData")
      save(summ_tt, file = "../data/summ_tt_brm_1.5.RData")
      ## ----end
      ## ---- q2_brm_summary_1.5_print
      load(file = "../data/summ_tt_brm_1.5.RData")
      summ_tt
      ## ----end
    }

    ## Partial plots -  Type 1. sector/shelf re_formula = NULL
    {
      ## ---- q2_brm_partial_effects_1.5
      load(file = "../data/modelled/mod_brm_1.5.RData")
      newdata <- brm_generate_newdata(data)
      pred <- t(brms::posterior_epred(mod_brm,
        newdata = newdata,
        re_formula = NULL,
        se.fit = TRUE,
        allow_new_levels = TRUE
      ))
      cellmeans_brm <- brm_generate_cellmeans(newdata, pred)

      save(cellmeans_brm, file = "../data/modelled/cellmeans_brm_1.5_type_1.RData")

      cellmeans_summ_brm <- cellmeans_brm |>
        dplyr::select(-.draws) |>
        dplyr::select(A_SECTOR, SHELF, Dist, Dist2, Dist3, cummulative, Pred) |>
        group_by(A_SECTOR, SHELF, Dist, Dist2, Dist3, cummulative) |>
        summarise_draws(median, HDInterval::hdi) |>
        dplyr::select(Dist, median, lower, upper)

      save(cellmeans_summ_brm, file = "../data/modelled/cellmeans_summ_brm_1.5_type_1.RData")

      p <- brm_partial_plot(cellmeans_summ_brm)
      ggsave(
        filename = "../docs/analysis_files/figure-html/partial_brm_1.5_type_1.png",
        p,
        width = 15, height = 12
      )
      ## ----end
    }
    ## Partial plots -  Type 2. sector/shelf re_formula = NA
    {
      ## ---- q2_brm_partial_effects_1.5_type_2
      load(file = "../data/modelled/mod_brm_1.5.RData")
      newdata <- brm_generate_newdata(data)
      pred <- t(brms::posterior_epred(mod_brm,
        newdata = newdata,
        re_formula = NA,
        se.fit = TRUE,
        allow_new_levels = TRUE
      ))
      cellmeans_brm <- brm_generate_cellmeans(newdata, pred)

      save(cellmeans_brm, file = "../data/modelled/cellmeans_brm_1.5_type_2.RData")

      cellmeans_summ_brm <- cellmeans_brm |>
        dplyr::select(-.draws) |>
        dplyr::select(A_SECTOR, SHELF, Dist, Dist2, Dist3, cummulative, Pred) |>
        group_by(A_SECTOR, SHELF, Dist, Dist2, Dist3, cummulative) |>
        summarise_draws(median, HDInterval::hdi) |>
        dplyr::select(Dist, median, lower, upper)

      save(cellmeans_summ_brm, file = "../data/modelled/cellmeans_summ_brm_1.5_type_2.RData")

      p <- brm_partial_plot(cellmeans_summ_brm)
      ggsave(
        filename = "../docs/analysis_files/figure-html/partial_brm_1.5_type_2.png",
        p,
        width = 15, height = 12
      )
      ## ----end
    }
    ## Partial plots -  Type 3. sector/shelf re_formula = (1|SecShelf:(s+b+c+d+u))
    {
      ## ---- q2_brm_partial_effects_1.5_type_3
      load(file = "../data/modelled/mod_brm_1.5.RData")
      newdata <- brm_generate_newdata(data)

      pred <- t(brms::posterior_epred(mod_brm,
        newdata = newdata,
        re_formula = ~ (1 | SecShelf:(s + c + d + b + u))
      ))
      cellmeans_brm <- brm_generate_cellmeans(newdata, pred)

      save(cellmeans_brm, file = "../data/modelled/cellmeans_brm_1.5_type_3.RData")

      cellmeans_summ_brm <- cellmeans_brm |>
        dplyr::select(-.draws) |>
        dplyr::select(A_SECTOR, SHELF, Dist, Dist2, Dist3, cummulative, Pred) |>
        group_by(A_SECTOR, SHELF, Dist, Dist2, Dist3, cummulative) |>
        summarise_draws(median, HDInterval::hdi) |>
        dplyr::select(Dist, median, lower, upper)

      save(cellmeans_summ_brm, file = "../data/modelled/cellmeans_summ_brm_1.5_type_3.RData")

      p <- brm_partial_plot(cellmeans_summ_brm)
      ggsave(
        filename = "../docs/analysis_files/figure-html/partial_brm_1.5_type_3.png",
        p,
        width = 15, height = 12
      )
      ## ----end
    }

    ## Contrasts - Type 1. sector/shelf re_formula = NULL
    {
      ## ---- q2_brm_contrasts_1.5_type_1
      load(file = "../data/modelled/cellmeans_brm_1.5_type_1.RData")

      eff_brm <- brm_calc_effect(cellmeans_brm)
      save(eff_brm, file = "../data/modelled/eff_brm_1.5_type_1.RData")

      p <- brm_effects_plot(eff_brm)
      ggsave(
        filename = "../docs/analysis_files/figure-html/contrasts_brm_1.5_type_1.png",
        p,
        width = 15, height = 12
      )
    }
    ## Contrasts - Type 2. sector/shelf re_formula = NA
    {
      ## ---- q2_brm_contrasts_1.5_type_2
      load(file = "../data/modelled/cellmeans_brm_1.5_type_2.RData")

      eff_brm <- brm_calc_effect(cellmeans_brm)
      save(eff_brm, file = "../data/modelled/eff_brm_1.5_type_2.RData")

      p <- brm_effects_plot(eff_brm)
      ggsave(
        filename = "../docs/analysis_files/figure-html/contrasts_brm_1.5_type_2.png",
        p,
        width = 15, height = 12
      )
    }
    ## Contrasts - Type 3. sector/shelf re_formula = (1|SecShelf:(s+b+c+d+u)) 
    {
      ## ---- q2_brm_contrasts_1.5_type_3
      load(file = "../data/modelled/cellmeans_brm_1.5_type_3.RData")

      eff_brm <- brm_calc_effect(cellmeans_brm)
      save(eff_brm, file = "../data/modelled/eff_brm_1.5_type_3.RData")

      p <- brm_effects_plot(eff_brm)
      ggsave(
        filename = "../docs/analysis_files/figure-html/contrasts_brm_1.5_type_3.png",
        p,
        width = 15, height = 12
      )
    }

    ## Contrast - GBR Type 1. sector/shelf re_formula = NULL
    {
      ## ---- q2_brm_contrasts_gbr_1.5_type_1
      load(file = "../data/modelled/mod_brm_1.5.RData")
      load(file = "../data/modelled/cellmeans_brm_1.5_type_1.RData")
      cellmeans_brm <- cellmeans_brm |>
        group_by(A_SECTOR, SHELF, Dist, .draws) |>
        summarise(value = mean(value)) |>
        ungroup() 

      save(cellmeans_brm, file = "../data/modelled/cellmeans_brm_1.5_gbr_type_1.RData")
      
      eff_brm <- brm_calc_effect_hier(cellmeans_brm)

      eff_brm <-
        eff_brm |>
        group_by(.draws, Dist) |>
        summarise(Values = mean(Values)) |>
        ungroup()
      save(eff_brm, file = "../data/modelled/eff_brm_1.5_gbr_type_1.RData")
      load(file = "../data/modelled/eff_brm_1.5_gbr_type_1.RData")

      dist_order <- tribble(
        ~levels, ~labels,
        "b",     "Bleaching",
        "c",     "COTS",
        "s",     "Storms",
        "d",     "Disease",
        "u",     "Unknown",
      )
      p <- brm_effects_plot_gbr(eff_brm, height = 1.3, p = 0.5,
        dist_order = dist_order, label_offset = 0.2
        ) +
        coord_cartesian(clip = "on", xlim = c(0.1, 2.1))

      ggsave(
        filename = "../docs/analysis_files/figure-html/contrasts_gbr_combined_disturbances_brm_1.5_type_1.png",
        p,
        width = 8, height = 8/1.6
      )

    }
    ## Contrast - GBR Type 2. sector/shelf re_formula = NA
    {
      ## ---- q2_brm_contrasts_gbr_1.5_type_2
      load(file = "../data/modelled/mod_brm_1.5.RData")
      load(file = "../data/modelled/cellmeans_brm_1.5_type_2.RData")
      cellmeans_brm <- cellmeans_brm |>
        group_by(A_SECTOR, SHELF, Dist, .draws) |>
        summarise(value = mean(value)) |>
        ungroup() 

      save(cellmeans_brm, file = "../data/modelled/cellmeans_brm_1.5_gbr_type_2.RData")
      
      eff_brm <- brm_calc_effect_hier(cellmeans_brm)

      eff_brm <-
        eff_brm |>
        group_by(.draws, Dist) |>
        summarise(Values = mean(Values)) |>
        ungroup()
      save(eff_brm, file = "../data/modelled/eff_brm_1.5_gbr_type_2.RData")
      load(file = "../data/modelled/eff_brm_1.5_gbr_type_2.RData")

      dist_order <- tribble(
        ~levels, ~labels,
        "b",     "Bleaching",
        "c",     "COTS",
        "s",     "Storms",
        "d",     "Disease",
        "u",     "Unknown",
      )
      p <- brm_effects_plot_gbr(eff_brm, height = 1.3, p = 0.5,
        dist_order = dist_order, label_offset = 0.2
        ) +
        coord_cartesian(clip = "on", xlim = c(0.04, 2.1))

      ggsave(
        filename = "../docs/analysis_files/figure-html/contrasts_gbr_combined_disturbances_brm_1.5_type_2.png",
        p,
        width = 8, height = 8/1.6
      )
    }
    ## Contrast - GBR Type 3. sector/shelf re_formula =  (1|SecShelf:(s+b+c+d+u)
    {
      ## ---- q2_brm_contrasts_gbr_1.5_type_3
      load(file = "../data/modelled/mod_brm_1.5.RData")
      load(file = "../data/modelled/cellmeans_brm_1.5_type_3.RData")
      cellmeans_brm <- cellmeans_brm |>
        group_by(A_SECTOR, SHELF, Dist, .draws) |>
        summarise(value = mean(value)) |>
        ungroup() 

      save(cellmeans_brm, file = "../data/modelled/cellmeans_brm_1.5_gbr_type_3.RData")
      
      eff_brm <- brm_calc_effect_hier(cellmeans_brm)

      eff_brm <-
        eff_brm |>
        group_by(.draws, Dist) |>
        summarise(Values = mean(Values)) |>
        ungroup()
      save(eff_brm, file = "../data/modelled/eff_brm_1.5_gbr_type_3.RData")
      load(file = "../data/modelled/eff_brm_1.5_gbr_type_3.RData")

      dist_order <- tribble(
        ~levels, ~labels,
        "b",     "Bleaching",
        "c",     "COTS",
        "s",     "Storms",
        "d",     "Disease",
        "u",     "Unknown",
      )
      p <- brm_effects_plot_gbr(eff_brm, height = 1.3, p = 0.5,
        dist_order = dist_order, label_offset = 0.2
        ) +
        coord_cartesian(clip = "on", xlim = c(0.1, 2.1))

      ggsave(
        filename = "../docs/analysis_files/figure-html/contrasts_gbr_combined_disturbances_brm_1.5_type_3.png",
        p,
        width = 8, height = 8/1.6
      )
    }

    ## OLD Contrast - GBR Type 1. sector/shelf re_formula = NULL
    {
      ## ---- Old q2_brm_contrasts_gbr_1.5_type_1
      if(1 != 2) {
        load(file = "../data/modelled/mod_brm_1.5.RData")
        newdata <- brm_generate_newdata(data, GBR = TRUE)
        pred <- t(brms::posterior_epred(mod_brm,
          newdata = newdata,
          re_formula = NULL,
          allow_new_levels = TRUE
        ))

        cellmeans_brm <- brm_generate_cellmeans(newdata, pred)
        save(cellmeans_brm, file = "../data/modelled/cellmeans_brm_1.5_gbr_type_1.RData")

        eff_brm <- brm_calc_effect(cellmeans_brm, GBR = TRUE)
        eff_brm <-
          eff_brm |>
          group_by(.draws, Dist) |>
          summarise(Values = mean(Values)) |>
          ungroup()
        save(eff_brm, file = "../data/modelled/eff_brm_1.5_gbr_type_1.RData")

        p <- brm_effects_plot_gbr(eff_brm)
        ggsave(
          filename = "../docs/analysis_files/figure-html/contrasts_gbr_combined_disturbances_brm_1.5_type_1.png",
          p,
          width = 8, height = 8/1.6
        )
      }
      ## ----end
    }
    ## OLD Contrast - GBR Type 2. sector/shelf re_formula = NA
    {
      ## ---- Old q2_brm_contrasts_gbr_1.5_type_2
      load(file = "../data/modelled/mod_brm_1.5.RData")
      newdata <- brm_generate_newdata(data, GBR = TRUE)
      pred <- t(brms::posterior_epred(mod_brm,
        newdata = newdata,
        re_formula = NA,
        allow_new_levels = TRUE
      ))

      cellmeans_brm <- brm_generate_cellmeans(newdata, pred)
      save(cellmeans_brm, file = "../data/modelled/cellmeans_brm_1.5_gbr_type_2.RData")

      eff_brm <- brm_calc_effect(cellmeans_brm, GBR = TRUE)
      eff_brm <-
        eff_brm |>
        group_by(.draws, Dist) |>
        summarise(Values = mean(Values)) |>
        ungroup()
      save(eff_brm, file = "../data/modelled/eff_brm_1.5_gbr_type_2.RData")

      p <- brm_effects_plot_gbr(eff_brm)

      ggsave(
        filename = "../docs/analysis_files/figure-html/contrasts_gbr_combined_disturbances_brm_1.5_type_2.png",
        p,
        width = 8, height = 8 / 1.6
      )
    }
    ## OLD Contrast - GBR Type 3. sector/shelf re_formula = (1|SecShelf:(s+b+c+d+u)
    {
      ## ---- Old q2_brm_contrasts_gbr_1.5_type_3
      load(file = "../data/modelled/mod_brm_1.5.RData")
      newdata <- brm_generate_newdata(data, GBR = TRUE)
      pred <- t(brms::posterior_epred(mod_brm,
        newdata = newdata,
        re_formula = ~(1|SecShelf:(s+b+c+d+u))
      ))

      cellmeans_brm <- brm_generate_cellmeans(newdata, pred)
      save(cellmeans_brm, file = "../data/modelled/cellmeans_brm_1.5_gbr_type_3.RData")

      eff_brm <- brm_calc_effect(cellmeans_brm, GBR = TRUE)
      eff_brm <-
        eff_brm |>
        group_by(.draws, Dist) |>
        summarise(Values = mean(Values)) |>
        ungroup()
      save(eff_brm, file = "../data/modelled/eff_brm_1.5_gbr_type_3.RData")

      p <- brm_effects_plot_gbr(eff_brm)

      ggsave(
        filename = "../docs/analysis_files/figure-html/contrasts_gbr_combined_disturbances_brm_1.5_type_3.png",
        p,
        width = 8, height = 8/1.6
      )
    }
    if (1 != 2) {
      {
        ## ---- q2_brm_contrasts_gbr_combined_1.5
        eff_summ_brm <-
     cellmeans_brm |>
        mutate(Dist = Dist3) |> 
        mutate(Values = Pred) |>
        nest(.by = c(.draws)) |>
        mutate(eff = map(
          .x = data,
          .f = ~ before_vs_afters(.x)
        )) |>
        dplyr::select(-data) |>
        unnest(c(eff)) |> 
        mutate(cummulative = ifelse(str_length(Dist) > 1, TRUE, FALSE)) |>
        mutate(Dist3 = Dist) |> 
        separate_longer_position(Dist, width = 1, keep_empty = TRUE) 

        eff_summ_brm <-
     eff_summ_brm |>
        group_by(.draws, Dist) |>
        summarise(Values = mean(Values)) |>
        ungroup()


        eff_summ_brm_lab <- eff_summ_brm |>
     dplyr::select(-.draws) |>
     group_by(Dist) |>
     summarise_draws(
       median,
       HDInterval::hdi,
       P = ~ mean(. < 1)
     )
        small_dist_palette <- dist_palette |>
     filter(lab %in% c("b", "c", "d", "s", "u")) |>
     droplevels() |>
     rename(Dist = lab)

        ## eff_summ_brm |>
        ##   left_join(small_dist_palette) |>
        ##   ## left_join(eff_summ_brm_lab |> dplyr::select(Dist, P)) |>
        ##   ## mutate(color = ifelse(P > 0.9, color, "#AOAOAOAO"))
        ##   ggplot(aes(x = Values, y = Dist)) +
        ##   geom_vline(xintercept = 1, linetype = "dashed") +
        ##   ## geom_density_ridges_gradient(aes(fill = Dist, alpha = after_stat(x) < 1),
        ##   ##   show.legend = FALSE, #alpha =  0.5,
        ##   ##   rel_min_height = 0.01,
        ##   ##   size = 0.25,
        ##   ##   scale = 2,
        ##   ##   from = -10, to = 4
        ##   ##   ) +
        ##   stat_density_ridges(
        ##     ## geom = "density_ridges_gradient",
        ##     ## aes(fill = Dist, alpha = after_stat(x) < 1),
        ##     aes(fill = Dist),
        ##     ## aes(fill = after_stat(x) < 1),
        ##     alpha = 0.6,
        ##     calc_ecdf = TRUE,
        ##     quantile_lines = TRUE, quantiles = 2,
        ##     show.legend = FALSE, 
        ##     rel_min_height = 0.01,
        ##     size = 0.25,
        ##     scale = 2,
        ##     from = -10, to = 4
        ##     ) +
        ##   scale_alpha_discrete(range = c(0.2, 0.7)) +
        ##   geom_label(data = eff_summ_brm_lab,
        ##     aes(y = as.numeric(as.factor(Dist))+0.2, x = 1.1, label = sprintf("P[exceed.]==%0.3f", P)),
        ##     nudge_y =  0.1, hjust = 0, parse = TRUE) +
        ##   scale_fill_viridis_d() +
        ##     ## scale_fill_cyclical(values = small_dist_palette$color) +
        ##   ## scale_fill_manual(breaks = c(FALSE, TRUE), values = c("#A0A0A0A0", "#ff7f00A0")) +
        ##   scale_x_continuous("Effect size (fold scale decline after disturbance)",
        ##     trans = scales::log2_trans(),
        ##     labels = function(x) 1 / x,
        ##     expand = c(0, 0)
        ##   ) +
        ##   scale_y_discrete("Disturbance type",
        ##     breaks = c("b", "c", "d", "s", "u"),
        ##     labels = c("Bleaching", "COTS", "Disease", "Storms", "Unknown"),
        ##     expand = c(0, 0.1)) +
        ##   ## scale_fill_manual("Disurbances",
        ##   ##   breaks = dist_palette |> filter(lab %in% c("b","c","d","s","u")) |> pull(lab),
        ##   ##   values = dist_palette |> filter(lab %in% c("b","c","d","s","u")) |> pull(color)) +
        ##   theme_bw()

        ## my_density <- function(x) {
        ##   D <- density(log(x))
        ##  exp(as.vector(D$x[D$y > 0.1][1]))
        ## }
        eff_summ_brm_lab <- eff_summ_brm |>
     dplyr::select(-.draws) |>
     group_by(Dist) |>
     summarise_draws(
       median,
       HDInterval::hdi,
       D1 = ~ my_density(.),
       D = ~ quantile(., prob = 0.05, names = FALSE),
       P = ~ mean(. < 1)
     )
        
        p <-
     eff_summ_brm |>
        left_join(small_dist_palette) |>
        mutate(Dist = forcats::fct_relevel(Dist, "b", "c", "s", "u", "d")) |>
        ggplot(aes(x = Values, y = Dist)) +
        geom_vline(xintercept = 1, linetype = "dashed") +
        stat_slabinterval(aes(color = Dist),
          position = position_dodge(width = 0.3, preserve = "single"),
          show.legend = FALSE) +
        ## geom_pointrange(data = eff_summ_brm_lab |>
        ##                   mutate(Dist = factor(Dist, levels = c("b", "c", "s", "u", "d"))),
        ##   aes(y = as.numeric(as.factor(Dist)) - 0.1,
        ##     x = median, xmin = lower, xmax = upper,
        ##     colour = Dist),
        ##   linewidth = 1.0,
        ##   show.legend = FALSE
        ## ) +
        stat_slab() +
        ## stat_slab(aes(fill = Dist, slab_alpha = after_stat(x) < 1), show.legend = FALSE) +
        stat_slab(aes(fill = Dist, slab_alpha = after_stat(x) < 1),
          fill_type = "segments",
          height = 2,
          ## scale = .5,
          expand = FALSE, trim = TRUE, density = "bounded",
          width = 0.95,
          colour = "grey10", alpha = 0.4, show.legend = FALSE) +
        geom_text(data = eff_summ_brm_lab |>
     mutate(Dist = factor(Dist, levels = c("b", "c", "s", "u", "d"))),
     aes(y = as.numeric(as.factor(Dist))+0.2,
       x = D1, #1.5,
       label = sprintf("P[decline]==%0.3f", P)
     ),
     nudge_y =  0.2, hjust = 1, parse = TRUE) +
        ## geom_label(data = eff_summ_brm_lab |>
        ##                   mutate(Dist = factor(Dist, levels = c("b", "c", "s", "u", "d"))),
        ##   aes(y = as.numeric(as.factor(Dist))+0.2,
        ##     x = 1.5,
        ##     label = sprintf("P[exceed.]==%0.3f", P)
        ##   ),
        ##   nudge_y =  0.1, hjust = 1, parse = TRUE) +
        geom_label(
          data = eff_summ_brm_lab |>
     mutate(Dist = factor(Dist, levels = c("b", "c", "s", "u", "d"))),
     aes(
       y = as.numeric(as.factor(Dist)) + 0.4,
       x = median,
       label = sprintf("%0.2f%%", 100 * (median - 1))
     )
     ) +
        scale_fill_viridis_d() +
        scale_colour_viridis_d() +
        scale_slab_alpha_discrete(range = c(0.2, 0.7))+
        ## scale_fill_manual(breaks = c(FALSE, TRUE), values = c("#A0A0A0A0", "#ff7f00A0")) +
        ## scale_x_continuous("Effect size (fold scale decline in Seriatopora cover after disturbance)",
        ##   trans = scales::log2_trans(),
        ##   labels = function(x) 1 / x,
        ##   expand = c(0, 0),
        ##   breaks = c(1 / 128, 1 / 64, 1 / 32, 1 / 16, 1 / 8, 1 / 4, 1 / 2, 1, 2, 4, 8),
        ##   #limits = c(0.001, 10)
        ## ) +
        scale_x_continuous("Effect size (% change in Seriatopora cover after disturbance)",
          trans = scales::log2_trans(),
          labels = function(x) round(100*(x-1), 2),
          expand = c(0, 0),
          ## breaks = c(1 / 200, 1 / 50, 1 / 20, 1 / 10, 1 / 5, 0.3, 1 / 2, 1, 2, 3, 6),
          breaks = c(1 / 200, 1 / 50, 1/20, 1 / 10, 1 / 4, 1 / 2, 1, 2, 4, 10),
          #limits = c(0.001, 10)
          ) +
        scale_y_discrete("Disturbance type",
          breaks = c("b", "c", "s", "u", "d"),
          labels = c("Bleaching", "COTS", "Storms", "Unknown", "Disease"),
          expand = c(0, 0.3)) +
        ## scale_fill_manual("Disurbances",
        ##   breaks = dist_palette |> filter(lab %in% c("b","c","d","s","u")) |> pull(lab),
        ##   values = dist_palette |> filter(lab %in% c("b","c","d","s","u")) |> pull(color)) +
        theme_bw() +
        theme(axis.title.y = element_text(size = rel(1.25), margin = margin(r = 1, unit = "char")),
          axis.title.x = element_text(size = rel(1.25), margin = margin(t = 1, unit = "char")),
          axis.text.y = element_text(size = rel(1.25))) +
        coord_cartesian(clip = "on", xlim = c(0.001, 11))
        p

        ggsave(
          filename = "../docs/analysis_files/figure-html/contrasts_gbr_combined_disturbances_pretty_brm_1.5.png",
          p,
          width = 8, height = 8/1.6
        )
        
        
        ## eff_summ_brm <-
        ##   eff_summ_brm |>
        ##   dplyr::select(-.draws) |>
        ##   group_by(Dist) |>
        ##   summarise_draws(
        ##     median,
        ##     HDInterval::hdi,
        ##     P = ~ mean(. < 1)
        ##   )
        
        ## p <- 
        ## eff_summ_brm |>
        ##   droplevels() |> 
        ##   ggplot(aes(y = Dist, x = median,
        ##   alpha = P > 0.9)) +
        ##   geom_vline(xintercept = 1, linetype = "dashed") +
        ##   geom_pointrange(aes(xmin = lower, xmax = upper),
        ##     position = position_dodge(width = 0.5)) +
        ##   scale_x_continuous("Effect size (fold scale decline after disturbance(x))",
        ##     trans = scales::log2_trans(),
        ##     labels = function(x) 1/x) +
        ##   scale_alpha_manual("P > 0.9", values = c(0.3, 1)) +
        ##     scale_y_discrete("Disturbance type") +
        ##   ## facet_wrap(A_SECTOR ~ SHELF,
        ##   ##   scales = "free",
        ##   ##   labeller = label_wrap_gen(multi_line = FALSE)
        ##   ## )
        ##   theme_bw()

        ## ggsave(
        ##   filename = "../docs/analysis_files/figure-html/contrasts_gbr_combined_disturbances_brm_1.5.png",
        ##   p,
        ##   width = 6, height = 4
        ## )
      }
      ## Contrast - GBR (primary disturbances)
      {
        ## ---- q2_brm_contrasts_gbr_single_1.5
        load(file = "../data/modelled/mod_brm_1.5.RData")

        newdata <-
     data |>
        dplyr::select(Dist, s, c, d, b, u) |>
        ## filter(SecShelf == "CA M") |> 
        distinct() |>
        rowwise() |>
        mutate(Dist2 = paste(c("s", "c", "d", "b", "u")[which(c(s, c, d, b, u) == 1)], collapse = "")) |>
        ungroup()

        newdata <- bind_rows(
          ## 1. based on the Dist2 field, duplicate the rows into their components
          ## 2. create a field to indicate whether the row represents a cummulative (multiple disturbances)
          newdata |>
     mutate(Dist3 = Dist2) |> 
     separate_longer_position(Dist2, width = 1, keep_empty = TRUE) |>
     mutate(cummulative = ifelse(Dist != "Before" & str_length(Dist)>1, TRUE, FALSE)),
     newdata |>
     mutate(Dist3 = Dist2) |> 
     separate_longer_position(Dist2, width = 1, keep_empty = TRUE) |>
     mutate(s = factor(ifelse(!str_detect(Dist2, "s") | is.na(Dist2), 0, 1))) |>
     mutate(c = factor(ifelse(!str_detect(Dist2, "c") | is.na(Dist2), 0, 1))) |>
     mutate(d = factor(ifelse(!str_detect(Dist2, "d") | is.na(Dist2), 0, 1))) |>
     mutate(b = factor(ifelse(!str_detect(Dist2, "b") | is.na(Dist2), 0, 1))) |>
     mutate(u = factor(ifelse(!str_detect(Dist2, "u") | is.na(Dist2), 0, 1))) |>
     mutate(cummulative = FALSE)) |>
     mutate(Dist2 = ifelse(is.na(Dist2), "Before", Dist2)) |>
     dplyr::select(-Dist) |>
     distinct() |>
     arrange(Dist2) |>
     mutate(AIMS_REEF_NAME = NA, SecShelf = NA, Site = NA, Transect = NA, event = NA, total.points = 1) |>
     ## mutate(cummulative = ifelse(s + c + d + b + u > 1), TRUE, FALSE) |>
     ## mutate(cummulative = ifelse(str_detect(Dist, "m") | str_length(Dist) > 1, TRUE, FALSE)) |>
     filter(!is.na(s)) |>
     mutate(Dist3 = ifelse(cummulative, Dist3, Dist2)) |>
     filter(cummulative == FALSE) |>
     droplevels()
        
        p <- t(brms::posterior_epred(mod_brm, newdata = newdata, se.fit = TRUE, allow_new_levels = TRUE))

        cellmeans_brm <-
     newdata |>
        cbind(p) |> 
        pivot_longer(cols = matches("^[0-9]*$"), names_to = ".draws") |>
        mutate(Pred = (value)) |>
        mutate(
          Dist = ifelse(Dist2 == "" | is.na(Dist2), "Before", Dist2),
          Dist = forcats::fct_relevel(Dist, "Before"),
          Dist3 = forcats::fct_relevel(Dist3, "Before")
        ) 

        eff_summ_brm <-
     cellmeans_brm |>
        mutate(Dist = Dist3) |> 
        mutate(Values = Pred) |>
        nest(.by = c(.draws)) |>
        mutate(eff = map(
          .x = data,
          .f = ~ before_vs_afters(.x)
        )) |>
        dplyr::select(-data) |>
        unnest(c(eff)) |> 
        mutate(cummulative = ifelse(str_length(Dist) > 1, TRUE, FALSE)) |>
        mutate(Dist3 = Dist) |> 
        separate_longer_position(Dist, width = 1, keep_empty = TRUE) 
        
        ## eff_summ_brm <-
        ##   eff_summ_brm |>
        ##   dplyr::select(-.draws) |>
        ##   group_by(Dist, cummulative, Dist3) |>
        ##   summarise_draws(
        ##     median,
        ##     HDInterval::hdi,
        ##     P = ~ mean(. < 1)
        ##   )

        eff_summ_brm_lab <- eff_summ_brm |>
     dplyr::select(-.draws) |>
     group_by(Dist, cummulative, Dist3) |>
     summarise_draws(
       median,
       HDInterval::hdi,
       D1 = ~ my_density(.),
       D = ~ quantile(., prob = 0.05, names = FALSE),
       P = ~ mean(. < 1)
     )

        ## dist_palette <- tribble(
        ##   ~lab, ~color,
        ##   "b",    "#e41a1c",
        ##   "c",    "#377eb8",
        ##   "cb",   "#8E4C6A",
        ##   "d",    "#000000",
        ##   "s",    "#4daf4a",
        ##   "sb",   "#996533",
        ##   "sc",   "#429781",
        ##   "scb",  "#786D5F",
        ##   "sd",   "#275825",
        ##   "su",   "#a69725",
        ##   "u",    "#ff7f00",
        ## )

        ## p <- 
        ## eff_summ_brm |>
        ##   droplevels() |> 
        ##   ggplot(aes(y = Dist, x = median,
        ##   color = Dist3, #shape = cummulative,
        ##   size = factor(str_length(Dist3)),
        ##   alpha = P > 0.9)) +
        ##   geom_vline(xintercept = 1, linetype = "dashed") +
        ##   geom_pointrange(aes(xmin = lower, xmax = upper),
        ##     position = position_dodge(width = 0.5)) +
        ##     scale_shape_manual(values = c(16, 21)) +
        ##     scale_size_manual(str_wrap("Number of disturbances types", 15), values = c(0.25,0.5,1)) +
        ##     ## scale_colour_discrete("Disurbances") +
        ##   scale_colour_manual("Disurbances",
        ##     breaks = dist_palette$lab,
        ##     values = dist_palette$color) +
        ##   scale_x_continuous("Effect size (fold scale decline after disturbance(x))",
        ##     trans = scales::log2_trans(),
        ##     labels = function(x) 1/x) +
        ##   scale_alpha_manual("P > 0.9", values = c(0.3, 1)) +
        ##     scale_y_discrete("Disturbance type") +
        ##   ## facet_wrap(A_SECTOR ~ SHELF,
        ##   ##   scales = "free",
        ##   ##   labeller = label_wrap_gen(multi_line = FALSE)
        ##   ## )
        ##   theme_bw()
        ## p

        p <-
     eff_summ_brm |>
        left_join(small_dist_palette) |>
        mutate(Dist = forcats::fct_relevel(Dist, "b", "c", "s", "u", "d")) |>
        ggplot(aes(x = Values, y = Dist)) +
        geom_vline(xintercept = 1, linetype = "dashed") +
        stat_slabinterval(aes(color = Dist),
          position = position_dodge(width = 0.3, preserve = "single"),
          show.legend = FALSE) +
        ## geom_pointrange(data = eff_summ_brm_lab |>
        ##                   mutate(Dist = factor(Dist, levels = c("b", "c", "s", "u", "d"))),
        ##   aes(y = as.numeric(as.factor(Dist)) - 0.1,
        ##     x = median, xmin = lower, xmax = upper,
        ##     colour = Dist),
        ##   linewidth = 1.0,
        ##   show.legend = FALSE
        ## ) +
        stat_slab() +
        ## stat_slab(aes(fill = Dist, slab_alpha = after_stat(x) < 1), show.legend = FALSE) +
        stat_slab(aes(fill = Dist, slab_alpha = after_stat(x) < 1),
          fill_type = "segments",
          height = 2,
          ## scale = .5,
          expand = FALSE, trim = TRUE, density = "bounded",
          width = 0.95,
          colour = "grey10", alpha = 0.4, show.legend = FALSE) +
        geom_text(data = eff_summ_brm_lab |>
     mutate(Dist = factor(Dist, levels = c("b", "c", "s", "u", "d"))),
     aes(y = as.numeric(as.factor(Dist))+0.2,
       x = D1, #1.5,
       label = sprintf("P[decline]==%0.3f", P)
     ),
     nudge_y =  0.2, hjust = 1, parse = TRUE) +
        ## geom_label(data = eff_summ_brm_lab |>
        ##                   mutate(Dist = factor(Dist, levels = c("b", "c", "s", "u", "d"))),
        ##   aes(y = as.numeric(as.factor(Dist))+0.2,
        ##     x = 1.5,
        ##     label = sprintf("P[exceed.]==%0.3f", P)
        ##   ),
        ##   nudge_y =  0.1, hjust = 1, parse = TRUE) +
        geom_label(
          data = eff_summ_brm_lab |>
     mutate(Dist = factor(Dist, levels = c("b", "c", "s", "u", "d"))),
     aes(
       y = as.numeric(as.factor(Dist)) + 0.2,
       x = median,
       label = sprintf("%0.2f%%", 100 * (median - 1))
     )
     ) +
        scale_fill_viridis_d() +
        scale_colour_viridis_d() +
        scale_slab_alpha_discrete(range = c(0.2, 0.7))+
        ## scale_fill_manual(breaks = c(FALSE, TRUE), values = c("#A0A0A0A0", "#ff7f00A0")) +
        ## scale_x_continuous("Effect size (fold scale decline in Seriatopora cover after disturbance)",
        ##   trans = scales::log2_trans(),
        ##   labels = function(x) 1 / x,
        ##   expand = c(0, 0),
        ##   breaks = c(1 / 128, 1 / 64, 1 / 32, 1 / 16, 1 / 8, 1 / 4, 1 / 2, 1, 2, 4, 8),
        ##   #limits = c(0.001, 10)
        ## ) +
        scale_x_continuous("Effect size (% change in Seriatopora cover after disturbance)",
          trans = scales::log2_trans(),
          labels = function(x) round(100*(x-1), 2),
          expand = c(0, 0),
          ## breaks = c(1 / 200, 1 / 50, 1 / 20, 1 / 10, 1 / 5, 0.3, 1 / 2, 1, 2, 3, 6),
          breaks = c(1 / 200, 1 / 50, 1/20, 1 / 10, 1 / 4, 1 / 2, 1, 2, 4, 10),
          #limits = c(0.001, 10)
          ) +
        scale_y_discrete("Disturbance type",
          breaks = c("b", "c", "s", "u", "d"),
          labels = c("Bleaching", "COTS", "Storms", "Unknown", "Disease"),
          expand = c(0, 0.3)) +
        ## scale_fill_manual("Disurbances",
        ##   breaks = dist_palette |> filter(lab %in% c("b","c","d","s","u")) |> pull(lab),
        ##   values = dist_palette |> filter(lab %in% c("b","c","d","s","u")) |> pull(color)) +
        theme_bw() +
        theme(axis.title.y = element_text(size = rel(1.25), margin = margin(r = 1, unit = "char")),
          axis.title.x = element_text(size = rel(1.25), margin = margin(t = 1, unit = "char")),
          axis.text.y = element_text(size = rel(1.25))) +
        coord_cartesian(clip = "on", xlim = c(0.001, 11))
        p
        ggsave(
          filename = "../docs/analysis_files/figure-html/contrasts_gbr_primary_disturbances_pretty_brm_1.5.png",
          p,
          width = 8, height = 8/1.6
        )

        
      }
    }
  }
}
