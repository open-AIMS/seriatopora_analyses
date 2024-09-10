## Retrieve data
## ---- q3_retrieve_data_1
load(file = "../data/processed/data_q3.RData")
full_data <- data
## ----end


## Abundance (density)

## add REPORT_YEAR and Sector/Shelf
{
  ## Prepare the data
  {
    ## ---- q3_prepare_data_for_model
    data <- data |>
      dplyr::select(
        AIMS_REEF_NAME, zone_depth, Site, abund,
        cREPORT_YEAR, A_SECTOR, SHELF, area.avail
      ) |>
      mutate(
        SecShelfYr = paste(A_SECTOR, SHELF, cREPORT_YEAR),
        SecShelf = paste(A_SECTOR, SHELF)
      )
    ## what combinations are missing
    data |>
      group_by(SecShelfYr) |>
      summarise(Abund = sum(abund)) |>
      arrange(desc(Abund))
    ## 492 missing - lets exclude those
    data_pre_excludes <- data
    excludes <- data |>
      group_by(SecShelfYr) |>
      summarise(Abund = sum(abund)) |>
      filter(Abund == 0) |>
      mutate(EXCLUDE = TRUE)
    save(excludes, file = "../data/modelled/excludes_q3.RData")
    data <- data |>
      left_join(excludes) |>
      filter(is.na(EXCLUDE)) |>
      droplevels() |>
      mutate(
        SecShelfYr = forcats::fct_relevel(SecShelfYr, "IN M 2016")
        ## SecShelfYr = forcats::fct_relevel(SecShelfYr, "SW M 2014")
      ) |>
      dplyr::select(-Abund, -EXCLUDE)
    ## ----end
  }
  ## glmmTMB
  {
    ## Fit model
    {
      if (rerun_models) {
        ## ---- q3_glmmTMB_fit
        mod_glmmTMB <- glmmTMB(
          abund ~ 1 + offset(log(area.avail)) +
            (1 | SecShelfYr) +
            (1 | AIMS_REEF_NAME) +
            (1 | zone_depth) +
            (1 | Site),
          ziformula = ~ (1|SecShelfYr) + (1| zone_depth) + (1 | Site),
          ## ziformula = ~ (1| zone_depth) + (1 | Site),
          data = data,
          ## family = "poisson",
          family = "nbinom2",
          REML = TRUE
        )

        save(mod_glmmTMB, file = "../data/modelled/mod_q3_glmmTMB.RData")
        ## ----end
      }
    }
    ## DHARMa
    {
      ## ---- q3_glmmTMB_validation
      load(file = "../data/modelled/mod_q3_glmmTMB.RData")
      resids <- simulateResiduals(mod_glmmTMB, plot = FALSE)
      wrap_elements(~ testUniformity(resids)) +
        wrap_elements(~ plotResiduals(resids)) +
        wrap_elements(~ testDispersion(resids))
      ## ----end
      testDispersion(resids)
      testZeroInflation(resids)
    }
    ## Summaries
    {
      ## ---- q3_glmmTMB_summaries_print
      load(file = "../data/modelled/mod_q3_glmmTMB.RData")
      summ <- bind_rows(
        model_summary(mod_glmmTMB) |> mutate(model = "Raw"),
        model_summary(mod_glmmTMB, exponentiate = TRUE) |> mutate(model = "Exponentiated")
      ) |>
        mutate(across(c(estimate, conf.low, conf.high), \(x) sprintf("%.3f", x))) |>
        mutate(estimate = paste0(estimate, "<br>[", conf.low, ",", conf.high, "]")) |>
        mutate(estimate = str_replace(estimate, "<br>\\[NA,NA\\]", "")) |>
        dplyr::select(-conf.low, -conf.high) |>
        pivot_longer(
          cols = c("estimate"),
          names_to = "names", values_to = "values"
        ) |>
        pivot_wider(id_cols = everything(), names_from = "model", values_from = "values") |>
        dplyr::select(-effect, -names)

      summ_tt <- summ |>
        tinytable::tt()
      save(summ, file = "../data/summ_glmmTMB_q3.RData")
      save(summ_tt, file = "../data/summ_tt_glmmTMB_q3.RData")
      ## ----end
      ## ---- q3_glmmTMB_summaries
      load(file = "../data/summ_glmmTMB_q3.RData")
      summ_tt
      ## ----end
    }
    ## Partial plots
    {
      ## ---- q3_glmmTMB_partial_effects
      load(file = "../data/modelled/mod_q3_glmmTMB.RData")
      load(file = "../data/modelled/excludes_q3.RData")
      excludes <- excludes |>
        mutate(rate = 0, asymp.LCL = 0, asymp.UCL = 0) |>
        dplyr::select(-Abund, -EXCLUDE) |>
        separate(SecShelfYr, into = c("A_SECTOR", "SHELF", "cREPORT_YEAR"), remove = FALSE)
      newdata <- brm_generate_newdata_no_dist(data) |>
        mutate(area.avail = 1) |>
        filter(str_detect(SecShelfYr, "2024", negate = TRUE))
      pred <- predict(mod_glmmTMB,
          newdata = newdata,
          re.form = NULL, allow.new.levels = TRUE,
          se.fit = TRUE
        )
      mod_em <-
        newdata |>
        mutate(asymp.LCL = exp(pred$fit - 2* pred$se.fit),
          asymp.UCL = exp(pred$fit + 2* pred$se.fit),
          rate = exp(pred$fit)) |>
        as.data.frame() |>
        separate(SecShelfYr, into = c("A_SECTOR", "SHELF", "cREPORT_YEAR"), remove = FALSE) |>
        full_join(excludes) |>
        mutate(A_SECTOR = factor(A_SECTOR, levels = c("CG", "PC", "CL", "CA", "IN", "TO", "WH", "PO", "SW", "CB")))
      p <- mod_em |>
        ggplot(aes(y = rate, x = as.numeric(as.character(cREPORT_YEAR)))) +
        geom_line() +
        geom_pointrange(aes(ymin = asymp.LCL, ymax = asymp.UCL)) +
        facet_grid(A_SECTOR ~ SHELF, scales = "free")
      ggsave(filename = "../docs/analysis_q3_files/figure-html/partial_glmmTMB_q3.png", width = 10, height = 10)
      ## ----end
      ## ---- q3_glmmTMB_partial_effects_old
      load(file = "../data/modelled/mod_q3_glmmTMB.RData")
      load(file = "../data/modelled/excludes_q3.RData")
      excludes <- excludes |>
        mutate(rate = 0, asymp.LCL = 0, asymp.UCL = 0) |>
        dplyr::select(-Abund, -EXCLUDE) |>
        separate(SecShelfYr, into = c("A_SECTOR", "SHELF", "cREPORT_YEAR"), remove = FALSE)
      mod_em <-
        mod_glmmTMB |>
        emmeans(~SecShelfYr, type = "response") |>
        as.data.frame() |>
        separate(SecShelfYr, into = c("A_SECTOR", "SHELF", "cREPORT_YEAR"), remove = FALSE) |>
        full_join(excludes) |>
        mutate(A_SECTOR = factor(A_SECTOR, levels = c("CG", "PC", "CL", "CA", "IN", "TO", "WH", "PO", "SW", "CB")))
      p <- mod_em |>
        ggplot(aes(y = rate, x = as.numeric(as.character(cREPORT_YEAR)))) +
        geom_line() +
        geom_pointrange(aes(ymin = asymp.LCL, ymax = asymp.UCL)) +
        facet_grid(A_SECTOR ~ SHELF, scales = "free")
      ggsave(filename = "../docs/analysis_q3_files/figure-html/partial_glmmTMB_q3.png", width = 10, height = 10)
      ## ----end
    }
  }
  ## brms
  {
    ## Fit model
    {
      if (rerun_models) {
        ## ---- q3_brm_fit
        form <- bf(
          abund ~ 1 + offset(log(area.avail)) +
            (1 | SecShelfYr) +
            (1 | AIMS_REEF_NAME) +
            (1 | zone_depth) +
            (1 | Site),
          zi ~ (1|SecShelfYr) + (1| zone_depth) + (1 | Site),
          family = zero_inflated_poisson()
        )
        data |>
          group_by(SecShelfYr) |>
          summarise(
            Mean = log(mean(abund/area.avail)),
            SD = log(sd(abund/area.avail))
          )
        priors <- prior(normal(0.5, 1), class = "Intercept") +
          ## prior(normal(0, 2), class = "b") +
          prior(student_t(3, 0, 1), class = "sd") +
          prior(logistic(0, 1), class = "Intercept", dpar = "zi") +
          ## prior(normal(0,1), class = 'b', dpar = 'zi') +
          prior(student_t(3, 0, 1), class = "sd", dpar = "zi") 
        mod_brm <- brm(
          formula = form,
          data = data,
          prior = priors,
          iter = 5000, warmup = 1000,
          chains = 3, cores = 3,
          sample_prior = "yes",
          thin = 10,
          backend = "cmdstanr",
          control = list(adapt_delta = 0.99)
        )
        save(mod_brm, file = "../data/modelled/mod_q3_brm_1.RData")
        ## ----end
      }
    }
    
    ## MCMC diagnostics
    {
      ## ---- q3_brm_trace_1
      load(file = "../data/modelled/mod_q3_brm_1.RData")
      rstan::stan_trace(mod_brm$fit)
      rstan::stan_ac(mod_brm$fit)
      rstan::stan_ess(mod_brm$fit) + rstan::stan_rhat(mod_brm$fit)
      ## ----end
      
      ## ---- q3_brm_trace
      load(file = "../data/modelled/mod_q3_brm_1.RData")
      ggsave(
        filename = "../docs/analysis_files/figure-html/stan_trace_q3_brm_1.png",
        rstan::stan_trace(mod_brm$fit),
        width = 15, height = 8
      )
      ggsave(
        filename = "../docs/analysis_files/figure-html/stan_ac_q3_brm_1.png",
        rstan::stan_ac(mod_brm$fit),
        width = 15, height = 8
      )
      ggsave(
        filename = "../docs/analysis_files/figure-html/stan_ess_q3_brm_1.png",
        rstan::stan_ess(mod_brm$fit) + rstan::stan_rhat(mod_brm$fit),
        width = 15, height = 8
      )
      ## ----end
    }
    ## Model validation
    {
      ## ---- q3_brm_validation_a
      load(file = "../data/modelled/mod_q3_brm_1.RData")
      mod_brm |> pp_check(type = "dens_overlay", ndraws = 100) + scale_x_log10() +
        mod_brm |> pp_check(type = "loo_pit_overlay", ndraws = 100)
      ## ----end

      ## ---- q3_brm_validation
      load(file = "../data/modelled/mod_q3_brm_1.RData")
      ggsave(
        filename = "../docs/analysis_files/figure-html/pcc_q3_brm_1.png",
        mod_brm |> pp_check(type = "dens_overlay", ndraws = 100) + scale_x_log10() +
          mod_brm |> pp_check(type = "loo_pit_overlay", ndraws = 100),
        width = 10,
        height = 5)
      ## ----end
    }
    ## DHARMa residuals
    {
      ## ---- q3_brm_dharma_a
      load(file = "../data/modelled/mod_q3_brm_1.RData")
      resids <- make_brms_dharma_res(
        mod_brm,
        integerResponse = TRUE
      )
      wrap_elements(~ testUniformity(resids)) +
        wrap_elements(~ plotResiduals(resids, form = factor(rep(1, nrow(mod_brm$data))))) +
        wrap_elements(~ plotResiduals(resids)) +
        wrap_elements(~ testDispersion(resids))
      ## ----end

      ## ---- q3_brm_dharma
      load(file = "../data/modelled/mod_q3_brm_1.RData")
      resids <- make_brms_dharma_res(
        mod_brm,
        integerResponse = FALSE
      )
      ggsave(
        filename = "../docs/analysis_files/figure-html/dharma_q3_brm_1.png",
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
      ## ---- q3_brm_summary
      load(file = "../data/modelled/mod_q3_brm_1.RData")

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
      save(summ, file = "../data/summ_q3_brm_2_1.RData")
      save(summ_tt, file = "../data/summ_tt_q3_brm_1.RData")
      ## ----end
      ## ---- q3_brm_summary_print
      load(file = "../data/summ_tt_q3_brm_1.RData")
      summ_tt
      ## ----end
    }
    ## Partial plots -  Type 1. sector/shelf re_formula = NULL
    {
      ## ---- q3_brm_partial_effects_1
      load(file = "../data/modelled/mod_q3_brm_1.RData")
      load(file = "../data/modelled/excludes_q3.RData")
      excludes <- excludes |>
        mutate(median = 0, lower = 0, upper = 0) |>
        dplyr::select(-Abund, -EXCLUDE) |>
        separate(SecShelfYr, into = c("A_SECTOR", "SHELF", "cREPORT_YEAR"), remove = FALSE)
      
      newdata <- brm_generate_newdata_no_dist(data) |>
        mutate(area.avail = 1)
      pred <- t(brms::posterior_epred(mod_brm,
        newdata = newdata,
        re_formula = NULL,
        se.fit = TRUE,
        allow_new_levels = TRUE
      ))
      cellmeans_brm <- brm_generate_cellmeans_no_dist(newdata, pred)
      save(cellmeans_brm, file = "../data/modelled/cellmeans_q3_brm_type_1_1.RData")

      cellmeans_summ_brm <- cellmeans_brm |>
        dplyr::select(-.draws) |>
        dplyr::select(A_SECTOR, SHELF, REPORT_YEAR, Pred) |>
        group_by(A_SECTOR, SHELF, REPORT_YEAR) |>
        summarise_draws(median, HDInterval::hdi) |>
        dplyr::select(A_SECTOR, SHELF, REPORT_YEAR, median, lower, upper) |> 
        full_join(excludes) |>
        mutate(A_SECTOR = factor(A_SECTOR, levels = c("CG", "PC", "CL", "CA", "IN", "TO", "WH", "PO", "SW", "CB")))


      save(cellmeans_summ_brm, file = "../data/modelled/cellmeans_summ_q3_brm_type_1_1.RData")

      p <- brm_partial_plot_no_dist(cellmeans_summ_brm)
      ggsave(
        filename = "../docs/analysis_files/figure-html/partial_q3_brm_type_1_1.png",
        p,
        width = 15, height = 12
      )


      
      ## cellmeans_summ_brm <- cellmeans_brm |>
      ##   dplyr::select(-.draws) |>
      ##   dplyr::select(A_SECTOR, SHELF, REPORT_YEAR, Pred) |>
      ##   group_by(A_SECTOR, SHELF, REPORT_YEAR) |>
      ##   summarise_draws(median, HDInterval::hdi) |>
      ##   dplyr::select(A_SECTOR, SHELF, REPORT_YEAR, median, lower, upper)

      ## cellmeans_brm <- newdata |>
      ##   mutate(pred <- t(brms::posterior_epred(mod_brm,
      ##   newdata = newdata,
      ##   re_formula = NULL,
      ##   se.fit = TRUE
      ##   ## allow_new_levels = TRUE
      ## ))
      ## ## cellmeans_brm <- brm_generate_cellmeans_no_dist(newdata, pred)
      ## cellmeans_brm <-
      ##   emmeans(mod_brm, ~SecShelfYr, type = "response") |>
      ##   gather_emmeans_draws() |>
      ##   separate(SecShelfYr, into = c("A_SECTOR", "SHELF", "REPORT_YEAR"), remove = TRUE) |>
      ##   mutate(Pred = exp(.value))

      ## save(cellmeans_brm, file = "../data/modelled/cellmeans_q3_brm_type_1.RData")

      ## cellmeans_summ_brm <- cellmeans_brm |>
      ##   dplyr::select(A_SECTOR, SHELF, REPORT_YEAR, Pred) |>
      ##   group_by(A_SECTOR, SHELF, REPORT_YEAR) |>
      ##   summarise_draws(median, HDInterval::hdi) |>
      ##   dplyr::select(A_SECTOR, SHELF, REPORT_YEAR, median, lower, upper)

      ## save(cellmeans_summ_brm, file = "../data/modelled/cellmeans_summ_q3_brm_type_1.RData")

      ## p <- brm_partial_plot_no_dist(cellmeans_summ_brm)
      ## ggsave(
      ##   filename = "../docs/analysis_files/figure-html/partial_q3_brm_type_1.png",
      ##   p,
      ##   width = 15, height = 12
      ## )
      ## ----end
    }
    ## Partial plots -  Type 2. sector/shelf re_formula = NA
    {
      ## ---- q3_brm_partial_effects_type_2
      load(file = "../data/modelled/mod_q3_brm_1.RData")
      load(file = "../data/modelled/excludes_q3.RData")
      excludes <- excludes |>
        mutate(median = 0, lower = 0, upper = 0) |>
        dplyr::select(-Abund, -EXCLUDE) |>
        separate(SecShelfYr, into = c("A_SECTOR", "SHELF", "cREPORT_YEAR"), remove = FALSE)
      newdata <- brm_generate_newdata_no_dist(data) |> 
        mutate(area.avail = 1)
      pred <- t(brms::posterior_epred(mod_brm,
        newdata = newdata,
        re_formula = ~ (1 | SecShelfYr),
        se.fit = TRUE,
        allow_new_levels = TRUE
      ))
      cellmeans_brm <- brm_generate_cellmeans_no_dist(newdata, pred)

      save(cellmeans_brm, file = "../data/modelled/cellmeans_q3_brm_type_2_1.RData")

      cellmeans_summ_brm <- cellmeans_brm |>
        dplyr::select(-.draws) |>
        dplyr::select(A_SECTOR, SHELF, REPORT_YEAR, Pred) |>
        group_by(A_SECTOR, SHELF, REPORT_YEAR) |>
        summarise_draws(median, HDInterval::hdi) |>
        dplyr::select(A_SECTOR, SHELF, REPORT_YEAR, median, lower, upper) |> 
        full_join(excludes) |>
        mutate(A_SECTOR = factor(A_SECTOR, levels = c("CG", "PC", "CL", "CA", "IN", "TO", "WH", "PO", "SW", "CB")))

      
      save(cellmeans_summ_brm, file = "../data/modelled/cellmeans_summ_q3_brm_type_2_1.RData")

      p <- brm_partial_plot_no_dist(cellmeans_summ_brm)
      p
      ggsave(
        filename = "../docs/analysis_files/figure-html/partial_q3_brm_type_2_1.png",
        p,
        width = 15, height = 12
      )
      ## load(file = "../data/modelled/mod_q3_brm.RData")
      ## newdata <- brm_generate_newdata_no_dist(data)
      ## pred <- t(brms::posterior_epred(mod_brm,
      ##   newdata = newdata,
      ##   re_formula = NA,
      ##   se.fit = TRUE,
      ##   allow_new_levels = TRUE
      ## ))
      ## cellmeans_brm <- brm_generate_cellmeans_no_dist(newdata, pred)

      ## save(cellmeans_brm, file = "../data/modelled/cellmeans_q3_brm_type_2.RData")

      ## cellmeans_summ_brm <- cellmeans_brm |>
      ##   dplyr::select(-.draws) |>
      ##   dplyr::select(A_SECTOR, SHELF, REPORT_YEAR, Pred) |>
      ##   group_by(A_SECTOR, SHELF, REPORT_YEAR) |>
      ##   summarise_draws(median, HDInterval::hdi) |>
      ##   dplyr::select(A_SECTOR, SHELF, REPORT_YEAR, median, lower, upper)

      ## save(cellmeans_summ_brm, file = "../data/modelled/cellmeans_summ_q3_brm_type_2.RData")

      ## p <- brm_partial_plot_no_dist(cellmeans_summ_brm)
      ## ggsave(
      ##   filename = "../docs/analysis_files/figure-html/partial_q3_brm_type_2.png",
      ##   p,
      ##   width = 15, height = 12
      ## )
      ## ----end
    }
  }

}
