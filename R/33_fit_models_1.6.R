## Retrieve the data
## ---- q2_data_1.6
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
  ## ---- q2_prepare_data_1.6
  ## Focus on only the necessary variables
  data <- data |>
    mutate(SecShelf = paste(A_SECTOR, SHELF)) |>
    mutate(zone_depth = factor(paste0(REEF_ZONE, DEPTH))) |>
    dplyr::select(
      n.points, total.points, Dist.time, s, c, d, b, u,
      AIMS_REEF_NAME, Site, Transect, SecShelf, Dist.number, Dist,
      SecShelf, SSDist, event, zone_depth
    ) |>
    ## filter(!SecShelf %in% c("CG M", "PC M", "PC O")) |>
    droplevels()
  ## mutate(across(c(s, c, d, b, u), \(x) ifelse(Dist.time == "Before", 0, x))) 
  ## ----end
}


## Fit the models
{
  ## Raw means
  {
    ## ---- q2_raw_fitted_1.6
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
    save(cellmeans_summ_raw, file = "../data/modelled/cellmeans_summ_raw_1.6.RData")
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
      filename = "../docs/analysis_files/figure-html/partial_raw_1.6.png",
      p,
      width = 15, height = 10
    )
    ## ----end
  }
  ## brms
  {
    ## fit model
    {
      if (rerun_models) {
        ## ---- q2_brm_fitted_1.6
        form <- bf(
          n.points | trials(total.points) ~ SecShelf * (s + c + d + b + u) +
            (1 | event) +
            (1 | AIMS_REEF_NAME) +
            (1 | zone_depth) +
            (1 | Site) +
            (1 | Transect),
          family = "beta_binomial"
        )
        data |>
          group_by(SecShelf, Dist) |> 
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
          chains = 3, cores = 3,
          sample_prior = "yes",
          thin = 2,
          backend = "cmdstanr",
          control = list(adapt_delta = 0.99)
        )
        summary(mod_brm)
        save(mod_brm, file = "../data/modelled/mod_brm_1.6.RData")
        ## ----end
      }
    }
    ## MCMC diagnostics
    {
      ## ---- q2_brm_trace_1.6a
      load(file = "../data/modelled/mod_brm_1.6.RData")
      rstan::stan_trace(mod_brm$fit)
      rstan::stan_ac(mod_brm$fit)
      rstan::stan_ess(mod_brm$fit) + rstan::stan_rhat(mod_brm$fit)
      ## ----end
      
      ## ---- q2_brm_trace_1.6
      load(file = "../data/modelled/mod_brm_1.6.RData")
      ggsave(
        filename = "../docs/analysis_files/figure-html/stan_trace_brm_1.6.png",
        rstan::stan_trace(mod_brm$fit),
        width = 15, height = 8
      )
      ggsave(
        filename = "../docs/analysis_files/figure-html/stan_ac_brm_1.6.png",
        rstan::stan_ac(mod_brm$fit),
        width = 15, height = 8
      )
      ggsave(
        filename = "../docs/analysis_files/figure-html/stan_ess_brm_1.6.png",
        rstan::stan_ess(mod_brm$fit) + rstan::stan_rhat(mod_brm$fit),
        width = 15, height = 8
      )
      ## ----end
    }
    ## Model validation
    {
      ## ---- q2_brm_validation_1.6a
      load(file = "../data/modelled/mod_brm_1.6.RData")
      mod_brm |> pp_check(type = "dens_overlay", ndraws = 100) + scale_x_log10() +
        mod_brm |> pp_check(type = "loo_pit_overlay", ndraws = 100)
      ## ----end

      ## ---- q2_brm_validation_1.6
      load(file = "../data/modelled/mod_brm_1.6.RData")
      ggsave(
        filename = "../docs/analysis_files/figure-html/pcc_brm_1.6.png",
        mod_brm |> pp_check(type = "dens_overlay", ndraws = 100) + scale_x_log10() +
          mod_brm |> pp_check(type = "loo_pit_overlay", ndraws = 100),
        width = 10,
        height = 5)
      ## ----end
    }
    ## DHARMa residuals
    {
      ## ---- q2_brm_dharma_1.6a
      load(file = "../data/modelled/mod_brm_1.6.RData")
      resids <- make_brms_dharma_res(
        mod_brm,
        integerResponse = FALSE
      )
      wrap_elements(~ testUniformity(resids)) +
        wrap_elements(~ plotResiduals(resids, form = factor(rep(1, nrow(mod_brm$data))))) +
        wrap_elements(~ plotResiduals(resids)) +
        wrap_elements(~ testDispersion(resids))
      ## ----end

      ## ---- q2_brm_dharma_1.6
      load(file = "../data/modelled/mod_brm_1.6.RData")
      resids <- make_brms_dharma_res(
        mod_brm,
        integerResponse = FALSE
      )
      ggsave(
        filename = "../docs/analysis_files/figure-html/dharma_brm_1.6.png",
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
      ## ---- q2_brm_summary_1.6
      load(file = "../data/modelled/mod_brm_1.6.RData")

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
      save(summ, file = "../data/summ_brm_1.6.RData")
      save(summ_tt, file = "../data/summ_tt_brm_1.6.RData")
      ## ----end
      ## ---- q2_brm_summary_1.6_print
      load(file = "../data/summ_tt_brm_1.6.RData")
      summ_tt
      ## ----end
    }
    ## Partial plots - Type 1. sector/shelf re_formula = NULL
    {
      ## ---- q2_brm_partial_effects_1.6_type_1
      load(file = "../data/modelled/mod_brm_1.6.RData")

      newdata <- brm_generate_newdata(data)
      pred <- t(brms::posterior_epred(mod_brm,
        newdata = newdata,
        re_formula = NULL
      ))
      cellmeans_brm <- brm_generate_cellmeans(newdata, pred)

      save(cellmeans_brm, file = "../data/modelled/cellmeans_brm_1.6_type_1.RData")

      cellmeans_summ_brm <- cellmeans_brm |>
        dplyr::select(-.draws) |>
        dplyr::select(A_SECTOR, SHELF, Dist, Dist2, Dist3, cummulative, Pred) |>
        group_by(A_SECTOR, SHELF, Dist, Dist2, Dist3, cummulative) |>
        summarise_draws(median, HDInterval::hdi) |>
        dplyr::select(Dist, median, lower, upper)

      save(cellmeans_summ_brm, file = "../data/modelled/cellmeans_summ_brm_1.6_type_1.RData")

      p <- brm_partial_plot(cellmeans_summ_brm)
      ggsave(
        filename = "../docs/analysis_files/figure-html/partial_brm_1.6_type_1.png",
        p,
        width = 15, height = 12
      )
      ## ----end
    }
    ## Partial plots - Type 2. sector/shelf re_formula = NA
    {
      ## ---- q2_brm_partial_effects_1.6_type_2
      load(file = "../data/modelled/mod_brm_1.6.RData")

      newdata <- brm_generate_newdata(data)
      pred <- t(brms::posterior_epred(mod_brm,
        newdata = newdata,
        re_formula = NA
      ))
      cellmeans_brm <- brm_generate_cellmeans(newdata, pred)

      save(cellmeans_brm, file = "../data/modelled/cellmeans_brm_1.6_type_2.RData")

      cellmeans_summ_brm <- cellmeans_brm |>
        dplyr::select(-.draws) |>
        dplyr::select(A_SECTOR, SHELF, Dist, Dist2, Dist3, cummulative, Pred) |>
        group_by(A_SECTOR, SHELF, Dist, Dist2, Dist3, cummulative) |>
        summarise_draws(median, HDInterval::hdi) |>
        dplyr::select(Dist, median, lower, upper)

      save(cellmeans_summ_brm, file = "../data/modelled/cellmeans_summ_brm_1.6_type_2.RData")

      p <- brm_partial_plot(cellmeans_summ_brm)
      ggsave(
        filename = "../docs/analysis_files/figure-html/partial_brm_1.6_type_2.png",
        p,
        width = 15, height = 12
      )
      ## ----end
    }

    ## Contrasts - Type 1. sector/shelf re_formula = NULL
    {
      ## ---- q2_brm_contrasts_1.6_type_1
      load(file = "../data/modelled/cellmeans_brm_1.6_type_1.RData")

      eff_brm <- brm_calc_effect(cellmeans_brm)
      save(eff_brm, file = "../data/modelled/eff_brm_1.6_type_1.RData")
      load(file = "../data/modelled/eff_brm_1.6_type_1.RData")

      p <- brm_effects_plot(eff_brm)
      ggsave(
        filename = "../docs/analysis_files/figure-html/contrasts_brm_1.6_type_1.png",
        p,
        width = 15, height = 12
      )
    }
    ## Contrasts - Type 2. sector/shelf re_formula = NA
    {
      ## ---- q2_brm_contrasts_1.6_type_1
      load(file = "../data/modelled/cellmeans_brm_1.6_type_2.RData")

      eff_brm <- brm_calc_effect(cellmeans_brm)
      save(eff_brm, file = "../data/modelled/eff_brm_1.6_type_2.RData")
      load(file = "../data/modelled/eff_brm_1.6_type_2.RData")

      p <- brm_effects_plot(eff_brm)
      ggsave(
        filename = "../docs/analysis_files/figure-html/contrasts_brm_1.6_type_2.png",
        p,
        width = 15, height = 12
      )
    }

    ## Contrast - GBR Type 1. sector/shelf re_formula = NULL
    {
      ## ---- q2_brm_contrasts_gbr_1.6_type_1
      load(file = "../data/modelled/mod_brm_1.6.RData")
      load(file = "../data/modelled/cellmeans_brm_1.6_type_1.RData")

      cellmeans_brm <- cellmeans_brm |>
        group_by(A_SECTOR, SHELF, Dist, .draws) |>
        summarise(value = mean(value)) |>
        ungroup() 

      save(cellmeans_brm, file = "../data/modelled/cellmeans_brm_1.6_gbr_type_1.RData")
      
      eff_brm <- brm_calc_effect_hier(cellmeans_brm)

      eff_brm <-
        eff_brm |>
        group_by(.draws, Dist) |>
        summarise(Values = mean(Values)) |>
        ungroup()
      save(eff_brm, file = "../data/modelled/eff_brm_1.6_gbr_type_1.RData")
      load(file = "../data/modelled/eff_brm_1.6_gbr_type_1.RData")

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
        filename = "../docs/analysis_files/figure-html/contrasts_gbr_combined_disturbances_brm_1.6_type_1.png",
        p,
        width = 8, height = 8/1.6
      )

    }
    ## Contrast - GBR Type 2. sector/shelf re_formula = NA
    {
      ## ---- q2_brm_contrasts_gbr_1.6_type_2
      load(file = "../data/modelled/mod_brm_1.6.RData")
      load(file = "../data/modelled/cellmeans_brm_1.6_type_2.RData")

      cellmeans_brm <- cellmeans_brm |>
        group_by(A_SECTOR, SHELF, Dist, .draws) |>
        summarise(value = mean(value)) |>
        ungroup() 

      save(cellmeans_brm, file = "../data/modelled/cellmeans_brm_1.6_gbr_type_2.RData")

      eff_brm <- brm_calc_effect_hier(cellmeans_brm)

      eff_brm <-
        eff_brm |>
        group_by(.draws, Dist) |>
        summarise(Values = mean(Values)) |>
        ungroup()
      save(eff_brm, file = "../data/modelled/eff_brm_1.6_gbr_type_2.RData")
      load(file = "../data/modelled/eff_brm_1.6_gbr_type_2.RData")

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
        filename = "../docs/analysis_files/figure-html/contrasts_gbr_combined_disturbances_brm_1.6_type_2.png",
        p,
        width = 8, height = 8/1.6
      )

    }
  }
}
