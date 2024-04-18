## Retrieve the data
load(file = "../data/processed/data_q2.RData")


## Strip out all unnecessary variables

full_data <- data
data <- full_data
rerun_models <- FALSE


## Start with an intercept only model
{
  ## Raw means
  {
    cellmeans_summ_raw <- data |>
     reframe(
       type = c("mean", "median"),
       Mean = c(mean(n.points / total.points), median(n.points / total.points)),
       SD = c(sd(n.points / total.points), MAD = mad(n.points / total.points)),
       N = c(n(), n())
     ) |>
     mutate(
       lower = Mean - 2 * (SD / sqrt(N)),
       upper = Mean + 2 * (SD / sqrt(N))
     )
    save(cellmeans_summ_raw, file = "../data/modelled/cellmeans_summ_raw_1.1.RData")
  }
  ## glmmTMB
  {
    ## fit model
    {
      if (rerun_models) {
        ## ---- q2_glmmTMB_fitted_1.1 
        mod_glmmTMB <- glmmTMB(
          cbind(n.points, total.points - n.points) ~ 1 +
            (1 | AIMS_REEF_NAME) +
            (1 | Site) +
            (1 | Transect),
          data = data,
          family = "betabinomial",
          REML = TRUE
        )
        save(mod_glmmTMB, file = "../data/modelled/mod_glmmTMB_1.1.RData")
        ## ----end
      }
    }
    ## Explore residuals (model validation)
    {
      ## ---- q2_glmmTMB_validation_1.1
      load(file = "../data/modelled/mod_glmmTMB_1.1.RData")
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
      ## ---- g2_glmmTMB_summary_1.1
      load(file = "../data/modelled/mod_glmmTMB_1.1.RData")
      ## summ <- list("Raw" = mod_glmmTMB, "Exponentiated" = mod_glmmTMB) |>
      ##   ## list("Exponentiated" = mod_glmmTMB) |>
      ##   modelsummary(
      ##     estimate = "{estimate} [{conf.low},{conf.high}]",
      ##     statistic = NULL,
      ##     shape = term + component ~ model,
      ##     group_map = "conditional",
      ##     exponentiate = c(FALSE, TRUE),
      ##     gof_map = "all",
      ##     fmt = function(x) format(x, digits = 3, nsmall = 2, scientific = FALSE, trim = TRUE)
      ##     ## fmt = function(x) format(x, digits = 3, nsmall = 2, scientific = TRUE, trim = TRUE)
      ##   )

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
      save(summ, file = "../data/summ_glmmTMB_1.1.RData")
      save(summ_tt, file = "../data/summ_tt_glmmTMB_1.1.RData")
        ## mutate(across(c(estimate, conf.low, conf.high), \(x) sprintf("%.3f", x))) |> 
        ## mutate(conf = paste0("[", conf.low, ",", conf.high, "]")) |>
        ## dplyr::select(-conf.low, -conf.high) |> 
        ## pivot_longer(cols = c("estimate", "conf"),
        ##   names_to = "names", values_to = "values") |> 
        ## pivot_wider(id_cols = everything(), names_from = "model", values_from = "values") |>
        ## dplyr::select(-effect, -names) |> 
        ## tinytable::tt() |>
        ## ## tinytable::style_tt(i = c(1, 3, 5, 7, 9), j = 1, rowspan = 2, alignv = "t")

      ## ----end
      ## ---- g2_glmmTMB_summary_1.1_print
      load(file = "../data/summ_tt_glmmTMB_1.1.RData")
      summ_tt
      ## ----end
    }
    ## Partial effects
    {
      ## ---- g2_glmmTMB_parial_effects_1.1
      load(file = "../data/modelled/mod_glmmTMB_1.1.RData")
      cellmeans_summ_glmmTMB <- mod_glmmTMB |>
        emmeans(~1, type = "response") |>
        as.data.frame()
      cellmeans_summ_glmmTMB |> ggplot(aes(y = prob, x = 1)) +
        geom_pointrange(aes(ymin = asymp.LCL, ymax = asymp.UCL))
      save(cellmeans_summ_glmmTMB, file = "../data/modelled/cellmeans_summ_glmmTMB_1.1.RData")
      ## ----end
    }
  }
  ## brms
  {
    ## fit model
    {
      if (rerun_models) {
        ## ---- q2_brm_fitted_1.1
        form <- bf(
          n.points | trials(total.points) ~ 1 +
            (1 | AIMS_REEF_NAME) +
            (1 | Site) +
            (1 | Transect),
          family = "beta_binomial"
        )
        data |>
          summarise(qlogis(mean(n.points / total.points)))
        priors <- prior(normal(-4.70, 1), class = "Intercept") +
          prior(student_t(3, 0, 1), class = "sd") +
          prior(gamma(0.01, 0.01), class = "phi")
        mod_brm <- brm(
          formula = form,
          data = data,
          prior = priors,
          iter = 5000, warmup = 1000,
          chains = 3, cores = 3,
          sample_prior = "yes",
          thin = 5,
          backend = "cmdstanr",
          control = list(adapt_delta = 0.99)
        )
        summary(mod_brm)
        save(mod_brm, file = "../data/modelled/mod_brm_1.1.RData")
        ## ----end
      }
    }
    ## MCMC diagnostics
    {
      ## ---- q2_brm_trace_1.1a
      load(file = "../data/modelled/mod_brm_1.1.RData")
      rstan::stan_trace(mod_brm$fit)
      rstan::stan_ac(mod_brm$fit)
      rstan::stan_ess(mod_brm$fit) + rstan::stan_rhat(mod_brm$fit)
      ## ----end
      
      ## ---- q2_brm_trace_1.1
      load(file = "../data/modelled/mod_brm_1.1.RData")
      ggsave(
        filename = "../docs/analysis_files/figure-html/stan_trace_brm_1.1.png",
        rstan::stan_trace(mod_brm$fit),
        width = 15, height = 8
      )
      ggsave(
        filename = "../docs/analysis_files/figure-html/stan_ac_brm_1.1.png",
        rstan::stan_ac(mod_brm$fit),
        width = 15, height = 8
      )
      ggsave(
        filename = "../docs/analysis_files/figure-html/stan_ess_brm_1.1.png",
        rstan::stan_ess(mod_brm$fit) + rstan::stan_rhat(mod_brm$fit),
        width = 15, height = 8
      )
      ## ----end
    }
    ## Model validation
    {
      ## ---- q2_brm_validation_1.1a
      load(file = "../data/modelled/mod_brm_1.1.RData")
      mod_brm |> pp_check(type = "dens_overlay", ndraws = 100) + scale_x_log10() +
        mod_brm |> pp_check(type = "loo_pit_overlay", ndraws = 100)
      ## ----end

      ## ---- q2_brm_validation_1.1
      load(file = "../data/modelled/mod_brm_1.1.RData")
      ggsave(
        filename = "../docs/analysis_files/figure-html/pcc_brm_1.1.png",
        mod_brm |> pp_check(type = "dens_overlay", ndraws = 100) + scale_x_log10() +
          mod_brm |> pp_check(type = "loo_pit_overlay", ndraws = 100),
        width = 10,
        height = 5)
      ## ----end
    }
    ## DHARMa residuals
    {
      ## ---- q2_brm_dharma_1.1a
      load(file = "../data/modelled/mod_brm_1.1.RData")
      resids <- make_brms_dharma_res(
        mod_brm,
        integerResponse = FALSE
      )
      wrap_elements(~ testUniformity(resids)) +
        wrap_elements(~ plotResiduals(resids, form = factor(rep(1, nrow(mod_brm$data))))) +
        wrap_elements(~ plotResiduals(resids)) +
        wrap_elements(~ testDispersion(resids))
      ## ----end

      ## ---- q2_brm_dharma_1.1
      load(file = "../data/modelled/mod_brm_1.1.RData")
      resids <- make_brms_dharma_res(
        mod_brm,
        integerResponse = FALSE
      )
      ggsave(
        filename = "../docs/analysis_files/figure-html/dharma_brm_1.1.png",
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
      ## ---- q2_brm_summary_1.1
      load(file = "../data/modelled/mod_brm_1.1.RData")

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
      save(summ, file = "../data/summ_brm_1.1.RData")
      save(summ_tt, file = "../data/summ_tt_brm_1.1.RData")
      ## mod_brm |>
      ##   tidybayes::tidy_draws() |>
      ##   tidybayes::summarise_draws(
      ##     median,
      ##     HDInterval::hdi,
      ##     rhat,
      ##     ess_bulk,
      ##     ess_tail
      ##   ) |>
      ##   filter(str_detect(variable, "^b_|^sd_|^sigma")) |>
      ##   knitr::kable()
      ## ----end
      ## ---- q2_brm_summary_1.1_print
      load(file = "../data/summ_tt_brm_1.1.RData")
      summ_tt
      ## ----end
    }
    ## Partial effects
    {
      ## ---- q2_brm_partial_1.1
      load(file = "../data/modelled/mod_brm_1.1.RData")
      cellmeans_summ_brm <- mod_brm |>
        emmeans(~1, type = "response") |>
        as.data.frame()
      save(cellmeans_summ_brm, file = "../data/modelled/cellmeans_summ_brm_1.1.RData")

      p <- 
      cellmeans_summ_brm |> ggplot(aes(y = prob, x = 1)) +
          geom_pointrange(aes(ymin = lower.HPD, ymax = upper.HPD))
      ggsave(
        filename = "../docs/analysis_files/figure-html/partial_brm_1.1.png",
        p,
        width = 6, height = 4
      )
      ## ----end

    }
  }
  ## INLA
  {
    ## Prepare data
    {
      ## ---- q2_inla_prepare_1.1
      data.pred <- data |>
        dplyr::select(
          n.points, total.points, AIMS_REEF_NAME,
          Site, Transect
        )
      newdata <- data.frame(AIMS_REEF_NAME = NA, Site = NA, Transect = NA)
      data_pred <- data |> bind_rows(newdata)
      i_newdata <- (nrow(data) + 1):nrow(data_pred) 
      ## ----end
    }
    ## Fit model
    {
      if (rerun_models) {
        ## ---- q2_inla_fitted_1.1
        mod <- inla(n.points ~ 1 +
                      f(model = "iid", AIMS_REEF_NAME) +
                      f(model = "iid", Site) +
                      f(model = "iid", Transect),
          data = data.pred,
          Ntrials = data.pred$total.points,
          family = "betabinomial", 
          control.predictor = list(link = 1, compute = TRUE),
          control.compute = list(config = TRUE,
            dic = TRUE, waic = TRUE, cpo = TRUE
          )
        )
        summary(mod)
        ## autoplot(mod)
        save(mod, file = "../data/modelled/mod_inla_1.1.RData")
        draws <- inla.posterior.sample(1000, result=mod, seed=123) %>% suppressWarnings()
        save(draws, file = "../data/modelled/draws_inla_1.1.RData")
        ## ----end
      }
    }
    ## Diagnostics
    {
      ## ---- q2_inla_diagnostics_1.1a
      load(file = "../data/modelled/mod_inla_1.1.RData")
      p <- 
      wrap_plots(
        pit_qq_plot(mod, i.mod = 1:nrow(data), logit_scale = TRUE),
        pit_plot(mod, i.mod = 1:nrow(data)),
        pit_resid_plot(mod, i.mod = 1:nrow(data)),
        pit_hist_plot(mod, i.mod = 1:nrow(data))
      ) +
        plot_layout(ncol = 2)

      ggsave(
        filename = "../docs/analysis_files/figure-html/diagnostics1a_inla_1.1.png",
        p,
        width = 10, height = 8
      )
      ## ----end

      ## ---- q2_inla_diagnostics_1.1b
      mod_cpo <- data.frame(cpo = mod$cpo$cpo,
        pit = mod$cpo$pit,
        failure = mod$cpo$failure) |>
        filter(failure == 0) |>
        mutate(row = 1:n())

      mod$cpo$failure

      g1 <- mod_cpo |>
        ggplot(aes(y = cpo, x = row)) +
        geom_point()
      g2 <- mod_cpo |>
        ggplot(aes(x = cpo)) +
        geom_histogram()
      g3 <- mod_cpo |>
        ggplot(aes(y = pit, x = row)) +
        geom_point()
      g4 <- mod_cpo |>
        ggplot(aes(x = pit)) +
        geom_histogram()
      p<-
      (g1 + g2)/(g3 + g4)

      ggsave(
        filename = "../docs/analysis_files/figure-html/diagnostics1b_inla_1.1.png",
        p,
        width = 10, height = 8
      )

      ## Ideally, we want to see that no CPO or PIT is very different in
      ## value from any others (often not informative in binomial models)
      ## and that the histograms are relatively flat (not really the case
      ## here).

      ## ggplot_inla_residuals(mod, observed = data$n.points)
      ## ggplot_inla_residuals2(mod, observed = data$n.points)
      ## ----end
    }
    ## DHARMa 
    {
      ## ---- q2_inla_dharma_1.1
      load(file = "../data/modelled/mod_inla_1.1.RData")
      preds <- posterior_predict.inla(mod, newdata = data)
      mod_resids <- createDHARMa(simulatedResponse = t(preds),
        observedResponse = data$n.points,
        fittedPredictedResponse = apply(preds, 2, mean),
        integerResponse = TRUE)
      p <- 
      wrap_elements(~testUniformity(mod_resids)) +
        wrap_elements(~plotResiduals(mod_resids)) +
        wrap_elements(~testDispersion(mod_resids)) 

      ggsave(
        filename = "../docs/analysis_files/figure-html/dharma_inla_1.1.png",
        p,
        width = 12,
        height = 8,
        dpi = 100
      )
      ## testDispersion(mod_resids)
      ## testZeroInflation(mod_resids)
      ## ----end
    }
    ## Summaries
    {
      ## ---- q2_inla_summary_1.1
      load(file = "../data/modelled/mod_inla_1.1.RData")
      load(file = "../data/modelled/draws_inla_1.1.RData")

      summ <- 
      bind_rows(
        model_summary(mod, draws) |> mutate(model = "Raw"),
        model_summary(mod, draws, exponentiate = TRUE) |> mutate(model = "Exponentiated")
      ) |>
        mutate(across(c(estimate, conf.low, conf.high), \(x) sprintf("%.3f", x))) |> 
        mutate(estimate = paste0(estimate, "<br>[", conf.low, ",", conf.high, "]")) |>
        mutate(estimate = str_replace(estimate, "<br>\\[NA,NA\\]","")) |> 
        dplyr::select(-conf.low, -conf.high) |> 
        pivot_longer(cols = c("estimate"),
          names_to = "names", values_to = "values") |> 
        pivot_wider(id_cols = everything(), names_from = "model", values_from = "values") |>
          dplyr::select(-effect, -names) |>
        mutate(term = ifelse(term == "dispersion", "phi", term))
      summ_tt <- summ |> 
        tinytable::tt() 
      save(summ, file = "../data/summ_inla_1.1.RData")
      save(summ_tt, file = "../data/summ_tt_inla_1.1.RData")
      ## ----end
      ## ---- q2_inla_summary_1.1_print
      load(file = "../data/summ_tt_inla_1.1.RData")
      summ_tt
      ## ----end
      ## ---- q2_inla_summaries_1.1
      load(file = "../data/modelled/mod_inla_1.1.RData")
      summary(mod)
      ## ----end
    }
    ## Partial effects
    {
      ## ---- q2_inla_partial_1.1
      load(file = "../data/modelled/mod_inla_1.1.RData")
      draws <- inla.posterior.sample(n=1000, result = mod)
      contents <- mod$misc$configs$contents
      i_newdata <- contents$start[contents$tag == "(Intercept)"]
      cellmeans <- newdata |> cbind(t(sapply(draws, function(x) x$latent[i_newdata])))
      cellmeans_summ_inla <- cellmeans |>
        pivot_longer(cols = matches("^[0-9]*$"), names_to = ".draws") |>
        posterior::as_draws() |>
        mutate(.draw = .draws) |>
        dplyr::select(-.draws, -AIMS_REEF_NAME, -Site, -Transect) |>
        mutate(value = plogis(value)) |> 
        summarise_draws(median, HDInterval::hdi)

      save(cellmeans_summ_inla, file = "../data/modelled/cellmeans_summ_inla_1.1.RData")
      p <- 
      cellmeans_summ_inla |> ggplot(aes(y = median, x = 1)) +
          geom_pointrange(aes(ymin = lower, ymax = upper))
      ggsave(
        filename = "../docs/analysis_files/figure-html/partial_inla_1.1.png",
        p,
        width = 6, height = 4
      )
      ## ----end
    }
  }
  ## Comparisons
  {
    ## ---- q2_comparison_table_1.1
    summ_tt <- 
      bind_rows(
        get(load(file = "../data/summ_glmmTMB_1.1.RData")) |> mutate(routine = "glmmTMB"),
        get(load(file = "../data/summ_brm_1.1.RData")) |> mutate(routine = "brm"),
        get(load(file = "../data/summ_inla_1.1.RData")) |> mutate(routine = "inla")
      ) |> 
      dplyr::select(-Raw) |> 
      pivot_wider(id_cols = everything(), names_from = "routine", values_from = "Exponentiated") |>
      tinytable::tt() 
    save(summ_tt, file = "../data/summ_tt_comparison_1.1.RData")
    ## ----end
    ## ---- q2_comparison_table_1.1_print
    load(file = "../data/summ_tt_comparison_1.1.RData")
    summ_tt
    ## ----end
    ## ---- q2_comparison_plot_1.1
    load(file = "../data/modelled/cellmeans_summ_raw_1.1.RData")
    load(file = "../data/modelled/cellmeans_summ_glmmTMB_1.1.RData")
    load(file = "../data/modelled/cellmeans_summ_brm_1.1.RData")
    load(file = "../data/modelled/cellmeans_summ_inla_1.1.RData")

    cellmeans_summ <- bind_rows(
      cellmeans_summ_raw |>
        mutate(method = "raw") |>
        dplyr::rename(median = Mean, lower = lower, upper = upper),
      cellmeans_summ_glmmTMB |>
        mutate(method = "glmmTMB", type = "median") |>
        dplyr::rename(median = prob, lower = asymp.LCL, upper = asymp.UCL),
      cellmeans_summ_brm |>
        mutate(method = "brm", type = "median") |>
        dplyr::rename(median = prob, lower = lower.HPD, upper = upper.HPD),
      cellmeans_summ_inla |>
        mutate(method = "inla", type = "median") |> 
        dplyr::rename(median = median, lower = lower, upper = upper)
    )
    cellmeans_summ

    p <-
    cellmeans_summ |> ggplot(aes(y = median, x = method, shape = type)) +
            geom_pointrange(aes(ymin = lower, ymax = upper))
    ggsave(
      filename = "../docs/analysis_files/figure-html/comparison_plot_1.1.png",
      p,
      width = 6, height = 4
    )

  }
}
