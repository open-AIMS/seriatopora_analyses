
library(tidyverse)
library(INLA)
library(glmmTMB)
library(brms)
library(posterior)
library(tidybayes)
library(DHARMa)
library(emmeans)
library(INLAutils)
library(patchwork)
library(standist)
library(inlatools)
library(bayesplot)
library(DT)
library(modelsummary)
config_modelsummary(factory_default = 'tinytable')
source("../R/helper_functions.R")
rerun_models <- TRUE

## Retrieve data
## ---- q1_retrieve_data_1
load(file = "../data/processed/data_q1.RData")
full_data <- data
## ----end

## Presence absence
{
  ## Prepare the data
  {
    ## ---- q1_prepare_data_for_model_1
    data <- data |>
      dplyr::select(
        P_CODE, AIMS_REEF_NAME, zone_depth, Site, Transect, n.points,
        total.points, cREPORT_YEAR, A_SECTOR, SHELF
      ) |>
      mutate(
        SecShelfYr = paste(A_SECTOR, SHELF, cREPORT_YEAR),
        SecShelf = paste(A_SECTOR, SHELF)
      )
    data <- data |>
      mutate(
        SecShelfYr = forcats::fct_relevel(SecShelfYr, "IN M 2016")
        ## SecShelfYr = forcats::fct_relevel(SecShelfYr, "SW M 2014")
      )
    data <- data |>
      mutate(PA = ifelse(n.points == 0, 0, 1))
    ## ----end
  }
  ## glmmTMB
  {
    ## Fit model
    {
      if (rerun_models) {
        ## ---- q1_glmmTMB_fit_1
        mod_glmmTMB <- glmmTMB(
          PA ~ SecShelfYr +
            ## (1 | SecShelfYr) +
            (1 | AIMS_REEF_NAME) +
            (1 | zone_depth) +
            (1 | Site) +
            (1 | Transect),
          data = data,
          family = "binomial",
          REML = TRUE
        )
        save(mod_glmmTMB, file = "../data/modelled/mod_q1_glmmTMB_1.RData")
        ## ----end
      }
    }
    ## DHARMa
    {
      ## ---- q1_glmmTMB_validation_1
      load(file = "../data/modelled/mod_q1_glmmTMB_1.RData")
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
      ## ---- q1_glmmTMB_summaries_print_1
      load(file = "../data/modelled/mod_q1_glmmTMB_1.RData")
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
      save(summ, file = "../data/summ_glmmTMB_q1_1.RData")
      save(summ_tt, file = "../data/summ_tt_glmmTMB_q1_1.RData")
      ## ----end
      ## ---- q1_glmmTMB_summaries_1
      load(file = "../data/summ_glmmTMB_q1_1.RData")
      summ_tt
      ## ----end
    }
    ## Partial plots
    {
      ## ---- q1_glmmTMB_partial_effects_1
      load(file = "../data/modelled/mod_q1_glmmTMB_1.RData")
      mod_em <-
        mod_glmmTMB |>
        emmeans(~SecShelfYr, type = "response") |>
        as.data.frame() |>
        separate(SecShelfYr, into = c("A_SECTOR", "SHELF", "cREPORT_YEAR"), remove = FALSE) |>
        mutate(A_SECTOR = factor(A_SECTOR, levels = c("CG", "PC", "CL", "CA", "IN", "TO", "WH", "PO", "SW", "CB")))
      p <- mod_em |>
        ggplot(aes(y = prob, x = as.numeric(as.character(cREPORT_YEAR)))) +
        geom_line() +
        geom_pointrange(aes(ymin = asymp.LCL, ymax = asymp.UCL)) +
        facet_grid(A_SECTOR ~ SHELF, scales = "free") +
        scale_y_continuous(limits = c(0, 1))
      ggsave(filename = "../docs/analysis_q1_files/figure-html/partial_glmmTMB_q1_1.png", p, width = 10, height = 10)
      ## ----end
    }
  }
  ## brms
  {
    ## Fit model
    ## Performed on the HPC (23_fit_models_HPC.R)
    {
      if (rerun_models) {
        ## ---- q1_brm_fit_1
        form <- bf(
          PA | trials(1) ~ SecShelfYr +
            (1 | AIMS_REEF_NAME) +
            (1 | zone_depth) +
            (1 | Site) +
            (1 | Transect),
          family = binomial()
        )
        data |>
          group_by(SecShelfYr) |>
          summarise(
            M = qlogis(mean(PA, na.rm = TRUE)),
            S = qlogis(sd(PA, na.rm = TRUE)),
          )
        priors <- prior(normal(0, 1), class = "Intercept") +
          prior(normal(0, 1), class = "b") +
          prior(student_t(3, 0, 1), class = "sd")
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
        save(mod_brm, file = "../data/modelled/mod_q1_brm_1.RData")
        ## ----end
      }
    }
    ## MCMC diagnostics
    {
      ## ---- q1_brm_trace_1
      ## load(file = "../data/modelled/mod_q1_brm.RData")
      ## system("scp mlogan@hpc-l001.aims.gov.au:~/Work/AIMS/LTMP/Seriatopora/data/modelled/mod_q1_brm_1.RData ../data/modelled/mod_q1_brm_1.RData")
      load(file = "../data/modelled/mod_q1_brm_1.RData")
      rstan::stan_trace(mod_brm$fit)
      rstan::stan_ac(mod_brm$fit)
      rstan::stan_ess(mod_brm$fit) + rstan::stan_rhat(mod_brm$fit)
      ## ----end

      ## ---- q1_brm_trace_1
      load(file = "../data/modelled/mod_q1_brm_1.RData")
      ggsave(
        filename = "../docs/analysis_files/figure-html/stan_trace_q1_brm_1.png",
        rstan::stan_trace(mod_brm$fit),
        width = 15, height = 8
      )
      ggsave(
        filename = "../docs/analysis_files/figure-html/stan_ac_q1_brm_1.png",
        rstan::stan_ac(mod_brm$fit),
        width = 15, height = 8
      )
      ggsave(
        filename = "../docs/analysis_files/figure-html/stan_ess_q1_brm_1.png",
        rstan::stan_ess(mod_brm$fit) + rstan::stan_rhat(mod_brm$fit),
        width = 15, height = 8
      )
      ## ----end
    }
    ## Model validation
    {
      ## ---- q1_brm_validation_1_a
      load(file = "../data/modelled/mod_q1_brm_1.RData")
      mod_brm |>
        pp_check(type = "dens_overlay", ndraws = 100) +
        mod_brm |>
        pp_check(type = "loo_pit_overlay", ndraws = 100)
      ## ----end

      ## ---- q1_brm_validation_1
      load(file = "../data/modelled/mod_q1_brm_1.RData")
      ggsave(
        filename = "../docs/analysis_files/figure-html/pcc_q1_brm_1.png",
        mod_brm |> pp_check(type = "dens_overlay", ndraws = 100) +
          mod_brm |> pp_check(type = "loo_pit_overlay", ndraws = 100),
        width = 10,
        height = 5
      )
      ## ----end
    }
      ## DHARMa residuals
    {
      ## ---- q1_brm_dharma_1_a
      load(file = "../data/modelled/mod_q1_brm_1.RData")
      resids <- make_brms_dharma_res(
        mod_brm,
        integerResponse = FALSE
      )
      wrap_elements(~ testUniformity(resids)) +
        wrap_elements(~ plotResiduals(resids, form = factor(rep(1, nrow(mod_brm$data))))) +
        wrap_elements(~ plotResiduals(resids)) +
        wrap_elements(~ testDispersion(resids))
      ## ----end

      ## ---- q1_brm_dharma_1
      load(file = "../data/modelled/mod_q1_brm_1.RData")
      resids <- make_brms_dharma_res(
        mod_brm,
        integerResponse = FALSE
      )
      ggsave(
        filename = "../docs/analysis_files/figure-html/dharma_q1_brm_1.png",
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
      ## ---- q1_brm_summary_1
      load(file = "../data/modelled/mod_q1_brm_1.RData")

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
      save(summ, file = "../data/summ_q1_brm_1.RData")
      save(summ_tt, file = "../data/summ_tt_q1_brm_1.RData")
      ## ----end
      ## ---- q1_brm_summary_print_1
      load(file = "../data/summ_tt_q1_brm_1.RData")
      summ_tt
      ## ----end
    }
    ## Partial plots -  Type 1. sector/shelf re_formula = NULL
    {
      ## ---- q1_brm_partial_effects_1
      load(file = "../data/modelled/mod_q1_brm_1.RData")
      newdata <- brm_generate_newdata_no_dist(data)
      pred <- t(brms::posterior_epred(mod_brm,
        newdata = newdata,
        re_formula = NULL,
        se.fit = TRUE,
        allow_new_levels = TRUE
      ))
      cellmeans_brm <- brm_generate_cellmeans_no_dist(newdata, pred)

      save(cellmeans_brm, file = "../data/modelled/cellmeans_q1_brm_type_1_1.RData")

      cellmeans_summ_brm <- cellmeans_brm |>
        dplyr::select(-.draws) |>
        dplyr::select(A_SECTOR, SHELF, REPORT_YEAR, Pred) |>
        group_by(A_SECTOR, SHELF, REPORT_YEAR) |>
        summarise_draws(median, HDInterval::hdi) |>
        dplyr::select(A_SECTOR, SHELF, REPORT_YEAR, median, lower, upper)

      save(cellmeans_summ_brm, file = "../data/modelled/cellmeans_summ_q1_brm_type_1_1.RData")

      p <- brm_partial_plot_no_dist(cellmeans_summ_brm) +
        scale_y_continuous("Probability of Seriatopora occurance",
          trans = scales::sqrt_trans(),
          breaks = c(0, 0.01, 0.1, 0.25, 0.5, 0.75, 1),
          limits = c(0, 1)
          ## trans = scales::log10_trans()
          ## trans = scales::boxcox_trans(p = 0.25)
          ## trans = scales::pseudo_log_trans(sigma = 1, base = 10)
          ## trans = scales::probit_trans()
          ## trans = scales::probability_trans()
        )
      ggsave(
        filename = "../docs/analysis_files/figure-html/partial_q1_brm_type_1_1.png",
        p,
        width = 15, height = 12
      )
      ## ----end
    }
    ## Partial plots -  Type 2. sector/shelf re_formula = NA
    {
      ## ---- q1_brm_partial_effects_type_2
      load(file = "../data/modelled/mod_q1_brm_1.RData")
      newdata <- brm_generate_newdata_no_dist(data)
      pred <- t(brms::posterior_epred(mod_brm,
        newdata = newdata,
        re_formula = ~SecShelfYr,
        se.fit = TRUE,
        allow_new_levels = FALSE
      ))
      cellmeans_brm <- brm_generate_cellmeans_no_dist(newdata, pred)

      save(cellmeans_brm, file = "../data/modelled/cellmeans_q1_brm_type_2_1.RData")

      cellmeans_summ_brm <- cellmeans_brm |>
        dplyr::select(-.draws) |>
        dplyr::select(A_SECTOR, SHELF, REPORT_YEAR, Pred) |>
        group_by(A_SECTOR, SHELF, REPORT_YEAR) |>
        summarise_draws(median, HDInterval::hdi) |>
        dplyr::select(A_SECTOR, SHELF, REPORT_YEAR, median, lower, upper)

      save(cellmeans_summ_brm, file = "../data/modelled/cellmeans_summ_q1_brm_type_2_1.RData")

      p <- brm_partial_plot_no_dist(cellmeans_summ_brm)
      ggsave(
        filename = "../docs/analysis_files/figure-html/partial_q1_brm_type_2_1.png",
        p,
        width = 15, height = 12
      )
      ## ----end
    }

  }
}

if (1 == 2) {
## Start with an intercept only model
{
  ## INLA
  {
    ## Data preparation
    {
      ## Focus on only the necessary variables
      data <- full_data |>
        dplyr::select(
          AIMS_REEF_NAME, zone_depth, Site, Transect, n.points,
          total.points
        )
      data.pred <- data
    }
    ## Fit model
    {
      mod <- inla(n.points ~ 1 +
                    f(model = "iid", AIMS_REEF_NAME) +
                    f(model = "iid", Site) +
                    f(model = "iid", Transect),
        data = data.pred,
        Ntrials = data.pred$total.points,
        family = "binomial", 
        control.predictor = list(link = 1, compute = TRUE),
        control.compute = list(config = TRUE,
          dic = TRUE, waic = TRUE, cpo = TRUE
        )
      )
    }
    ## Diagnostics
    {
      wrap_plots(
        pit_qq_plot(mod, i.mod = 1:nrow(data), logit_scale = TRUE),
        pit_plot(mod, i.mod = 1:nrow(data)),
        pit_resid_plot(mod, i.mod = 1:nrow(data)),
        pit_hist_plot(mod, i.mod = 1:nrow(data))
      ) +
      plot_layout(ncol = 3)

      mod.cpo <- data.frame(cpo = mod$cpo$cpo,
        pit = mod$cpo$pit,
        failure = mod$cpo$failure) |>
        filter(failure == 0) |>
        mutate(row = 1:n())

      mod$cpo$failure

      g1 <- mod.cpo |>
        ggplot(aes(y = cpo, x = row)) +
        geom_point()
      g2 <- mod.cpo |>
        ggplot(aes(x = cpo)) +
        geom_histogram()
      g3 <- mod.cpo |>
        ggplot(aes(y = pit, x = row)) +
        geom_point()
      g4 <- mod.cpo |>
        ggplot(aes(x = pit)) +
        geom_histogram()

      (g1 + g2)/(g3 + g4)

   ## Ideally, we want to see that no CPO or PIT is very different in
   ## value from any others (often not informative in binomial models)
   ## and that the histograms are relatively flat (not really the case
   ## here).

   ## ggplot_inla_residuals(mod, observed = data$n.points)
   ## ggplot_inla_residuals2(mod, observed = data$n.points)
    }
    ## DHARMa version 1
    {
      preds <- posterior_predict.inla(mod, newdata = data)
      mod_resids <- createDHARMa(simulatedResponse = t(preds),
        observedResponse = data$n.points,
        fittedPredictedResponse = apply(preds, 2, mean),
        integerResponse = TRUE)
      wrap_elements(~testUniformity(mod_resids)) +
        wrap_elements(~plotResiduals(mod_resids)) +
        wrap_elements(~testDispersion(mod_resids)) 

      mod_resids |> testDispersion()
      mod_resids |> testZeroInflation()
      }
  }
  ## glmmTMB
  {
    library(glmmTMB)
    mod_glmmTMB <- glmmTMB(cbind(n.points, total.points-n.points) ~ 1 +
                             (1|AIMS_REEF_NAME) +
                             (1|Site) +
                             (1|Transect),
      data = data,
      family = "binomial", 
      REML = TRUE
    )
    resids <- simulateResiduals(mod_glmmTMB, plot = TRUE)
    resids |> testDispersion()
    resids |> testZeroInflation()
  }
}

## add REPORT_YEAR
{
  ## INLA
  {
    ## Prepare data
    {
      data <- full_data

      data <- data |>
              dplyr::select(
                AIMS_REEF_NAME, Site, Transect, n.points,
                total.points, cREPORT_YEAR
              )
      ## Cellmeans
      newdata <- data.frame(cREPORT_YEAR = sort(unique(data$cREPORT_YEAR)))

      data_pred <- data |> bind_rows(newdata)
      i_newdata <- (nrow(data) + 1):nrow(data_pred)
    }
    ## Fit model
    {
      mod <- inla(n.points ~ cREPORT_YEAR +
                    f(model = "iid", AIMS_REEF_NAME) +
                    f(model = "iid", Site) +
                    f(model = "iid", Transect),
        data = data_pred,
        Ntrials = data_pred$total.points,
        family = "zeroinflatedbinomial2", #"binomial", 
        control.predictor = list(link = 1, compute = TRUE),
        control.compute = list(config = TRUE,
          dic = TRUE, waic = TRUE, cpo = TRUE
        )
      )
      ## summary(mod)
      ## autoplot(mod)
    }

    ## Diagnostics
    {
      ggplot_inla_residuals(mod,
        observed = data$n.points / data$total.points, CI = TRUE
      )
      ggplot_inla_residuals2(mod,
        observed = data$n.points / data$total.points, CI = TRUE
      )
      wrap_plots(
        pit_qq_plot(mod, i.mod = 1:nrow(data), logit_scale = TRUE),
        pit_plot(mod, i.mod = 1:nrow(data)),
        pit_resid_plot(mod, i.mod = 1:nrow(data)),
        pit_hist_plot(mod, i.mod = 1:nrow(data))
      ) +
      plot_layout(ncol = 3)

      mod.cpo <- data.frame(cpo = mod$cpo$cpo,
        pit = mod$cpo$pit,
        failure = mod$cpo$failure) |>
        filter(failure == 0) |>
        mutate(row = 1:n())

      mod$cpo$failure

      g1 <- mod.cpo |>
        ggplot(aes(y = cpo, x = row)) +
        geom_point()
      g2 <- mod.cpo |>
        ggplot(aes(x = cpo)) +
        geom_histogram()
      g3 <- mod.cpo |>
        ggplot(aes(y = pit, x = row)) +
        geom_point()
      g4 <- mod.cpo |>
        ggplot(aes(x = pit)) +
        geom_histogram()

      (g1 + g2)/(g3 + g4)

      ## Ideally, we want to see that no CPO or PIT is very different in
      ## value from any others (often not informative in binomial models)
      ## and that the histograms are relatively flat (not really the case
      ## here).

      ## ggplot_inla_residuals(mod, observed = data$n.points)
      ## ggplot_inla_residuals2(mod, observed = data$n.points)
    }
    ## DHARMa version 1
    {
      preds <- posterior_predict.inla(mod, newdata = data)
      mod_resids <- createDHARMa(simulatedResponse = t(preds),
        observedResponse = data$n.points,
        fittedPredictedResponse = apply(preds, 2, mean),
        integerResponse = TRUE)
      mod_resids |> plot()
      mod_resids |> testDispersion()
      mod_resids |> testZeroInflation()
    }
    ## Partial plots version 1
    {
      newdata_pred <- newdata |>
        bind_cols(mod$summary.fitted.values[i_newdata, ])
      newdata_pred 

      newdata_pred |>
        ggplot(aes(y = `0.5quant`,
          x = as.numeric(as.character(cREPORT_YEAR)))) +
        geom_line() +
                geom_pointrange(aes(
                        ymin = `0.025quant`,
                        ymax = `0.975quant`
                )) +
        scale_x_continuous("")
    }
    ## Partial plots version 2
    {
      newdata_fitted <- posterior_fitted.inla(mod, newdata = newdata)
      newdata_fitted <- newdata_fitted |>
              group_by(cREPORT_YEAR) |>
              summarise_draws(
                      median,
                      HDInterval::hdi
              )
      newdata_fitted 

      newdata_fitted |>
        ggplot(aes(y = median,
          x = as.numeric(as.character(cREPORT_YEAR)))) +
        geom_line() +
                geom_pointrange(aes(
                        ymin = lower,
                        ymax = upper
                )) +
        scale_x_continuous("")
    }
  }
  ## glmmTMB
  {
    ## Fit model
    {
      library(glmmTMB)
      mod_glmmTMB <- glmmTMB(cbind(n.points, total.points-n.points) ~ cREPORT_YEAR +
                               (1|AIMS_REEF_NAME) +
                               (1|Site) +
                               (1|Transect),
        ziformula = ~1,
        data = data,
        family = "binomial", 
        REML = TRUE
      )

      summary(mod_glmmTMB)
    }
    ## DHARMa
    {
      resids <- simulateResiduals(mod_glmmTMB, plot = TRUE)
      wrap_elements(~testUniformity(resids)) +
        wrap_elements(~plotResiduals(resids)) +
        wrap_elements(~testDispersion(resids)) 
      testDispersion(resids)
      testZeroInflation(resids)
    }
    ## Partial plots
    {
      mod_glmmTMB |>
        emmeans(~cREPORT_YEAR, type = "response") |>
        as.data.frame() |>
        ggplot(aes(y = prob, x = as.numeric(as.character(cREPORT_YEAR)))) +
        geom_line() +
        geom_pointrange(aes(ymin = asymp.LCL, ymax = asymp.UCL))
    }
  }
}

}

data <- full_data
## add REPORT_YEAR and Sector/Shelf
## Abundance
{
  ## Prepare the data
  {
    ## ---- q1_prepare_data_for_model
    data <- data |>
      dplyr::select(P_CODE, AIMS_REEF_NAME, zone_depth, Site, Transect, n.points,
        total.points, cREPORT_YEAR, A_SECTOR, SHELF) |>
      mutate(
        SecShelfYr = paste(A_SECTOR, SHELF, cREPORT_YEAR),
        SecShelf = paste(A_SECTOR, SHELF)
      )
    ## what combinations are missing
    data |>
      group_by(SecShelfYr) |>
      summarise(Points = sum(n.points)) |>
      arrange(desc(Points))
    ## 492 missing - lets exclude those
    data_with_excludes <- data
    save(data_with_excludes, file = "../data/processed/data_q1_with_excludes.RData")
    excludes <- data |>
      group_by(SecShelfYr) |>
      summarise(Points = sum(n.points)) |>
      filter(Points == 0) |>
      mutate(EXCLUDE = TRUE)
    data <- data |>
      left_join(excludes) |>
      filter(is.na(EXCLUDE)) |>
      droplevels() |>
      mutate(
        SecShelfYr = forcats::fct_relevel(SecShelfYr, "IN M 2016")
        ## SecShelfYr = forcats::fct_relevel(SecShelfYr, "SW M 2014")
      ) |>
      dplyr::select(-Points, -EXCLUDE)
    ## ----end
  }
  ## Intercept only model
  {
    
    ## glmmTMB
    {
      ## Fit model.1
      {
        if (rerun_models) {
          ## ---- q1_glmmTMB_fit.1
          mod_glmmTMB <- glmmTMB(
            cbind(n.points, total.points - n.points) ~ 0 +
              (1 | SecShelfYr) +
              (1 | AIMS_REEF_NAME) +
              (1 | zone_depth) +
              (1 | Site) +
              (1 | Transect),
            ziformula = ~ (1 |SecShelfYr) + (1 | zone_depth) + (1| Site),
            data = data,
            family = "binomial",
            REML = TRUE
          )

          save(mod_glmmTMB, file = "../data/modelled/mod_q1_glmmTMB.1.RData")
          ## ----end
        }
      }
      ## DHARMa
      {
        ## ---- q1_glmmTMB_validation
        load(file = "../data/modelled/mod_q1_glmmTMB.1.RData")
        resids <- simulateResiduals(mod_glmmTMB, plot = FALSE)
        wrap_elements(~testUniformity(resids)) +
          wrap_elements(~plotResiduals(resids)) +
          wrap_elements(~testDispersion(resids)) 
        ## ----end
        testDispersion(resids)
        testZeroInflation(resids)
      }
      ## Summaries
      {
        ## ---- q1_glmmTMB_summaries_print.1
        load(file = "../data/modelled/mod_q1_glmmTMB.1.RData")
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
        save(summ, file = "../data/summ_glmmTMB_q1.RData")
        save(summ_tt, file = "../data/summ_tt_glmmTMB_q1.RData")
        ## ----end
        ## ---- q1_glmmTMB_summaries
        load(file = "../data/summ_glmmTMB_q1.RData")
        summ_tt
        ## ----end
      }
      ## Partial plots
      {
        ## ---- q1_glmmTMB_partial_effects
        load(file = "../data/modelled/mod_q1_glmmTMB.1.RData")

        newdata <- brm_generate_newdata_no_dist(data) |>
          mutate(unit = NA)
        excludes <- excludes |>
          mutate(prob = 0, asymp.LCL = 0, asymp.UCL = 0) |>
          dplyr::select(-Points, -EXCLUDE) |>
          separate(SecShelfYr, into = c("A_SECTOR", "SHELF", "cREPORT_YEAR"), remove = FALSE) 
        fit <- mod_glmmTMB |>
          predict(
            newdata = newdata, re.form = NULL, allow.new.levels = TRUE,
            type = "link", se.fit = TRUE
          )
        mod_em <- newdata |> mutate(
          asymp.LCL = plogis(fit$fit - 2 * fit$se.fit),
          asymp.UCL = plogis(fit$fit + 2 * fit$se.fit),
          prob =  plogis(fit$fit)
        ) |> 
        ## mod_em <- 
        ##   mod_glmmTMB |>
        ##   emmeans(~SecShelfYr, type = "response") |>
        ##   as.data.frame() |>
          separate(SecShelfYr, into = c("A_SECTOR", "SHELF", "cREPORT_YEAR"), remove = FALSE) |>
          full_join(excludes) |>
          mutate(A_SECTOR = factor(A_SECTOR, levels = c("CG", "PC", "CL", "CA", "IN", "TO", "WH", "PO", "SW", "CB")))
        p <- mod_em |> 
          ggplot(aes(y = prob, x = as.numeric(as.character(cREPORT_YEAR)))) +
          geom_line() +
          geom_pointrange(aes(ymin = asymp.LCL, ymax = asymp.UCL)) +
          facet_grid(A_SECTOR~SHELF, scales = "free")
        ggsave(filename = "../docs/analysis_q1_files/figure-html/partial_glmmTMB_q1.1.png", width = 10, height = 10)
        ## ----end
      }
    }
  }
  
  ## Full model
  {
    ## glmmTMB
    {
      ## Fit model
      {
        if (rerun_models) {
          ## ---- q1_glmmTMB_fit
          mod_glmmTMB <- glmmTMB(
            cbind(n.points, total.points - n.points) ~ SecShelfYr +
              (1 | SecShelfYr) +
              (1 | AIMS_REEF_NAME) +
              (1 | zone_depth) +
              (1 | Site) +
              (1 | Transect),
            ziformula = ~ (1 |SecShelfYr) + (1 | zone_depth) + (1| Site),
            data = data,
            family = "binomial",
            REML = TRUE
          )

          save(mod_glmmTMB, file = "../data/modelled/mod_q1_glmmTMB.RData")
          ## ----end
        }
      }
      ## DHARMa
      {
        ## ---- q1_glmmTMB_validation
        load(file = "../data/modelled/mod_q1_glmmTMB.RData")
        resids <- simulateResiduals(mod_glmmTMB, plot = FALSE)
        wrap_elements(~testUniformity(resids)) +
          wrap_elements(~plotResiduals(resids)) +
          wrap_elements(~testDispersion(resids)) 
        ## ----end
        testDispersion(resids)
        testZeroInflation(resids)
      }
      ## Summaries
      {
        ## ---- q1_glmmTMB_summaries_print
        load(file = "../data/modelled/mod_q1_glmmTMB.RData")
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
        save(summ, file = "../data/summ_glmmTMB_q1.RData")
        save(summ_tt, file = "../data/summ_tt_glmmTMB_q1.RData")
        ## ----end
        ## ---- q1_glmmTMB_summaries
        load(file = "../data/summ_glmmTMB_q1.RData")
        summ_tt
        ## ----end
      }
      ## Partial plots
      {
        ## ---- q1_glmmTMB_partial_effects
        load(file = "../data/modelled/mod_q1_glmmTMB.RData")
        excludes <- excludes |>
          mutate(prob = 0, asymp.LCL = 0, asymp.UCL = 0) |>
          dplyr::select(-Points, -EXCLUDE) |>
          separate(SecShelfYr, into = c("A_SECTOR", "SHELF", "cREPORT_YEAR"), remove = FALSE) 
        mod_em <- 
          mod_glmmTMB |>
          emmeans(~SecShelfYr, type = "response") |>
          as.data.frame() |>
          separate(SecShelfYr, into = c("A_SECTOR", "SHELF", "cREPORT_YEAR"), remove = FALSE) |>
          full_join(excludes) |>
          mutate(A_SECTOR = factor(A_SECTOR, levels = c("CG", "PC", "CL", "CA", "IN", "TO", "WH", "PO", "SW", "CB")))
        p <- mod_em |> 
          ggplot(aes(y = prob, x = as.numeric(as.character(cREPORT_YEAR)))) +
          geom_line() +
          geom_pointrange(aes(ymin = asymp.LCL, ymax = asymp.UCL)) +
          facet_grid(A_SECTOR~SHELF, scales = "free")
        ggsave(filename = "../docs/analysis_q1_files/figure-html/partial_glmmTMB_q1.png", width = 10, height = 10)
        ## ----end
      }
    }
    if (1 == 2) {
      ## INLA
      {
        ## Prepare data
        {
          data <- full_data

          data <- data |>
            dplyr::select(
              AIMS_REEF_NAME, Site, Transect, n.points,
              total.points, cREPORT_YEAR, A_SECTOR, SHELF
            ) |>
            mutate(
              SecShelfYr = paste(A_SECTOR, SHELF, cREPORT_YEAR),
              SecShelf = paste(A_SECTOR, SHELF)
            )
          ## what combinations are missing
          data |>
            group_by(SecShelfYr) |>
            summarise(Points = sum(n.points)) |>
            arrange(desc(Points))
          ## 492 missing - lets exclude those

          excludes <- data |>
            group_by(SecShelfYr) |>
            summarise(Points = sum(n.points)) |>
            filter(Points == 0) |>
            mutate(EXCLUDE = TRUE)
          data <- data |>
            left_join(excludes) |>
            filter(is.na(EXCLUDE)) |>
            droplevels() |>
            mutate(
              SecShelfYr = forcats::fct_relevel(SecShelfYr, "IN M 2016")
              ## SecShelfYr = forcats::fct_relevel(SecShelfYr, "SW M 2014")
            ) |>
            dplyr::select(-Points, -EXCLUDE)

          ## Cellmeans
          newdata <- crossing(SecShelfYr = data$SecShelfYr) |>
            separate(SecShelfYr,
              into = c("A_SECTOR", "SHELF", "cREPORT_YEAR"),
              remove = FALSE
            )

          data_pred <- data |> bind_rows(newdata)
          i_newdata <- (nrow(data) + 1):nrow(data_pred)
        }
        ## Fit model
        {
          mod <- inla(
            n.points ~ SecShelfYr +
              f(model = "iid", AIMS_REEF_NAME) +
              f(model = "iid", Site) +
              f(model = "iid", Transect),
            data = data_pred,
            Ntrials = data_pred$total.points,
            family = "zeroinflatedbinomial1", # "zeroinflatedbetabinomial2", #"binomial",
            control.predictor = list(link = 1, compute = TRUE),
            control.compute = list(
              config = TRUE,
              dic = TRUE, waic = TRUE, cpo = TRUE
            ) # ,
            ## control.fixed = control.fixed
          )

          ## summary(mod)
          ## autoplot(mod)
        }

        ## Diagnostics
        {
          ggplot_inla_residuals(mod,
            observed = data$n.points / data$total.points, CI = TRUE
          )
          ggplot_inla_residuals2(mod,
            observed = data$n.points / data$total.points, CI = TRUE
          )
          wrap_plots(
            pit_qq_plot(mod, i.mod = 1:nrow(data), logit_scale = TRUE),
            pit_plot(mod, i.mod = 1:nrow(data)),
            pit_resid_plot(mod, i.mod = 1:nrow(data)),
            pit_hist_plot(mod, i.mod = 1:nrow(data))
          ) +
            plot_layout(ncol = 3)

          mod.cpo <- data.frame(
            cpo = mod$cpo$cpo,
            pit = mod$cpo$pit,
            failure = mod$cpo$failure
          ) |>
            filter(failure == 0) |>
            mutate(row = 1:n())

          mod$cpo$failure
          sum(mod$cpo$cop, na.omit = TRUE)
          -mean(log(na.omit(mod$cpo$cpo)))

          g1 <- mod.cpo |>
            ggplot(aes(y = cpo, x = row)) +
            geom_point()
          g2 <- mod.cpo |>
            ggplot(aes(x = cpo)) +
            geom_histogram()
          g3 <- mod.cpo |>
            ggplot(aes(y = pit, x = row)) +
            geom_point()
          g4 <- mod.cpo |>
            ggplot(aes(x = pit)) +
            geom_histogram()

          (g1 + g2) / (g3 + g4)

          ## Ideally, we want to see that no CPO or PIT is very different in
          ## value from any others (often not informative in binomial models)
          ## and that the histograms are relatively flat (not really the case
          ## here).

          ## ggplot_inla_residuals(mod, observed = data$n.points)
          ## ggplot_inla_residuals2(mod, observed = data$n.points)
        }
        ## DHARMa version 1
        {
          preds <- posterior_predict.inla(mod, newdata = data)
          mod_resids <- createDHARMa(
            simulatedResponse = t(preds),
            observedResponse = data$n.points,
            fittedPredictedResponse = apply(preds, 2, mean),
            integerResponse = TRUE
          )
          mod_resids |> plot()
          mod_resids |> testDispersion()
          mod_resids |> testZeroInflation()

          x <- list(
            y = na.omit(data$n.points / data$total.points),
            yrep = (preds)
          )
          class(x) <- "foo"
          pp_check.foo(x, type = "overlaid")
        }
        ## Partial plots
        {
          newdata_pred <- newdata |>
            bind_cols(mod$summary.fitted.values[i_newdata, ])
          newdata_pred

          ## Put the zeros back
          excludes <- excludes |>
            mutate(`0.5quant` = 0, `0.025quant` = 0, `0.975quant` = 0) |>
            dplyr::select(-Points, -EXCLUDE) |>
            separate(SecShelfYr,
              into = c("A_SECTOR", "SHELF", "cREPORT_YEAR"),
              remove = FALSE
            )

          newdata_pred <- newdata_pred |>
            full_join(excludes) |>
            mutate(A_SECTOR = factor(A_SECTOR,
              levels = c("CG", "PC", "CL", "CA", "IN", "TO", "WH", "PO", "SW", "CB")
            ))

          newdata_pred |>
            ggplot(aes(
              y = `0.5quant`,
              x = as.numeric(as.character(cREPORT_YEAR))
            )) +
            geom_line() +
            geom_pointrange(aes(ymin = `0.025quant`, ymax = `0.975quant`)) +
            facet_grid(A_SECTOR ~ SHELF, scales = "free")
        }
        ## Partial plots 2
        {
          newdata_fitted <- mod |>
            posterior_fitted.inla(newdata)

          newdata_fitted <- newdata_fitted |>
            group_by(A_SECTOR, SHELF, cREPORT_YEAR) |>
            dplyr::select(-SecShelfYr) |>
            posterior::summarise_draws(
              median,
              HDInterval::hdi
            )

          excludes <- excludes |>
            mutate(median = 0, lower = 0, upper = 0) |>
            dplyr::select(-any_of(c("Points", "EXCLUDE"))) |>
            separate(SecShelfYr,
              into = c("A_SECTOR", "SHELF", "cREPORT_YEAR"),
              remove = TRUE
            )

          newdata_fitted <- newdata_fitted |>
            full_join(excludes) |>
            mutate(A_SECTOR = factor(A_SECTOR,
              levels = c("CG", "PC", "CL", "CA", "IN", "TO", "WH", "PO", "SW", "CB")
            ))

          newdata_fitted |>
            ggplot(aes(y = median, x = as.numeric(as.character(cREPORT_YEAR)))) +
            geom_line() +
            geom_pointrange(aes(ymin = lower, ymax = upper)) +
            facet_grid(A_SECTOR ~ SHELF, scales = "free")
        }
      }
    }
    ## brms
    {
      ## Fit model
      ## Performed on the HPC (23_fit_models_HPC.R)
      {
        if (rerun_models) {
          ## ---- q1_brm_fit
          form <- bf(
            n.points | trials(total.points) ~ SecShelfYr +
              (1 | AIMS_REEF_NAME) +
              (1 | zone_depth) +
              (1 | Site) +
              (1 | Transect),
            ## family = "beta_binomial"
            zi ~ SecShelfYr + (1| zone_depth) + (1 | Site),
            family = zero_inflated_binomial()
          )
          data |>
            group_by(SecShelfYr) |>
            summarise(
              M = median(qlogis(n.points / total.points), na.rm = TRUE),
              S = mad(qlogis(n.points / total.points), na.rm = TRUE),
              Mean = qlogis(mean(n.points / total.points)),
              SD = qlogis(sd(n.points / total.points))
            )
          priors <- prior(normal(-4, 1), class = "Intercept") +
            prior(normal(0, 2), class = "b") +
            prior(student_t(3, 0, 1), class = "sd") +
            prior(logistic(0, 1), class = "Intercept", dpar = "zi") +
            prior(normal(0,1), class = 'b', dpar = 'zi') +
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
          ## save(mod_brm, file = "../data/modelled/mod_q1_brm.RData")
          save(mod_brm, file = "../data/modelled/mod_q1_brm.RData")
          ## ----end
        }
      }
      ## MCMC diagnostics
      {
        ## ---- q1_brm_trace
        ## load(file = "../data/modelled/mod_q1_brm.RData")
        ## system("scp mlogan@hpc-l001.aims.gov.au:~/Work/AIMS/LTMP/Seriatopora/data/modelled/mod_q1_brm_2.RData ../data/modelled/mod_q1_brm.RData")
        load(file = "../data/modelled/mod_q1_brm.RData")
        rstan::stan_trace(mod_brm$fit)
        rstan::stan_ac(mod_brm$fit)
        rstan::stan_ess(mod_brm$fit) + rstan::stan_rhat(mod_brm$fit)
        ## ----end
        
        ## ---- q1_brm_trace
        load(file = "../data/modelled/mod_q1_brm.RData")
        ggsave(
          filename = "../docs/analysis_files/figure-html/stan_trace_q1_brm.png",
          rstan::stan_trace(mod_brm$fit),
          width = 15, height = 8
        )
        ggsave(
          filename = "../docs/analysis_files/figure-html/stan_ac_q1_brm.png",
          rstan::stan_ac(mod_brm$fit),
          width = 15, height = 8
        )
        ggsave(
          filename = "../docs/analysis_files/figure-html/stan_ess_q1_brm.png",
          rstan::stan_ess(mod_brm$fit) + rstan::stan_rhat(mod_brm$fit),
          width = 15, height = 8
        )
        ## ----end
      }
      ## Model validation
      {
        ## ---- q1_brm_validation_a
        load(file = "../data/modelled/mod_q1_brm.RData")
        mod_brm |> pp_check(type = "dens_overlay", ndraws = 100) + scale_x_log10() +
          mod_brm |> pp_check(type = "loo_pit_overlay", ndraws = 100)
        ## ----end

        ## ---- q1_brm_validation
        load(file = "../data/modelled/mod_q1_brm.RData")
        ggsave(
          filename = "../docs/analysis_files/figure-html/pcc_q1_brm.png",
          mod_brm |> pp_check(type = "dens_overlay", ndraws = 100) + scale_x_log10() +
            mod_brm |> pp_check(type = "loo_pit_overlay", ndraws = 100),
          width = 10,
          height = 5)
        ## ----end
      }
      ## DHARMa residuals
      {
        ## ---- q1_brm_dharma_a
        load(file = "../data/modelled/mod_q1_brm.RData")
        resids <- make_brms_dharma_res(
          mod_brm,
          integerResponse = FALSE
        )
        wrap_elements(~ testUniformity(resids)) +
          wrap_elements(~ plotResiduals(resids, form = factor(rep(1, nrow(mod_brm$data))))) +
          wrap_elements(~ plotResiduals(resids)) +
          wrap_elements(~ testDispersion(resids))
        ## ----end

        ## ---- q1_brm_dharma
        load(file = "../data/modelled/mod_q1_brm.RData")
        resids <- make_brms_dharma_res(
          mod_brm,
          integerResponse = FALSE
        )
        ggsave(
          filename = "../docs/analysis_files/figure-html/dharma_q1_brm.png",
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
        ## ---- q1_brm_summary
        load(file = "../data/modelled/mod_q1_brm.RData")

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
        save(summ, file = "../data/summ_q1_brm.RData")
        save(summ_tt, file = "../data/summ_tt_q1_brm.RData")
        ## ----end
        ## ---- q1_brm_summary_print
        load(file = "../data/summ_tt_q1_brm.RData")
        summ_tt
        ## ----end
      }
      ## Partial plots -  Type 1. sector/shelf re_formula = NULL
      {
        ## ---- q1_brm_partial_effects
        load(file = "../data/modelled/mod_q1_brm.RData")
        newdata <- brm_generate_newdata_no_dist(data)
        pred <- t(brms::posterior_epred(mod_brm,
          newdata = newdata,
          re_formula = NULL,
          se.fit = TRUE,
          allow_new_levels = TRUE
        ))
        cellmeans_brm <- brm_generate_cellmeans_no_dist(newdata, pred)

        save(cellmeans_brm, file = "../data/modelled/cellmeans_q1_brm_type_1.RData")

        cellmeans_summ_brm <- cellmeans_brm |>
          dplyr::select(-.draws) |>
          dplyr::select(A_SECTOR, SHELF, REPORT_YEAR, Pred) |>
          group_by(A_SECTOR, SHELF, REPORT_YEAR) |>
          summarise_draws(median, HDInterval::hdi) |>
          dplyr::select(A_SECTOR, SHELF, REPORT_YEAR, median, lower, upper)

        save(cellmeans_summ_brm, file = "../data/modelled/cellmeans_summ_q1_brm_type_1.RData")

        p <- brm_partial_plot_no_dist(cellmeans_summ_brm)
        ggsave(
          filename = "../docs/analysis_files/figure-html/partial_q1_brm_type_1.png",
          p,
          width = 15, height = 12
        )
        ## ----end
      }
      ## Partial plots -  Type 2. sector/shelf re_formula = NA
      {
        ## ---- q1_brm_partial_effects_type_2
        load(file = "../data/modelled/mod_q1_brm.RData")
        newdata <- brm_generate_newdata_no_dist(data)
        pred <- t(brms::posterior_epred(mod_brm,
          newdata = newdata,
          re_formula = NA,
          se.fit = TRUE,
          allow_new_levels = TRUE
        ))
        cellmeans_brm <- brm_generate_cellmeans_no_dist(newdata, pred)

        save(cellmeans_brm, file = "../data/modelled/cellmeans_q1_brm_type_2.RData")

        cellmeans_summ_brm <- cellmeans_brm |>
          dplyr::select(-.draws) |>
          dplyr::select(A_SECTOR, SHELF, REPORT_YEAR, Pred) |>
          group_by(A_SECTOR, SHELF, REPORT_YEAR) |>
          summarise_draws(median, HDInterval::hdi) |>
          dplyr::select(A_SECTOR, SHELF, REPORT_YEAR, median, lower, upper)

        save(cellmeans_summ_brm, file = "../data/modelled/cellmeans_summ_q1_brm_type_2.RData")

        p <- brm_partial_plot_no_dist(cellmeans_summ_brm)
        ggsave(
          filename = "../docs/analysis_files/figure-html/partial_q1_brm_type_2.png",
          p,
          width = 15, height = 12
        )
        ## ----end
      }
    }
    
    {
      ## Fit model
      ## Performed on the HPC (23_fit_models_HPC.R)
      {
        if (rerun_models) {
          ## ---- q1_brm_fit
          form <- bf(
            n.points | trials(total.points) ~ 1 +
              (1 | SecShelfYr) +
              (1 | AIMS_REEF_NAME) +
              (1 | zone_depth) +
              (1 | Site) +
              (1 | Transect),
            ## family = "beta_binomial"
            zi ~ (1 | SecShelfYr) + (1| zone_depth) + (1 | Site),
            family = zero_inflated_binomial()
          )
          data |>
            group_by(SecShelfYr) |>
            summarise(
              M = median(qlogis(n.points / total.points), na.rm = TRUE),
              S = mad(qlogis(n.points / total.points), na.rm = TRUE),
              Mean = qlogis(mean(n.points / total.points)),
              SD = qlogis(sd(n.points / total.points))
            )
          priors <- prior(normal(-4, 1), class = "Intercept") +
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
          save(mod_brm, file = "../data/modelled/mod_q1_brm_3.RData")
          ## ----end
        }
      }
      ## Performed on the HPC (23_fit_models_HPC_with_excludes.R)
      {
        if (rerun_models) {
          ## ---- q1_brm_fit_with_excludes
          form <- bf(
            n.points | trials(total.points) ~ 0 + SecShelfYr +
              (1 | SecShelfYr) +
              (1 | AIMS_REEF_NAME) +
              (1 | zone_depth) +
              (1 | Site) +
              (1 | Transect),
            ## family = "beta_binomial"
            zi ~ (1 | SecShelfYr) + (1| zone_depth) + (1 | Site),
            family = zero_inflated_binomial()
          )
          save(data_with_excludes, file = "../data/processed/data_q1_with_excludes.RData")
          data <- dat_with_excludes
          data |>
            group_by(SecShelfYr) |>
            summarise(
              M = median(qlogis(n.points / total.points), na.rm = TRUE),
              S = mad(qlogis(n.points / total.points), na.rm = TRUE),
              Mean = qlogis(mean(n.points / total.points)),
              SD = qlogis(sd(n.points / total.points))
            )
          priors <- prior(normal(0, 3), class = "b") +
          ## priors <- prior(normal(-4, 1), class = "Intercept") +
            ## prior(normal(0, 2), class = "b") +
            prior(student_t(3, 0, 1), class = "sd") +
            prior(logistic(0, 1), class = "Intercept", dpar = "zi") +
            ## prior(normal(0,1), class = 'b', dpar = 'zi') +
            prior(student_t(3, 0, 1), class = "sd", dpar = "zi") 

          for (j in 1:nrow(excludes)) { 
            priors <- priors +   
              prior_string("normal(-4, 0.01)", class = "b",
                coef = paste0("SecShelfYr", gsub(" ", "", excludes[j,'SecShelfYr'])))
          }
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
          save(mod_brm, file = "../data/modelled/mod_q1_brm_3_with_excludes.RData")
          load(file = "../data/modelled/mod_q1_brm_3_with_excludes.RData")
          ## ----end
        }
      }
      ## MCMC diagnostics
      {
        ## ---- q1_brm_trace
        ## load(file = "../data/modelled/mod_q1_brm.RData")
        ## system("scp mlogan@hpc-l001.aims.gov.au:~/Work/AIMS/LTMP/Seriatopora/data/modelled/mod_q1_brm_3.RData ../data/modelled/mod_q1_brm_3.RData")
        ## load(file = "../data/modelled/mod_q1_brm_3.RData")
        load(file = "../data/modelled/mod_q1_brm_3_with_excludes.RData")
        rstan::stan_trace(mod_brm$fit)
        rstan::stan_ac(mod_brm$fit)
        rstan::stan_ess(mod_brm$fit) + rstan::stan_rhat(mod_brm$fit)
        ## ----end
        
        ## ---- q1_brm_trace
        ## load(file = "../data/modelled/mod_q1_brm_3.RData")
        load(file = "../data/modelled/mod_q1_brm_3_with_excludes.RData")
        ggsave(
          filename = "../docs/analysis_files/figure-html/stan_trace_q1_brm.png",
          rstan::stan_trace(mod_brm$fit),
          width = 15, height = 8
        )
        ggsave(
          filename = "../docs/analysis_files/figure-html/stan_ac_q1_brm.png",
          rstan::stan_ac(mod_brm$fit),
          width = 15, height = 8
        )
        ggsave(
          filename = "../docs/analysis_files/figure-html/stan_ess_q1_brm.png",
          rstan::stan_ess(mod_brm$fit) + rstan::stan_rhat(mod_brm$fit),
          width = 15, height = 8
        )
        ## ----end
      }
      ## Model validation
      {
        ## ---- q1_brm_validation_a
        ## load(file = "../data/modelled/mod_q1_brm_3.RData")
        load(file = "../data/modelled/mod_q1_brm_3_with_excludes.RData")
        mod_brm |> pp_check(type = "dens_overlay", ndraws = 100) + scale_x_log10() +
          mod_brm |> pp_check(type = "loo_pit_overlay", ndraws = 100)
        ## ----end

        ## ---- q1_brm_validation
        ## load(file = "../data/modelled/mod_q1_brm_3.RData")
        load(file = "../data/modelled/mod_q1_brm_3_with_excludes.RData")
        ggsave(
          filename = "../docs/analysis_files/figure-html/pcc_q1_brm.png",
          mod_brm |> pp_check(type = "dens_overlay", ndraws = 100) + scale_x_log10() +
            mod_brm |> pp_check(type = "loo_pit_overlay", ndraws = 100),
          width = 10,
          height = 5)
        ## ----end
      }
      ## DHARMa residuals
      {
        ## ---- q1_brm_dharma_a
        load(file = "../data/modelled/mod_q1_brm_3.RData")
        resids <- make_brms_dharma_res(
          mod_brm,
          integerResponse = FALSE
        )
        wrap_elements(~ testUniformity(resids)) +
          wrap_elements(~ plotResiduals(resids, form = factor(rep(1, nrow(mod_brm$data))))) +
          wrap_elements(~ plotResiduals(resids)) +
          wrap_elements(~ testDispersion(resids))
        ## ----end

        ## ---- q1_brm_dharma
        load(file = "../data/modelled/mod_q1_brm_3.RData")
        resids <- make_brms_dharma_res(
          mod_brm,
          integerResponse = FALSE
        )
        ggsave(
          filename = "../docs/analysis_files/figure-html/dharma_q1_brm.png",
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
        ## ---- q1_brm_summary
        ## load(file = "../data/modelled/mod_q1_brm_3.RData")
        load(file = "../data/modelled/mod_q1_brm_3_with_excludes.RData")

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
        save(summ, file = "../data/summ_q1_brm_3.RData")
        save(summ_tt, file = "../data/summ_tt_q1_brm_3.RData")
        ## ----end
        ## ---- q1_brm_summary_print
        load(file = "../data/summ_tt_q1_brm_3.RData")
        summ_tt
        ## ----end
      }
      ## Partial plots -  Type 1. sector/shelf re_formula = NULL
      {
        ## ---- q1_brm_partial_effects
        ## load(file = "../data/modelled/mod_q1_brm_3.RData")
        load(file = "../data/modelled/mod_q1_brm_3_with_excludes.RData")
        newdata <- brm_generate_newdata_no_dist(data)
        pred <- t(brms::posterior_epred(mod_brm,
          newdata = newdata,
          re_formula = NULL,
          se.fit = TRUE,
          allow_new_levels = TRUE
        ))
        cellmeans_brm <- brm_generate_cellmeans_no_dist(newdata, pred)

        save(cellmeans_brm, file = "../data/modelled/cellmeans_q1_brm_type_1_3.RData")

        cellmeans_summ_brm <- cellmeans_brm |>
          dplyr::select(-.draws) |>
          dplyr::select(A_SECTOR, SHELF, REPORT_YEAR, Pred) |>
          group_by(A_SECTOR, SHELF, REPORT_YEAR) |>
          summarise_draws(median, HDInterval::hdi) |>
          dplyr::select(A_SECTOR, SHELF, REPORT_YEAR, median, lower, upper)

        save(cellmeans_summ_brm, file = "../data/modelled/cellmeans_summ_q1_brm_type_1_3.RData")

        p <- brm_partial_plot_no_dist(cellmeans_summ_brm)
        ggsave(
          filename = "../docs/analysis_files/figure-html/partial_q1_brm_type_1_3.png",
          p,
          width = 15, height = 12
        )
        ## ----end
      }
      ## Partial plots -  Type 2. sector/shelf re_formula = NA
      {
        ## ---- q1_brm_partial_effects_type_2
        ## load(file = "../data/modelled/mod_q1_brm_3.RData")
        load(file = "../data/modelled/mod_q1_brm_3_with_excludes.RData")
        newdata <- brm_generate_newdata_no_dist(data)
        pred <- t(brms::posterior_epred(mod_brm,
          newdata = newdata,
          re_formula = ~ (1|SecShelfYr),
          se.fit = TRUE,
          allow_new_levels = TRUE
        ))
        cellmeans_brm <- brm_generate_cellmeans_no_dist(newdata, pred)

        save(cellmeans_brm, file = "../data/modelled/cellmeans_q1_brm_type_2_3.RData")
        load(file = "../data/modelled/cellmeans_q1_brm_type_2_3.RData")

        cellmeans_summ_brm <- cellmeans_brm |>
          dplyr::select(-.draws) |>
          dplyr::select(A_SECTOR, SHELF, REPORT_YEAR, Pred) |>
          group_by(A_SECTOR, SHELF, REPORT_YEAR) |>
          summarise_draws(median, HDInterval::hdi) |>
          dplyr::select(A_SECTOR, SHELF, REPORT_YEAR, median, lower, upper)

        save(cellmeans_summ_brm, file = "../data/modelled/cellmeans_summ_q1_brm_type_2_3.RData")
        load(file = "../data/modelled/cellmeans_summ_q1_brm_type_2_3.RData")

        p <- brm_partial_plot_no_dist(cellmeans_summ_brm)
        p <- p + scale_y_continuous("Seriatopora cover (%)", labels = scales::percent_format(accuracy = 0.1, suffix = "")) +
          scale_x_continuous("Year") +
          theme(
            axis.text.x = element_text(size = 12),
            axis.text.y = element_text(size = 12),
            axis.title.x = element_text(size = 16),
            axis.title.y = element_text(size = 16),
            strip.text = element_text(size = 14)
          )
        p
        ggsave(
          filename = "../docs/analysis_files/figure-html/partial_q1_brm_type_2_3.png",
          p,
          width = 15, height = 13, dpi = 300
        )
        ## ----end
      }
    }
    
  }
}
