rerun_models <- FALSE

## Model 1.1  logit(P()) = β₀ + γ₀ᵢ + γ₁ⱼ + γ₂ₖ  (Binomial)
## source("33_fit_models_1.1.R")

## Model 1.2  logit(P()) = β₀ + γ₀ᵢ + γ₁ⱼ + γ₂ₖ, logit(P()) = α₀  (Zero-inflated Binomial)
## source("33_fit_models_1.2.R")

## Model 1.3  logit(P()) = β₀ + γ₀ᵢ + γ₁ⱼ + γ₂ₖ, logit(P()) = α₀ + γ₀ᵢᵣ + γ₁ⱼᵣ + γ₂ₖᵣ (Zero-inflated Binomial)
## source("33_fit_models_1.3.R")

## Model 2.1  logit(P()) = β₀ + β₁*Dist.time + γ₀ᵢ + γ₁ⱼ + γ₂ₖ, logit(P()) = α₀ + γ₀ᵢᵣ + γ₁ⱼᵣ + γ₂ₖᵣ (Zero-inflated Binomial)
## source("33_fit_models_2.1.R")

## Model 3.1  logit(P()) = β₀ + β₁*Dist.time + β₃₋₁₁*SecShelf + γ₀ᵢ + γ₁ⱼ + γ₂ₖ, logit(P()) = α₀ + γ₀ᵢᵣ + γ₁ⱼᵣ + γ₂ₖᵣ (Zero-inflated Binomial)
## source("33_fit_models_3.1.R")

## Model 4.1  logit(P()) = β₀ + β₁*Dist.time + γ₀ᵢ + γ₁ⱼ + γ₂ₖ, logit(P()) = α₀  (Zero-inflated Binomial)
## source("33_fit_models_4.1.R")

## Model 5.1  logit(P()) = β₀ + β₁*Dist.time + γ₀ᵢ + γ₁ⱼ + γ₂ₖ, logit(P()) = α₀ + γ₀ᵢᵣ + γ₁ⱼᵣ + γ₂ₖᵣ  (Zero-inflated Binomial)
## source("33_fit_models_5.1.R")

## Model 6.1  logit(P()) = β₀ +  γₐᵢ + γ₀ᵢ + γ₁ⱼ + γ₂ₖ, logit(P()) = α₀ (Zero-inflated Binomial)
## source("33_fit_models_6.1.R")

## Model 7.1  logit(P()) = β₀ +  γₐᵢ + γ₀ᵢ + γ₁ⱼ + γ₂ₖ, logit(P()) = α₀ + γ₀ᵢᵣ + γ₁ⱼᵣ + γ₂ₖᵣ  (Zero-inflated Binomial)
## source("33_fit_models_7.1.R")
## source("33_fit_models_7.1_raw.R")
## source("33_fit_models_1.1.R")


source("33_fit_models_1.5.R")
source("33_fit_models_1.6.R")

## data <- data |>
##   dplyr::select(n.points, total.points, A_SECTOR, SHELF, SecShelf,
##     AIMS_REEF_NAME, Site, Transect, Dist.time)


##     ## Prepare data
##     {
##       newdata <- data.frame(Dist.time = factor(c("Before", "After"),
##         levels = c("Before", "After"))) |>
##         crossing(A_SECTOR = data$A_SECTOR, SHELF = data$SHELF) |>
##         mutate(SecShelf = factor(paste(A_SECTOR, SHELF)))

##       data.pred <- data |> bind_rows(newdata)
##       i.newdata <- (nrow(data) + 1):nrow(data.pred)
      
##     }












## SecShelf, Dist.time and Disturbances (zi intercept only) 
{
  ## Prepare data
  {
    data <- full_data

    ## Focus on only the necessary variables
    data <- data |>
      dplyr::select(
        n.points, total.points, Dist.time, s, c, d, b, u,
        AIMS_REEF_NAME, Site, Transect, SecShelf, Dist.number,
        SecShelf
      ) |>
      filter(!SecShelf %in% c("CG M", "PC M", "PC O")) |>
              mutate(across(c(s, c, d, b, u), \(x) ifelse(Dist.time == "Before", 0, x))) 
  }
  ## Raw means
  {
    data_summ_s <- data |>
      dplyr::select(n.points, total.points, Dist.time, SecShelf, s, c, d, b, u) |>
      group_by(Dist.time, SecShelf, s) |>
      summarise(
        Median = median(n.points / total.points),
        Mean = mean(n.points / total.points),
        SD = sd(n.points / total.points),
        N = n()
      ) |>
      mutate(
        lower = Mean - 2 * (SD / N),
        upper = Mean + 2 * (SD / N)
      ) |>
      separate(SecShelf, into = c("A_SECTOR", "SHELF"), sep = " ") |>
              filter((Dist.time == "Before" & s == 0) | Dist.time == "After" & s == 1)

    data_summ_c <- data |>
      dplyr::select(n.points, total.points, Dist.time, SecShelf, s, c, d, b, u) |>
      group_by(Dist.time, SecShelf, c) |>
      summarise(
        Median = median(n.points / total.points),
        Mean = mean(n.points / total.points),
        SD = sd(n.points / total.points),
        N = n()
      ) |>
      mutate(
        lower = Mean - 2 * (SD / N),
        upper = Mean + 2 * (SD / N)
      ) |>
      separate(SecShelf, into = c("A_SECTOR", "SHELF"), sep = " ") |>
              filter((Dist.time == "Before" & c == 0) | Dist.time == "After" & c == 1)

    data_summ_d <- data |>
      dplyr::select(n.points, total.points, Dist.time, SecShelf, s, c, d, b, u) |>
      group_by(Dist.time, SecShelf, d) |>
      summarise(
        Median = median(n.points / total.points),
        Mean = mean(n.points / total.points),
        SD = sd(n.points / total.points),
        N = n()
      ) |>
      mutate(
        lower = Mean - 2 * (SD / N),
        upper = Mean + 2 * (SD / N)
      ) |>
      separate(SecShelf, into = c("A_SECTOR", "SHELF"), sep = " ") |>
              filter((Dist.time == "Before" & d == 0) | Dist.time == "After" & d == 1)

    data_summ_b <- data |>
      dplyr::select(n.points, total.points, Dist.time, SecShelf, s, c, d, b, u) |>
      group_by(Dist.time, SecShelf, b) |>
      summarise(
        Median = median(n.points / total.points),
        Mean = mean(n.points / total.points),
        SD = sd(n.points / total.points),
        N = n()
      ) |>
      mutate(
        lower = Mean - 2 * (SD / N),
        upper = Mean + 2 * (SD / N)
      ) |>
      separate(SecShelf, into = c("A_SECTOR", "SHELF"), sep = " ") |>
              filter((Dist.time == "Before" & b == 0) | Dist.time == "After" & b == 1)
    
    data_summ_u <- data |>
      dplyr::select(n.points, total.points, Dist.time, SecShelf, s, c, d, b, u) |>
      group_by(Dist.time, SecShelf, u) |>
      summarise(
        Median = median(n.points / total.points),
        Mean = mean(n.points / total.points),
        SD = sd(n.points / total.points),
        N = n()
      ) |>
      mutate(
        lower = Mean - 2 * (SD / N),
        upper = Mean + 2 * (SD / N)
      ) |>
      separate(SecShelf, into = c("A_SECTOR", "SHELF"), sep = " ") |>
              filter((Dist.time == "Before" & u == 0) | Dist.time == "After" & u == 1)

    data_summ <- bind_rows(data_summ_s, data_summ_c, data_summ_d,
      data_summ_b, data_summ_u) |> 
      mutate(Dist = case_when(
        s == 1 ~ "s",
        c == 1 ~ "c",
        d == 1 ~ "d",
        b == 1 ~ "b",
        u == 1 ~ "u",
        .default = "Before"
      ))
      
    ## data_summ |>
    ##   ggplot(aes(y = Mean, x = Dist.time)) +
    ##   geom_pointrange(aes(ymin = lower, ymax = upper)) +
    ##   facet_grid(A_SECTOR ~ SHELF, scales = "free")

    data_summ |> ggplot(aes(y = Mean, x = Dist, colour = factor(Dist.time))) +
      geom_pointrange(aes(ymin = lower, ymax = upper), position = position_dodge(width = 0.5)) +
      facet_grid(A_SECTOR ~ SHELF, scales = "free")
  
  }
  ## glmmTMB
  {
    ## Fit model
    {
      mod_glmmTMB <- glmmTMB(cbind(n.points, total.points-n.points) ~ 1+
                               (s + c + d + b + u) +
                               (1 | Dist.time:SecShelf:(s + c + d + b + u)) +
                               (1|AIMS_REEF_NAME) +
                               (1|Site) +
                               (1|Transect),
        ziformula = ~1 + (1|AIMS_REEF_NAME) +
                               (1|Site) +
                               (1|Transect),
        data = data,
        family = "binomial", 
        REML = TRUE
      )

      summary(mod_glmmTMB)

    }
    ## DHARMa
    {
      resids <- simulateResiduals(mod_glmmTMB, plot = TRUE, integerResponse = TRUE)
      wrap_elements(~testUniformity(resids)) +
        wrap_elements(~plotResiduals(resids)) +
        wrap_elements(~testDispersion(resids)) 
      testDispersion(resids)
      testZeroInflation(resids)
    }
    ## Partial plots
    {
      ## Since this model operates mainly through random effects, we cannot
      ## use emmeans to produce the cell means or contrasts
      ## Instead we just have to perform predictions.
      
      newdata <- crossing(Dist.time = data$Dist.time) |>
        crossing(s = c(0,1), d = c(0, 1), c = c(0, 1), b = c(0, 1), u = c(0, 1),
          SecShelf = as.character(unique(data$SecShelf))) |>
        mutate(AIMS_REEF_NAME = NA, Site = NA, Transect = NA, Dist.number =  NA)
      newdata <- newdata |>
        filter((Dist.time == "Before" & s == 0 & d == 0 & c == 0 & b == 0 & u == 0) |
                 (Dist.time == "After" & s == 1 & d == 0 & c == 0 & b == 0 & u == 0) |
                 (Dist.time == "After" & s == 0 & d == 1 & c == 0 & b == 0 & u == 0) |
                 (Dist.time == "After" & s == 0 & d == 0 & c == 1 & b == 0 & u == 0) |
                 (Dist.time == "After" & s == 0 & d == 0 & c == 0 & b == 1 & u == 0) |
                 (Dist.time == "After" & s == 0 & d == 0 & c == 0 & b == 0 & u == 1) 
        )

      p <- predict(mod_glmmTMB, newdata = newdata, se.fit = TRUE)
      newdata <- newdata |>
        bind_cols(fit = p$fit, se = p$se.fit) |>
        mutate(Pred = plogis(fit), lower = plogis(fit - 2 * se), upper = plogis(fit + 2 * se))
      newdata <- newdata |>
        mutate(Dist = case_when(
          s == 1 ~ "s",
          c == 1 ~ "c",
          d == 1 ~ "d",
          b == 1 ~ "b",
          u == 1 ~ "u",
          .default = "Before"
        )) |>
        separate(SecShelf, into = c("A_SECTOR", "SHELF"), sep = " ")

      newdata |> dplyr::select(Dist.time, Dist, A_SECTOR, SHELF, Pred, lower, upper)
      
      newdata |>
        ggplot(aes(y = Pred, x = Dist, colour = factor(Dist.time))) +
        geom_pointrange(aes(ymin = lower, ymax = upper),
          position = position_dodge(width = 0.5)) +
        facet_grid(A_SECTOR ~ SHELF, scales = "free")
    }
    ## Contrasts
    {
      ## These are not really possible - at least not with confidence intervals
    }
    
  }
  ## brms
  {
    form <- bf(
      n.points | trials(total.points) ~ 1+ (s + c + d + b + u) + 
        (1 | Dist.time:SecShelf:(s + c + d + b + u)) +
        (1 | AIMS_REEF_NAME) +
        (1 | Site) +
        (1 | Transect),
      zi = ~ 1 +,
      family = "zero_inflated_binomial"
    )
    priors <- prior(normal(0, 1), class = "Intercept") +
            prior(normal(0, 1), class = "b") +
            prior(student_t(3, 0, 1), class = "sd") +
            prior(logistic(0, 1), class = "Intercept", dpar = "zi")
    mod_brm <- brm(form,
      data = data,
      prior = priors,
      iter = 5000, warmup = 1000,
      chains = 3, cores = 3,
      thin = 4,
      backend = "cmdstanr",
      control = list(adapt_delta = 0.99),
      silent =  0#,
      ## refresh = 100
    )
    summary(mod_brm)
    save(mod_brm, file = "../data/modelled/mod_brm_6.1.RData")
    load(file = "../data/modelled/mod_brm_6.1.RData")
  }
  ## INLA
  {
    ## Prepare data OLD
    {
      if (1 != 2) {
        data <- full_data

        ## Focus on only the necessary variables
        data <- data |>
          dplyr::select(
            n.points, total.points, Dist.time, s, c, d, b, u,
            AIMS_REEF_NAME, Site, Transect, SecShelf, Dist.number
          ) |>
          mutate(across(c(s, c, d, b, u), \(x) ifelse(Dist.time == "Before", 0, x)))
        ## Cellmeans
        newdata <- data.frame(Dist.time = unique(data$Dist.time)) |>
          crossing(
            SecShelf = data$SecShelf,
            s = data$s, c = data$c, d = data$d, b = data$b, u = data$u
          )

        ## Restrict this to only where the sum of the rows is one so
        ## that our newdata is one row per disturbance (for Before and
        ## After)
        newdata <- newdata |>
          rowwise() |>
          filter(sum(c_across(where(is.numeric))) == 1)
        ## Further compress this to just a single Before (no disturbances)
        ## and single disturbance types for each After
        newdata <- newdata |>
          filter(Dist.time == "After") |>
          bind_rows(data.frame(
            Dist.time = "Before",
            s = 0, c = 0, d = 0, b = 0, u = 0
          ) |> crossing(SecShelf = newdata$SecShelf))

        data_pred <- data |>
          bind_rows(newdata) |>
          mutate(SecShelf = forcats::fct_relevel(SecShelf, "CA M"))
        i_newdata <- (nrow(data) + 1):nrow(data_pred)

        ## Pre-defined contrasts
        ## Compare Afters (for each disturbance) to Before (no disturbance)
        ## for each sector/shelf
        Xmat <- model.matrix(~ SecShelf*Dist.time * (s + b + c + d + u), data = newdata)

        nd <- newdata |>
          mutate(Dist = case_when(s == 1 ~ "s",
            c == 1 ~ "c",
            d == 1 ~ "d",
            b == 1 ~ "b",
            u == 1 ~ "u")) |>
          dplyr::select(-any_of(c("s", "c", "d", "b", "u"))) |>
          bind_cols(Xmat) |>
          filter(!is.na(Dist)) |>
          mutate(`(Intercept)` = 0) |>
          dplyr::select(-Dist.time, -Dist, -SecShelf) |>
          as.matrix()

        lincomb <- inla.make.lincombs(as.data.frame(nd))

        nd1 <- newdata |>
          mutate(Dist = case_when(s == 1 ~ "s",
            c == 1 ~ "c",
            d == 1 ~ "d",
            b == 1 ~ "b",
            u == 1 ~ "u")) |>
          dplyr::select(-any_of(c("s", "c", "d", "b", "u"))) |> 
          filter(!is.na(Dist)) 
      }
    }
    ## Prepare data
    {
      newdata <- crossing(Dist.time = data$Dist.time) |>
        crossing(s = c(0,1), d = c(0, 1), c = c(0, 1), b = c(0, 1), u = c(0, 1),
          SecShelf = as.character(unique(data$SecShelf))) |>
        separate(SecShelf, into = c("A_SECTOR", "SHELF"), sep = " ", remove = FALSE) |> 
        mutate(AIMS_REEF_NAME = NA, Site = NA, Transect = NA, Dist.number =  NA)
      newdata <- newdata |>
        filter((Dist.time == "Before" & s == 0 & d == 0 & c == 0 & b == 0 & u == 0) |
                 (Dist.time == "After" & s == 1 & d == 0 & c == 0 & b == 0 & u == 0) |
                 (Dist.time == "After" & s == 0 & d == 1 & c == 0 & b == 0 & u == 0) |
                 (Dist.time == "After" & s == 0 & d == 0 & c == 1 & b == 0 & u == 0) |
                 (Dist.time == "After" & s == 0 & d == 0 & c == 0 & b == 1 & u == 0) |
                 (Dist.time == "After" & s == 0 & d == 0 & c == 0 & b == 0 & u == 1) 
        )

      data_pred <- data |>
        bind_rows(newdata |>
                    dplyr::select(Dist.time, SecShelf, s, b, c, d, u)) |> 
        mutate(SS = paste(Dist.time, SecShelf, s, b, c, d, u))
      i_newdata <- (nrow(data) + 1):nrow(data_pred)
    }
    ## Fit model
    {
      mod <- inla(n.points ~ Dist.time + (s + b + c + d + u) +
                    f(model = "iid", SS) +
                    f(model = "iid", AIMS_REEF_NAME) +
                    f(model = "iid", Site) +
                    f(model = "iid", Transect),
        data = data_pred,
        Ntrials = data_pred$total.points,
        family = "zeroinflatedbinomial1", #"binomial" 
        control.predictor = list(link = 1, compute = TRUE),
        control.compute = list(config = TRUE,
          dic = TRUE, waic = TRUE, cpo = TRUE
        )
      )
      summary(mod)
      ## autoplot(mod)
    }
    ## Fit model
    {
      mod <- inla(n.points ~ SecShelf * Dist.time * (s + b + c + d + u) +
                    f(model = "iid", Dist.number) +
                    f(model = "iid", AIMS_REEF_NAME) +
                    f(model = "iid", Site) +
                    f(model = "iid", Transect),
        data = data_pred,
        Ntrials = data_pred$total.points,
        ## lincomb =  lincomb,
        family = "zeroinflatedbinomial1", #"binomial" 
        control.predictor = list(link = 1, compute = TRUE),
        control.compute = list(config = TRUE,
          dic = TRUE, waic = TRUE, cpo = TRUE
        )
      )
      summary(mod)
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
    ## Partial plots
    {
      newdata_pred <- newdata |>
        bind_cols(mod$summary.fitted.values[i_newdata, ]) |> 
        mutate(Dist = case_when(
          s == 1 ~ "s",
          c == 1 ~ "c",
          d == 1 ~ "d",
          b == 1 ~ "b",
          u == 1 ~ "u",
          .default = "Before"
        )) |>
        separate(SecShelf, into = c("A_SECTOR", "SHELF"), sep = " ")
      newdata_pred 
      
      newdata_pred |>
        ggplot(aes(y = `0.5quant`, x = Dist, colour = factor(Dist.time))) +
        geom_pointrange(aes(ymin = `0.025quant`, ymax = `0.975quant`),
          position = position_dodge(width = 0.5)) +
        facet_grid(A_SECTOR ~ SHELF, scales = "free")


      newdata_pred |> ggplot(aes(y = `0.5quant`, x = Dist, colour = factor(Dist.time))) +
        geom_pointrange(aes(ymin = `0.025quant`, ymax = `0.975quant`, shape = "INLA"), position = position_dodge(width = 0.5)) +
        geom_pointrange(data = data_summ, aes(y = Mean, ymin = lower, ymax = upper, shape = "Raw"), position = position_nudge(x = 0.2)) +
        geom_pointrange(data =  newdata, aes(y =  Pred, ymin = lower, ymax = upper, shape = "glmmTMB"), position = position_nudge(x = -0.2)) +
        facet_wrap(A_SECTOR ~ SHELF, scales = "free")

      
      ## newdata_pred <- newdata_pred |>
      ##   mutate(Dist = case_when(
      ##     s == 1 ~ "s",
      ##     c == 1 ~ "c",
      ##     d == 1 ~ "d",
      ##     b == 1 ~ "b",
      ##     u == 1 ~ "u"
      ##   ))
      ## newdata_pred <- newdata_pred |>
      ##   separate(SecShelf, into = c("A_SECTOR", "SHELF"), remove = FALSE)
      ## ## newdata_before <- data.frame(Dist = unique(newdata_pred$Dist)) |>
      ## ##   crossing(
      ## ##     A_SECTOR = newdata_pred$A_SECTOR,
      ## ##     SHELF = newdata_pred$SHELF
      ## ##   ) |>
      ## ##   bind_cols(newdata_pred |> filter(Dist.time == "Before") |> dplyr::select(-SecShelf, -A_SECTOR, -SHELF, -s, -c, -d, -b, -u, -Dist))
      
      ##   newdata_pred |>
      ##           ## bind_rows(newdata_before) |> 
      ##         separate(SecShelf, into = c("A_SECTOR", "SHELF"), remove = FALSE) |> 
      ##         ggplot(aes(y = `0.5quant`, x = factor(Dist), color = Dist.time)) +
      ##         geom_pointrange(
      ##                 aes(
      ##                         ymin = `0.025quant`,
      ##                         ymax = `0.975quant`
      ##                 ),
      ##                 position = position_dodge(width = 0.5)
      ##         ) +
      ## facet_grid(A_SECTOR ~ SHELF)
    }
    ## Partial plots - version 2
    {
      draws <- inla.posterior.sample(n=1000, result = mod)
      cellmeans <- newdata |>
        mutate(Dist = case_when(
          s == 1 ~ "s",
          c == 1 ~ "c",
          d == 1 ~ "d",
          b == 1 ~ "b",
          u == 1 ~ "u",
          .default = "Before"
        )) |>
        separate(SecShelf, into = c("A_SECTOR", "SHELF"), sep = " ") |> 
        bind_cols(sapply(draws, function(x) x$latent[i_newdata]))

      cellmeans <- cellmeans |>
              pivot_longer(cols = matches("^[\\.]{3}[0-9]*"), names_to = ".draws") |>
              dplyr::select(Dist.time, A_SECTOR, SHELF, Dist, .draws, value) |>
              posterior::as_draws() |>
              mutate(.draw = .draws) |>
        mutate(value = plogis(value)) |> 
        dplyr::select(-.draws)
      cellmeans_summ <- cellmeans |> 
        group_by(Dist.time, A_SECTOR, SHELF, Dist) |>
        summarise_draws(median, HDInterval::hdi)

      cellmeans_summ |>
        ggplot(aes(y = median, x = Dist, colour = factor(Dist.time))) +
        geom_pointrange(aes(ymin = lower, ymax = upper),
          position = position_dodge(width = 0.5)) +
        facet_grid(A_SECTOR ~ SHELF, scales = "free")
    }
    ## Contrasts - version 2
    {
      eff <-
              cellmeans |>
        ## filter(.draw == "...14", A_SECTOR == "CA", SHELF == "I") |>
              nest(.by = c(A_SECTOR, SHELF, .draw)) |>
              mutate(data1 = map(
                      .x = data,
                      .f = ~ {
                              .x <- .x |>
                                      mutate(Dist = factor(Dist,
                                              levels = c("Before", "s", "c", "d", "b", "u")
                                      )) |>
                                      arrange(Dist)
                              xmat <- cbind(-1, 1 * contr.treatment(6, base = 1, contrast = TRUE))
                              xmat <- xmat[-1, ]
                              ## print(.x)
                              ## print(xmat)
                              x <- log(as.vector(as.vector(.x$value)))
                              ## print(x)
                              ## print(as.vector(x %*% t(xmat)))
                              ## print(exp(as.vector(x %*% t(xmat))))
                              data.frame(
                                      Dist = .x$Dist[-1],
                                      Values = exp(as.vector(x %*% t(xmat)))
                              )
                      }, .progress = TRUE
              ))

      eff_summ <- eff |> dplyr::select(-data) |> 
                unnest(c(data1)) |>
                ungroup() |>
                group_by(A_SECTOR, SHELF, Dist) |> 
        posterior::summarise_draws(
          median,
          HDInterval::hdi,
          Pl = ~ mean(.x < 1),
          Pg = ~ mean(.x > 1)
        )


      eff_summ |>
        ggplot(aes(x = median, y = Dist)) +
        geom_vline(xintercept = 1, linetype = "dashed") +
        geom_pointrange(aes(xmin = lower, xmax = upper)) +
        scale_x_continuous("Effect (Before - After) on a fold scale",
          trans = "log2",
          breaks = c(0.5, 1, 2)
          ## breaks = scales::log_breaks(n = 8, base = 2)
          ## breaks = scales::pretty_breaks(n = 8)
          ## breaks = scales::extended_breaks(n = 8)
        ) +
        facet_grid(A_SECTOR ~ SHELF, scales = "free")
    }
  }
}



#####


