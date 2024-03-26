## Retrieve the data
load(file = "../data/processed/data_q2.R")


## Strip out all unnecessary variables

full_data <- data
data <- full_data

## SecShelf, Dist.time and Disturbances (zi full RE) 
{
  ## Prepare data
  {
    data <- full_data

    ## Focus on only the necessary variables
    data <- data |>
      mutate(SecShelf = paste(A_SECTOR, SHELF)) |> 
      dplyr::select(
        n.points, total.points, Dist.time, s, c, d, b, u,
        AIMS_REEF_NAME, Site, Transect, SecShelf, Dist.number,
        SecShelf, SSDist
      ) |>
      filter(!SecShelf %in% c("CG M", "PC M", "PC O")) |>
      mutate(across(c(s, c, d, b, u), \(x) ifelse(Dist.time == "Before", 0, x))) 
  }
  ## Raw means
  {
    ## Cellmeans
    {
      ## ---- 7_1_cellmeans
      cellmeans_summ_raw <-
        data |>
        mutate(Dist.time2 = ifelse(Dist == "Before", "Before", "After")) |> 
        group_by(A_SECTOR, SHELF, Dist, Dist.time) |>
        reframe(
          Mean = c(mean(n.points / total.points)), 
          SD = c(sd(n.points / total.points)), 
          N = c(n(), n())
        ) |>
        mutate(
          lower = Mean - 2 * (SD / sqrt(N)),
          upper = Mean + 2 * (SD / sqrt(N))
        )
      
      save(cellmeans_summ_raw, file = "../data/modelled/cellmeans_summ_raw_7.1.RData")
      cellmeans_summ_raw |>
        ggplot(aes(y = Mean, x = Dist, colour = Dist.time)) +
        geom_pointrange(aes(ymin = lower, ymax = upper)) +
        facet_wrap(A_SECTOR ~ SHELF,
          scales = "free",
          labeller = label_wrap_gen(multi_line = FALSE)
        )
      ## ----end
    }
    ## Effects
    {
      eff_summ_raw <- cellmeans_summ_raw |>
        mutate(Values = Mean) |>
        nest(.by = c(A_SECTOR, SHELF, type)) |>
        mutate(eff = map(
          .x = data,
          .f = ~ {
            before_vs_afters(.x)
          }
        )) |>
        dplyr::select(-data) |> 
        unnest(c(eff))
      save(eff_summ_raw, file = "../data/modelled/eff_summ_raw_7.1.RData")
      eff_summ_raw |> ggplot(aes(y = Dist, x = Values, colour = type)) +
        geom_vline(xintercept = 1, linetype = "dashed") +
        geom_point() +
        facet_grid(A_SECTOR ~ SHELF, scales = "free")
    }
  }
  ## glmmTMB
  {
    ## Fit model
    {
      if (rerun) {
        mod_glmmTMB <- glmmTMB(
          cbind(n.points, total.points - n.points) ~ 1 + (s + c + d + b + u) +
            ## (1 | Dist.time:SecShelf:(s + c + d + b + u)) +
            (1 | SecShelf:(s + c + d + b + u)) +
            (1 | AIMS_REEF_NAME) +
            (1 | Site) +
            (1 | Transect),
          ziformula = ~ 1 + (1 | AIMS_REEF_NAME) +
            (1 | Site) +
            (1 | Transect),
          data = data,
          family = "binomial",
          REML = TRUE
        )

        summary(mod_glmmTMB)
        save(mod_glmmTMB, file = "../data/modelled/mod_glmmTMB_7.1.RData")
      }
    }
    ## Partial plots
    {
      load(file = "../data/modelled/mod_glmmTMB_7.1.RData")

      newdata <- crossing(s = c(0, 1), c = c(0, 1), d = c(0, 1), b = c(0, 1), u = c(0, 1)) |>
        crossing(SecShelf = data$SecShelf) |> 
        filter((s == 0 & d == 0 & c == 0 & b == 0 & u == 0) |
                 (s == 1 & d == 0 & c == 0 & b == 0 & u == 0) |
                 (s == 0 & d == 1 & c == 0 & b == 0 & u == 0) |
                 (s == 0 & d == 0 & c == 1 & b == 0 & u == 0) |
                 (s == 0 & d == 0 & c == 0 & b == 1 & u == 0) |
                 (s == 0 & d == 0 & c == 0 & b == 0 & u == 1) 
        ) |> 
        mutate(Dist = case_when(
                s == 1 ~ "s",
                c == 1 ~ "c",
                d == 1 ~ "d",
                b == 1 ~ "b",
                u == 1 ~ "u",
                .default = "Before"
        ), Dist = factor(Dist, levels = c("Before", "s", "c", "d", "b", "u"))) |>
        arrange(Dist) |>
        mutate(AIMS_REEF_NAME = NA, Site = NA, Transect = NA)

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
          )) 

      cellmeans_summ_glmmTMB <- newdata |>
        dplyr::select(SecShelf, Dist, Pred, lower, upper) |>
        separate(SecShelf, into = c("A_SECTOR", "SHELF"), sep = " ") 
      save(cellmeans_summ_glmmTMB, file = "../data/modelled/cellmeans_summ_glmmTMB_7.1.RData")
      cellmeans_summ_glmmTMB |> ggplot(aes(y = Pred, x = Dist)) +
        geom_pointrange(aes(ymin = lower, ymax = upper)) +
        facet_wrap(A_SECTOR ~ SHELF, scales = "free")
    }
    ## Contrasts
    {
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
    }
    ## Fit model
    {
      data_pred <- data |>
        bind_rows(newdata |>
                    mutate(SecShelf =  paste(A_SECTOR, SHELF)) |>
                    dplyr::select(Dist.time, SecShelf, s, b, c, d, u)) |> 
        mutate(SS = paste(Dist.time, SecShelf, s, b, c, d, u))
      i_newdata <- (nrow(data) + 1):nrow(data_pred)
      mod <- inla(n.points ~ Dist.time + (s + b + c + d + u) +
                    f(model = "iid", SS) +
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
      newdata_pred <- newdata |>
        bind_cols(mod$summary.fitted.values[i_newdata, ])
      newdata_pred 
      newdata_pred |> ggplot(aes(y = `0.5quant`, x = Dist, colour = factor(Dist.time))) +
        geom_pointrange(aes(ymin = `0.025quant`, ymax = `0.975quant`), position = position_dodge(width = 0.5)) +
        facet_grid(A_SECTOR ~ SHELF, scales = "free")


      
  draws <- inla.posterior.sample(n=1000, result = mod)
  cellmeans <- newdata |> bind_cols(sapply(draws, function(x) x$latent[i_newdata]))
cellmeans <- cellmeans |>
        pivot_longer(cols = matches("^[\\.]{3}[0-9]*"), names_to = ".draws") |>
        dplyr::select(Dist.time, A_SECTOR, SHELF, Dist, .draws, value) |>
        posterior::as_draws() |>
        mutate(.draw = .draws) |>
        dplyr::select(-.draws) |>
        mutate(value = plogis(value)) |> 
        group_by(Dist.time, A_SECTOR, SHELF, Dist) |>
        summarise_draws(median, HDInterval::hdi)

      cellmeans |> ggplot(aes(y = median, x = Dist, colour = factor(Dist.time))) +
        geom_pointrange(aes(ymin = lower, ymax = upper), position = position_dodge(width = 0.5)) +
        facet_grid(A_SECTOR ~ SHELF, scales = "free")
      
  contents <- mod$misc$configs$contents
      preds <- posterior_predict.inla(mod, newdata = data_pred[i_newdata, ])
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
      zi = ~ 1 + (1 | AIMS_REEF_NAME) + (1 | Site) + (1 | Transect),
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
    save(mod_brm, file = "../data/modelled/mod_brm_7.1.RData")
    load(file = "../data/modelled/mod_brm_7.1.RData")
    load(file = "../data/modelled/mod_brm_6.1.RData")
  }
    ## Partial plot
    {
      load(file = "../data/modelled/mod_brm_5.1.RData")

      newdata <- crossing(SS = data$SS) |>
        separate(SS, into = c("A_SECTOR", "SHELF", "s", "c", "d", "b", "u"), remove = FALSE) |>
        mutate(across(c(s, c, d, b, u), as.numeric)) |> 
        mutate(SecShelf = paste(A_SECTOR, SHELF)) |> 
        mutate(AIMS_REEF_NAME = NA, Site = NA, Transect = NA, Dist.number = NA, total.points = 1) |>
        mutate(Dist.time = factor(ifelse(s == 0 & d == 0 & c == 0 & b == 0 & u == 0, "Before", "After"), levels = c("Before", "After")))
        p <- t(brms::posterior_epred(mod_brm, newdata = newdata, se.fit = TRUE, allow_new_levels = TRUE))
cellmeans_brm <-
          newdata |>
          cbind(p) |> 
          pivot_longer(cols = matches("^[0-9]*$"), names_to = ".draws") |>
          mutate(Dist = case_when(
            s == 1 ~ "s",
            c == 1 ~ "c",
            d == 1 ~ "d",
            b == 1 ~ "b",
            u == 1 ~ "u",
            .default = "Before"
          )) |>
          dplyr::select(-any_of(c(
            "s", "d", "c", "b", "u", "AIMS_REEF_NAME",
            "Dist.time", "Site", "Transect", "Dist.number", "total.points"
          )))

        cellmeans_summ_brm <- cellmeans_brm |>
          dplyr::select(-.draws, -SS, -SecShelf) |> 
          group_by(A_SECTOR, SHELF, Dist) |>
          summarise_draws(median, HDInterval::hdi) |> 
          dplyr::select(Dist, median, lower, upper)
      cellmeans_summ_brm |>
        ggplot(aes(y = median, x = Dist)) +
        geom_pointrange(aes(ymin = lower, ymax = upper), position = position_dodge(width = 0.5)) +
        facet_wrap(A_SECTOR~SHELF, scales = "free")

      
        newdata <- crossing(Dist.time = data$Dist.time) |>
          crossing(s = c(0, 1), d = c(0, 1), c = c(0, 1), b = c(0, 1), u = c(0, 1)) |>
          crossing(SecShelf = data$SecShelf) |> 
          mutate(AIMS_REEF_NAME = NA, Site = NA, Transect = NA, Dist.number =  NA)
        newdata <- newdata |>
          filter((Dist.time == "Before" & s == 0 & d == 0 & c == 0 & b == 0 & u == 0) |
                   (Dist.time == "After" & s == 1 & d == 0 & c == 0 & b == 0 & u == 0) |
                   (Dist.time == "After" & s == 0 & d == 1 & c == 0 & b == 0 & u == 0) |
                   (Dist.time == "After" & s == 0 & d == 0 & c == 1 & b == 0 & u == 0) |
                   (Dist.time == "After" & s == 0 & d == 0 & c == 0 & b == 1 & u == 0) |
                   (Dist.time == "After" & s == 0 & d == 0 & c == 0 & b == 0 & u == 1)) |>
          mutate(total.points = 1)

        p <- t(brms::posterior_epred(mod_brm, newdata = newdata, se.fit = TRUE, allow_new_levels = TRUE))
        cellmeans_brm <-
          newdata |>
          cbind(p) |>
          pivot_longer(cols = matches("^[0-9]*$"), names_to = ".draws") |>
          mutate(Dist = case_when(
            s == 1 ~ "s",
            c == 1 ~ "c",
            d == 1 ~ "d",
            b == 1 ~ "b",
            u == 1 ~ "u",
            .default = "Before"
          )) |>
          dplyr::select(-any_of(c(
            "s", "d", "c", "b", "u", "AIMS_REEF_NAME",
            "Dist.time", "Site", "Transect", "Dist.number", "total.points"
          )))

        cellmeans_summ_brm <- cellmeans_brm |>
          dplyr::select(-.draws) |> 
          group_by(SecShelf, Dist) |>
          summarise_draws(median, HDInterval::hdi) |> 
          dplyr::select(Dist, median, lower, upper)


      cellmeans_summ_brm |>
        ggplot(aes(y = median, x = Dist)) +
        geom_pointrange(aes(ymin = lower, ymax = upper), position = position_dodge(width = 0.5)) +
        facet_wrap(~SecShelf, scales = "free")



      
      newdata <- crossing(s = c(0, 1), c = c(0, 1), d = c(0, 1), b = c(0, 1), u = c(0, 1)) |>
        crossing(SecShelf = data$SecShelf) |> 
        filter((s == 0 & d == 0 & c == 0 & b == 0 & u == 0) |
                 (s == 1 & d == 0 & c == 0 & b == 0 & u == 0) |
                 (s == 0 & d == 1 & c == 0 & b == 0 & u == 0) |
                 (s == 0 & d == 0 & c == 1 & b == 0 & u == 0) |
                 (s == 0 & d == 0 & c == 0 & b == 1 & u == 0) |
                 (s == 0 & d == 0 & c == 0 & b == 0 & u == 1) 
        ) |> 
        mutate(Dist = case_when(
          s == 1 ~ "s",
          c == 1 ~ "c",
          d == 1 ~ "d",
          b == 1 ~ "b",
          u == 1 ~ "u",
          .default = "Before"
        ), Dist = factor(Dist, levels = c("Before", "s", "c", "d", "b", "u"))) |>
        arrange(Dist) |> 
        dplyr::select(-any_of(c("s", "c", "d", "b", "u")))
      Xmat <- cbind(1, contr.treatment(6, contrasts = TRUE))

      fit <- mod_brm |>
        as.data.frame() |>
        dplyr::select(matches("^b_[^z].*", perl = TRUE))

      cellmeans_brm <- cbind(newdata, t(as.matrix(fit) %*% t(Xmat))) |>
        pivot_longer(cols = matches("^[0-9]*$"), names_to = ".draws")
      cellmeans_summ_brm <- cellmeans_brm |> 
        dplyr::select(-.draws) |>
        group_by(SecShelf, Dist) |>
        mutate(value =  plogis(value)) |>
        summarise_draws(median, HDInterval::hdi)

      save(cellmeans_summ_brm, file = "../data/modelled/cellmeans_summ_brm_5.1.RData")
      
      cellmeans_summ_brm |>
        ggplot(aes(y = median, x = Dist)) +
        geom_pointrange(aes(ymin = lower, ymax = upper), position = position_dodge(width = 0.5)) +
        facet_wrap(~SecShelf, scales = "free")
    }
  ## INLA
  {
    ## Prepare data
    {
      data <- full_data

      ## Focus on only the necessary variables
      data <- data |>
        dplyr::select(
          n.points, total.points, Dist.time, s, c, d, b, u,
          AIMS_REEF_NAME, Site, Transect, SecShelf, Dist.number
        ) |>
        mutate(across(c(s, c, d, b, u), \(x) ifelse(Dist.time == "Before", 0, x))) |> 
        mutate(SS = paste(SecShelf, s, b, c, d, u))


      newdata <- crossing(SecShelf = data$SecShelf, s = c(0, 1), c = c(0, 1),
        d = c(0, 1), b = c(0, 1), u = c(0, 1)) |> 
        mutate(AIMS_REEF_NAME = NA, Site = NA, Transect = NA, Dist.number =  NA)
      newdata <- newdata |>
        filter((s == 0 & d == 0 & c == 0 & b == 0 & u == 0) |
                 (s == 1 & d == 0 & c == 0 & b == 0 & u == 0) |
                 (s == 0 & d == 1 & c == 0 & b == 0 & u == 0) |
                 (s == 0 & d == 0 & c == 1 & b == 0 & u == 0) |
                 (s == 0 & d == 0 & c == 0 & b == 1 & u == 0) |
                 (s == 0 & d == 0 & c == 0 & b == 0 & u == 1) 
        ) |> 
        mutate(Dist = case_when(
          s == 1 ~ "s",
          c == 1 ~ "c",
          d == 1 ~ "d",
          b == 1 ~ "b",
          u == 1 ~ "u",
          .default = "Before"
        )) |> 
        mutate(SS = paste(SecShelf, s, b, c, d, u))

      data_pred <- data |>
        bind_rows(newdata |>
                    dplyr::select(SS, s, b, c, d, u)) 
      i_newdata <- (nrow(data) + 1):nrow(data_pred)


      ## Cellmeans
      ## newdata <- data.frame(Dist.time = unique(data$Dist.time)) |>
      ##   crossing(
      ##     SecShelf = data$SecShelf,
      ##     s = data$s, c = data$c, d = data$d, b = data$b, u = data$u
      ##   )

      ## ## Restrict this to only where the sum of the rows is one so
      ## ## that our newdata is one row per disturbance (for Before and
      ## ## After)
      ## newdata <- newdata |>
      ##         rowwise() |>
      ##         filter(sum(c_across(where(is.numeric))) == 1)
      ## ## Further compress this to just a single Before (no disturbances)
      ## ## and single disturbance types for each After
      ## newdata <- newdata |>
      ##   filter(Dist.time == "After") |>
      ##   bind_rows(data.frame(
      ##     Dist.time = "Before",
      ##     s = 0, c = 0, d = 0, b = 0, u = 0
      ##   ) |> crossing(SecShelf = newdata$SecShelf))

      ## data_pred <- data |>
      ##         bind_rows(newdata) |>
      ##         mutate(SecShelf = forcats::fct_relevel(SecShelf, "CA M"))
      ## i_newdata <- (nrow(data) + 1):nrow(data_pred)

    }
    ## Fit model
    {
      if (rerun) {
        mod <- inla(n.points ~ (s + b + c + d + u) +
                      f(model = "iid", SS) +
        mod <- inla(n.points ~ SecShelf:(s+c+b+d+u)+
                      ## f(model = "iid", Dist.number) +
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
        save(mod, file = "../data/modelled/mod_inla_7.1.RData")
      }
    }
    ## Partial plots
    {
      load(file = "../data/modelled/mod_inla_7.1.RData")

cellmeans_inla <- newdata |> cbind(mod$summary.fitted.values[i_newdata, ])
      cellmeans_inla |>
        ggplot(aes(y = `0.5quant`, x = Dist)) +
        geom_pointrange(aes(ymin = `0.025quant`, ymax = `0.975quant`), position = position_dodge(width = 0.5)) +
        facet_wrap(~SecShelf, scales = "free")
      
      draws <- inla.posterior.sample(n=1000, result = mod)
      cellmeans <- newdata |> cbind(sapply(draws, function(x) x$latent[i_newdata]))
      cellmeans_inla <- cellmeans |>
        pivot_longer(cols = matches("^[0-9]*$"), names_to = ".draws") |>
        dplyr::select(SecShelf, Dist, .draws, value) |>
        posterior::as_draws() |>
        mutate(.draw = .draws) |>
        dplyr::select(-.draws) |>
        separate(SecShelf, into = c("A_SECTOR", "SHELF"), sep = " ") |> 
        mutate(value = plogis(value))
      cellmeans_summ_inla <- cellmeans_inla |> 
        group_by(A_SECTOR, SHELF, Dist) |>
        summarise_draws(median, HDInterval::hdi)

      save(cellmeans_summ_inla, file = "../data/modelled/cellmeans_summ_inla_7.1.RData")
      cellmeans_summ_inla |>
        ggplot(aes(y = median, x = Dist)) +
        geom_pointrange(aes(ymin = lower, ymax = upper), position = position_dodge(width = 0.5)) +
        facet_wrap(A_SECTOR ~ SHELF, scales = "free")
    }
    ## Contrasts
    {
      eff_inla <- cellmeans_inla |>
        mutate(Values = value) |>
        nest(.by = c(A_SECTOR, SHELF, .draw)) |>
        mutate(eff = map(
          .x = data,
          .f = ~ before_vs_afters(.x)
        )) |>
        dplyr::select(-data) |>
        unnest(c(eff))

      eff_summ_inla <- eff_inla |> 
        dplyr::select(-.draw) |>
        group_by(A_SECTOR, SHELF, Dist) |> 
        summarise_draws(
          median,
          HDInterval::hdi
        )

      save(eff_summ_inla, file = "../data/modelled/eff_summ_inla_7.1.RData")
      eff_summ_inla |> ggplot(aes(y = Dist, x = median)) +
        geom_vline(xintercept = 1, linetype = "dashed") +
        geom_pointrange(aes(xmin = lower, xmax = upper)) +
        facet_wrap(A_SECTOR ~ SHELF, scales = "free")
    }
    ## Partial plots - old
    {
      if (1 == 2) {
        newdata_pred <- newdata |>
          bind_cols(mod$summary.fitted.values[i_newdata, ])
        newdata_pred 

        newdata_pred <- newdata_pred |>
          mutate(Dist = case_when(
            s == 1 ~ "s",
            c == 1 ~ "c",
            d == 1 ~ "d",
            b == 1 ~ "b",
            u == 1 ~ "u"
          ))
        newdata_pred <- newdata_pred |>
          separate(SecShelf, into = c("A_SECTOR", "SHELF"), remove = FALSE)
        ## newdata_before <- data.frame(Dist = unique(newdata_pred$Dist)) |>
        ##   crossing(
        ##     A_SECTOR = newdata_pred$A_SECTOR,
        ##     SHELF = newdata_pred$SHELF
        ##   ) |>
        ##   bind_cols(newdata_pred |> filter(Dist.time == "Before") |> dplyr::select(-SecShelf, -A_SECTOR, -SHELF, -s, -c, -d, -b, -u, -Dist))
        
        newdata_pred |>
          ## bind_rows(newdata_before) |> 
          separate(SecShelf, into = c("A_SECTOR", "SHELF"), remove = FALSE) |> 
          ggplot(aes(y = `0.5quant`, x = factor(Dist), color = Dist.time)) +
          geom_pointrange(
            aes(
              ymin = `0.025quant`,
              ymax = `0.975quant`
            ),
            position = position_dodge(width = 0.5)
          ) +
          facet_grid(A_SECTOR ~ SHELF)
      }
    }
    ## Contrasts - old
    {
      if (1 == 2) {
        nd.pred <- nd1 |>
          bind_cols(mod$summary.lincomb.derived) |>
          mutate(across(where(is.numeric), exp))
        nd.pred |> ggplot(aes(x = `0.5quant`, y = Dist)) +
          geom_vline(xintercept = 1, linetype = "dashed") +
          geom_pointrange(aes(xmin = `0.025quant`, xmax = `0.975quant`)) +
          scale_x_continuous("Effect (Before - After) on a fold scale",
            trans = "log2", breaks = scales::extended_breaks(n = 8)
          )
      }
    }
  }
  ## Comparisons
  {
    
    load(file = "../data/modelled/cellmeans_summ_raw_7.1.RData")
    load(file = "../data/modelled/cellmeans_summ_glmmTMB_7.1.RData")
    load(file = "../data/modelled/cellmeans_summ_brm_7.1.RData")
    load(file = "../data/modelled/cellmeans_summ_inla_7.1.RData")

    cellmeans_summ <- bind_rows(
      cellmeans_summ_raw |>
        mutate(method = "raw") |>
        dplyr::rename(median = Mean, lower = lower, upper = upper),
      cellmeans_summ_glmmTMB |>
        mutate(method = "glmmTMB", type = "median") |>
        dplyr::rename(median = Pred, lower = lower, upper = upper),
      cellmeans_summ_brm |> mutate(method = "brm", type = "median") |>
        dplyr::rename(median = median, lower = lower, upper = upper),
      cellmeans_summ_inla |>
        mutate(method = "inla", type = "median") |> 
        dplyr::rename(median = median, lower = lower, upper = upper)
    )
    cellmeans_summ

    cellmeans_summ |>
      ggplot(aes(y = median, x = method, shape = type, colour = Dist)) +
      geom_pointrange(aes(ymin = lower, ymax = upper),
        position = position_dodge(width = 0.5)
      ) +
      facet_wrap(A_SECTOR ~ SHELF, scales = "free")

    load(file = "../data/modelled/eff_summ_raw_7.1.RData")
    load(file = "../data/modelled/eff_summ_glmmTMB_7.1.RData")
    load(file = "../data/modelled/eff_summ_brm_7.1.RData")
    load(file = "../data/modelled/eff_summ_inla_7.1.RData")

    eff_summ <- bind_rows(
      eff_summ_raw |>
        mutate(method = "raw") |>
        dplyr::rename(median = Values),
      eff_summ_glmmTMB |>
        mutate(method = "glmmTMB", type = "median") |>
        dplyr::rename(median = Values),
      eff_summ_brm |>
        mutate(method = "brm", type = "median") |>
        dplyr::rename(median = median, lower = lower, upper = upper),
      eff_summ_inla |>
        mutate(method = "inla", type = "median") |> 
        dplyr::rename(median = median, lower = lower, upper = upper)
    )
    eff_summ

    eff_summ |>
      ggplot(aes(x = median, y = Dist, colour = method, shape = type)) +
      geom_vline(xintercept = 1, linetype = "dashed") +
      geom_pointrange(aes(xmin = lower, xmax = upper),
        position = position_dodge(width = 0.5)) +
      scale_x_continuous("Effect (Before - After) on a fold scale",
        trans = "log2", breaks = scales::extended_breaks(n = 8)
      ) +
      facet_wrap(A_SECTOR ~ SHELF, scales = "free")
  }
}
