
## ---- setup
assignInNamespace('.sep.label',  "^\\ *(#|--)+\\s*(@knitr|----+)(.*?)-*\\s*$", ns='knitr')

tidyverse_style_with_comments_removed <- function() {
  remove_comments <- function(pd) {
    is_comment <- pd$token == "COMMENT"
    pd <- pd[!is_comment,]
    pd
  }
  tidyverse_with_comments_removed <- styler::tidyverse_style()
  tidyverse_with_comments_removed$token$remove_comments <- remove_comments
  tidyverse_with_comments_removed
}

knitr::opts_chunk$set(echo = TRUE, message=FALSE, warning=FALSE,cache.lazy = FALSE, tidy='styler',
                      tidy.opts=list(transformers = tidyverse_style_with_comments_removed()))
## knitr::opts_chunk$set(echo = TRUE, warning=FALSE, message=FALSE)

# save the built-in output hook
hook_output <- knitr::knit_hooks$get("output")

# set a new output hook to truncate text output
knitr::knit_hooks$set(output = function(x, options) {
  if (!is.null(n <- options$out.lines)) {
    x <- xfun::split_lines(x)
    if (length(x) > n) {
      # truncate the output
      x <- c(head(x, n), "....\n")
    }
    x <- paste(x, collapse = "\n")
  }
  hook_output(x, options)
})
options(tinytex.engine = 'xelatex')

## ----end

## ---- before_after_tests_function
before_after_tests <- function(.x) {
  .y <- .x |> mutate(before_flag = NA, after_flag = NA, n_flag = NA, cover_flag = NA)
  before <- .x |>
    filter(Dist.time == "Before") |>
    distinct()
  if (nrow(before) < 1) {
    .y <- .y |> mutate(before_flag = "no before")
  }
  if (nrow(before) > 1) {
    .y <- .y |> mutate(before_flag = "multiple before")
  }
  not_before <- .x |> filter(Dist.time == "After")
  if ("DISTURBANCE_TYPE" %in% colnames(not_before)) {
    if (any(not_before$DISTURBANCE_TYPE == "n")) {
      .y <- .y |> mutate(n_flag = "after disturb n")
    }
  }
  if (nrow(not_before) < 1) {
    .y <- .y |> mutate(after_flag = "no after")
  } else {
    ## no_before <- no_before |> mutate(
    .y <- .y |> mutate(
      before_cover = mean(before$COVER),
      cover_diff = COVER - before_cover
    )
    .y |> mutate(cover_flag = ifelse(cover_diff > 0, "cover increased", NA))
  }
}

## ----end


## ---- model_summary_functions
model_summary <- function(mod, ...) {
   UseMethod("model_summary")
}

model_summary.glmmTMB <- function(mod, exponentiate = FALSE) {
  mod |>
    parameters::model_parameters(exponentiate = exponentiate) |>
    as.data.frame() |> 
    ## sanitise
    dplyr::rename(
      term = Parameter,
      estimate = Coefficient,
      conf.low = CI_low,
      conf.high = CI_high,
      effect = Effects
    ) |>
    mutate(term = str_replace(term, "SD \\((.*)\\)", paste0("SD (Intercept ", Group, ")"))) |>
      mutate(term = ifelse(Component == "dispersion", "phi", term)) |>
      arrange(effect) |>
      dplyr::select(term, estimate, conf.low, conf.high, effect)
  }
model_summary.brmsfit <- function(mod, exponentiate = FALSE) {
  mod |>
    tidybayes::tidy_draws() |>
    dplyr::select(matches("^b_|^sd_|^phi")) |>
    {\(.)
      if (exponentiate) mutate(., across(matches("^b_"), \(x) exp(x)))
      else .
    }() |> 
    tidybayes::summarise_draws(
      median,
      HDInterval::hdi
    ) |>
    ## sanitise
    dplyr::rename(
      term = variable,
      estimate = median,
      conf.low = lower,
      conf.high = upper
    ) |>
    mutate(effect = ifelse(str_detect(term, "sd_"), "random", "fixed")) |>
    mutate(term = str_replace(term, "b_Intercept", "(Intercept)")) |> 
    mutate(term = str_replace(term, "b_(.*)", "\\1")) |> 
    mutate(term = str_replace(term, "sd_(.*)__Intercept", "SD (Intercept \\1)")) |> 
    arrange(effect) |> 
      dplyr::select(term, estimate, conf.low, conf.high, effect)
}
model_summary.inla <- function(mod, draws, exponentiate = FALSE) {

  contents <- mod$misc$configs$contents

  form <- mod$.args$formula
  gf <- INLA:::inla.interpret.formula(form)

  fixed_terms <- "(Intercept)"
  i_params <- contents$start[match(fixed_terms, contents$tag)]
  fixed_params <- t(sapply(draws, function(x) x$latent[i_params])) |>
    matrix(ncol = gf$n.fix) |>
    as.data.frame() |>
    setNames(fixed_terms) 
  if (exponentiate) fixed_params <- exp(fixed_params)

  random_terms <- sapply(mod$all.hyper$random, function(x) x$hyperid)
  i_params <- str_which(names(draws[[1]]$hyperpar), paste(random_terms, collapse = "|"))
  random_params <- t(sapply(draws, function(x) x$hyperpar[i_params])) |>
    matrix(ncol = gf$n.random) |>
    as.data.frame() |>
    setNames(paste0("SD (Intercept ", random_terms, ")")) |>
    mutate(across(everything(), \(x) 1 / sqrt(x))) 

  family_terms <- sapply(mod$all.hyper$family, function(x) x$label)
  family_terms <- str_subset(names(draws[[1]]$hyperpar), family_terms) |>
    str_replace("(\\w+).*", "\\1")
  i_params <- str_which(names(draws[[1]]$hyperpar), paste(family_terms, collapse = "|"))
  family_params <- t(sapply(draws, function(x) x$hyperpar[i_params])) |>
    matrix(ncol = length(family_terms)) |>
    as.data.frame() |>
    setNames(str_replace(family_terms, "over", "")) |> 
    mutate(across(everything(), \(x) 1 / sqrt(x))) 

  params <-
    cbind(fixed_params, family_params, random_params) |>
    as.data.frame() |>
    posterior::as_draws() |>
    posterior::summarise_draws(
      median,
      HDInterval::hdi
    ) |>
    dplyr::rename(term = variable, estimate = median, conf.low = lower, conf.high = upper) |>
    mutate(effect = ifelse(str_detect(term, "^SD "), "random", "fixed")) |>
    arrange(effect)
}
## ----end




## ---- make_brms_dharma_res_functions
make_brms_dharma_res <- function(brms_model, seed = 10, ...) {
                                        # equivalent to `simulateResiduals(lme4_model, use.u = FALSE)`
                                        # cores are set to 1 just to ensure reproducibility
    options(mc.cores = 1)
    on.exit(options(mc.cores = parallel::detectCores()))
    response <- brms::standata(brms_model)$Y
    ndraws <- nrow(as_draws_df(brms_model))
    manual_preds_brms <- matrix(0, ndraws, nrow(brms_model$data))
    random_terms <- insight::find_random(
                                 brms_model, split_nested = TRUE, flatten = TRUE
                             )
                                        # for this to have a similar output to `glmmTMB`'s default, we need to
                                        #   create new levels in the hierarchical variables, so then we can
                                        #   use `allow_new_levels = TRUE` and `sample_new_levels = "gaussian"` in
                                        #   `brms::posterior_epred`. This is equivalent to
                                        #   `simulateResiduals(lme4_model, use.u = FALSE)`. See details in
                                        #   `lme4:::simulate.merMod` and `glmmTMB:::simulate.glmmTMB`
                                        ## random_terms <- unlist(str_split(random_terms, ".\\+."))
  random_terms <- str_subset(random_terms, "\\+", negate = TRUE)
    new_data <- brms_model$data |>
        dplyr::mutate(across(
                   all_of(random_terms), \(x)paste0("NEW_", x) |> as.factor()
               ))
    set.seed(seed)
    brms_sims <- brms::posterior_predict(
                           brms_model, re_formula = NULL, newdata = new_data,
                           allow_new_levels = TRUE, sample_new_levels = "gaussian"
                       ) |>
        t()
    fitted_median_brms <- apply(brms_sims, 1, median)
    ## fitted_median_brms <- apply(
    ##     t(brms::posterior_epred(brms_model, ndraws = ndraws, re.form = NA)),
    ##     1,
    ##     mean)
    DHARMa::createDHARMa(
                simulatedResponse = brms_sims,
                observedResponse = response,
                fittedPredictedResponse = fitted_median_brms,
                ...
            )
}
## ----end


## ---- posterior_fitted.inla function
posterior_fitted.inla <- function(object, newdata = NULL, ndraws = 1000) {
  draws <- inla.posterior.sample(n=ndraws, result = object)
  contents <- object$misc$configs$contents

  form <- object$.args$formula
  gf <- INLA:::inla.interpret.formula(form)
  Xmat <- model.matrix(update(gf$fixf, NULL ~ .), newdata)
  nms <- colnames(Xmat)
  i_lp <- contents$start[contents$tag %in% nms]
  lp <- t(sapply(draws, function(x) x$latent[i_lp]))
  b <- tcrossprod(Xmat, lp)
  link <- object$misc$linkfunctions$names
  out <- my_ilink(b, link) 
  cellmeans.full <- newdata |> bind_cols(as.data.frame(out)) |>
    pivot_longer(cols=matches("^V[0-9]*$"),
      names_to='.draws',
      values_to='Values') |> 
    group_by(across(all_of(colnames(newdata))))
  cellmeans <-
          cellmeans.full |>
          posterior::as_draws() |>
          dplyr::select(-.draw) |>
          dplyr::mutate(.draw = as.integer(str_replace(.draws, "V", ""))) |> 
          dplyr::select(-.draws) 
  
}
## ----end


## ---- posterior_predict.inla function
posterior_predict.inla <- function(object, newdata = NULL, ndraws=250, include_re = TRUE, new_random_levels = FALSE) {
  draws <- inla.posterior.sample(n=ndraws, result = object)
  contents <- object$misc$configs$contents

  ## get the linear predictor
  if (is.null(newdata)) {
    mm <- object$model.matrix
    ## wch <- which(contents$tag=='Predictor')
    ## i_data <- contents$start[wch]:(1+contents$length[wch])
    ## mm <- mm[i_data,]
  } else {
      ## Fixed effects
      contrasts <- object$.args$contrasts
      form <- object$.args$formula
      gf <- INLA:::inla.interpret.formula(form)
      Xmat <- model.matrix(update(gf$fixf, NULL ~ .), newdata, contrasts = contrasts)

      ## Random effects
      if (!is.null(gf$randf)) {
          nms <- unlist(lapply(gf$random.spec, `[[`, "label"))
          ## newdata <- newdata %>% mutate(across(!!nms, ~ paste0('NEW',.x) %>% as.factor()))  
          Zmat <- lapply(nms, function(n) model.matrix(as.formula(paste0("~ 0 + ",n)), newdata))
          Zmat <- do.call('cbind', Zmat)
      }
  }
  ## Fixed effects
  nms <- colnames(Xmat)
  i_lp <- contents$start[contents$tag %in% nms]
  lp <- t(sapply(draws, function(x) x$latent[i_lp]))
  if (nrow(lp)==1) lp <- t(lp)
  b <- tcrossprod(Xmat, lp)
  link <- object$misc$linkfunctions$names
  out <- my_ilink(b, link) 

  ## family theta
  if ((ntheta <- object$misc$configs$ntheta) > 0) {
    theta <- vector('list', ntheta)
    for (t in 1:ntheta) {
      theta[[t]] <- sapply(draws, function(x) 1/sqrt(x$hyperpar[[t]]))
    }
  }

  ## Random effects
  if (gf$n.random>0 & include_re) {
      nms <- unlist(lapply(gf$random.spec, `[[`, "label"))
      i_lp <- contents$start[contents$tag %in% nms]
      lp <- ii_lp <- vector('list', length(nms))
      for (i in 1:length(nms)) {
          ii_lp[[i]] <- i_lp[[i]]:(i_lp[[i]] + (contents$length[contents$tag %in% nms] - 1)[i])
          lp[[i]] <- t(sapply(draws, function(x) x$latent[ii_lp[[i]]]))
          if (new_random_levels) lp[[i]] <- apply(lp[[i]], 2,
                                                  function(x) rnorm(n = length(x), mean = 0, sd = theta[[1+i]]))
      }
      lp <- do.call(`cbind`, lp)
      ## i_lp <- unlist(lapply(1:length(i_lp),
      ##                       function(i)
      ##                           i_lp[i]:(i_lp[i] + (contents$length[contents$tag %in% nms] - 1)[i])))
      
      ## i_lp <- i_lp:(i_lp + contents$length[contents$tag %in% nms] - 1)
      ## lp <- t(sapply(draws, function(x) x$latent[i_lp]))
      ## if (new_random_levels) lp <- apply(lp, 2, function(x) rnorm(n = length(x), mean = 0, sd = theta[[2]]))
      r <- tcrossprod(Zmat, lp)
      link <- object$misc$linkfunctions$names
      out <- my_ilink(b+r, link) 
  }

  fam <- object$.args$family
  N <- nrow(out)
  if (fam == 'gaussian') {
    rdist <- function(N, x, sigma) rnorm(N, x, sigma)
    out <- apply(out, 2, rdist, N=N, sigma=theta[[1]])
  }
  if (fam=='poisson') {
    rdist <- function(N, x) rpois(N, x)
    out <- apply(out, 2, rdist, N=N)
  }
  if (fam=="nbinomial") {
      wch <- grep(
          "size for the nbinomial observations",
          names(draws[[1]]$hyperpar)
      )
      size <- sapply(draws, function(x) x$hyperpar[[wch]])
      ## size <- inla.hyperpar.sample(n = ndraws, result = mod.inla)[, wch]
      ## rdist <- function(N, x, size) MASS::rnegbin(N, mu = x, size) #rnbinom(N, mu = x, prob)
      rdist <- function(N, x, size) rnbinom(N, mu = x, size) #rnbinom(N, mu = x, prob)
      out <- apply(out, 2, rdist, N=N, size=size)
  }
  if (fam=="binomial") {
    ## size <- ifelse(is.null(object$.args$Ntrials), 1, object$.args$Ntrials)
    if (is.null(object$.args$Ntrials)) {
        size <- rep(1, N)
    } else if (length(object$.args$Ntrials) == 1) {
        size <- rep(object$.args$Ntrials, N)
    } else {
        size <- object$.args$Ntrials[1:N]
    }
    rdist <- function(N, x, size) rbinom(N, prob = x, size)
    out <- apply(out, 2, rdist, N=N, size=size)
  }
  if (fam=="beta") {
      wch <- grep(
          "precision parameter for the beta observations",
          names(draws[[1]]$hyperpar)
      )
      phi <- sapply(draws, function(x) x$hyperpar[[wch]])
      rdist <- function(N, x, phi) rbeta(N, shape1 = x*phi, shape2 = (1 - x)*phi)
      ## out <- apply(out, 2, rdist, N=N, phi=phi)
      out <- sapply(seq_len(ncol(out)), function(i)
          rdist(N = N, x = out[,i], phi = phi[i])) %>%
          matrix(ncol = ncol(out), byrow = TRUE)
  }
   if (fam=="betabinomial") {
      wch <- grep(
          "overdispersion for the betabinomial observations",
          names(draws[[1]]$hyperpar)
      )
      phi <-1/sapply(draws, function(x) x$hyperpar[[wch]])
      if (is.null(object$.args$Ntrials)) {
          size <- rep(1, N)
      } else {
          size <- object$.args$Ntrials[1:N]
      }
      ## size <- ifelse(is.null(object$.args$Ntrials), 1, object$.args$Ntrials)
      rdist <- function(N, x, phi, size) {
          rbinom(N,
                 prob = rbeta(N, shape1 = x*phi, shape2 = (1 - x)*phi),
                 size = size)
          }
      out <- sapply(seq_len(ncol(out)), function(i)
          rdist(N = N, x = out[,i], phi = phi[i], size = size)) %>%
          matrix(ncol = ncol(out), byrow = TRUE)
  }
  if (fam=="zeroinflatedpoisson1") {
      wch <- grep(
          "zero-probability parameter for zero-inflated poisson_1",
          names(draws[[1]]$hyperpar)
      )
      phi <- sapply(draws, function(x) x$hyperpar[[wch]])
      rdist <- function(N, x, phi) rbinom(N, size = 1, prob = 1 - phi) * rpois(N, lambda = x)
      out <- apply(out, 2, rdist, N=N, phi=phi)
  }
  
  
  t(out)
}



## ----end

## ---- my_ilink function
my_ilink <- function(x, link) {
  switch(link, identity = x, log = exp(x), logm1 = expp1(x), 
         log1p = expm1(x), inverse = 1/x, sqrt = x^2, `1/mu^2` = 1/sqrt(x), 
         tan_half = 2 * atan(x), logit = inv_logit(x), probit = pnorm(x), 
         cauchit = pcauchy(x), cloglog = inv_cloglog(x), probit_approx = pnorm(x), 
        softplus = log1p_exp(x), stop2("Link '", link, "' not supported."))
  }
## ----end

## posterior_predict.inla <- function(object, newdata = NULL, ndraws=250, include_re = TRUE, new_random_levels = FALSE) {
##   draws <- inla.posterior.sample(n=ndraws, result = object)
##   contents <- object$misc$configs$contents

##   ## get the linear predictor
##   if (is.null(newdata)) {
##     mm <- object$model.matrix
##     ## wch <- which(contents$tag=='Predictor')
##     ## i_data <- contents$start[wch]:(1+contents$length[wch])
##     ## mm <- mm[i_data,]
##   } else {
##       ## Fixed effects
##       contrasts <- object$.args$contrasts
##       form <- object$.args$formula
##       gf <- INLA:::inla.interpret.formula(form)
##       Xmat <- model.matrix(update(gf$fixf, NULL ~ .), newdata, contrasts = contrasts)

##       ## Random effects
##       if (!is.null(gf$randf)) {
##           nms <- unlist(lapply(gf$random.spec, `[[`, "label"))
##           ## newdata <- newdata %>% mutate(across(!!nms, ~ paste0('NEW',.x) %>% as.factor()))  
##           Zmat <- lapply(nms, function(n) model.matrix(as.formula(paste0("~ 0 + ",n)), newdata))
##           Zmat <- do.call('cbind', Zmat)
##       }
##   }
##   ## Fixed effects
##   nms <- colnames(Xmat)
##   i_lp <- contents$start[contents$tag %in% nms]
##   lp <- t(sapply(draws, function(x) x$latent[i_lp]))
##   if (nrow(lp)==1) lp <- t(lp)
##   b <- tcrossprod(Xmat, lp)
##   link <- object$misc$linkfunctions$names
##   out <- my_ilink(b, link) 

##   ## family theta
##   if ((ntheta <- object$misc$configs$ntheta) > 0) {
##     theta <- vector('list', ntheta)
##     for (t in 1:ntheta) {
##       theta[[t]] <- sapply(draws, function(x) 1/sqrt(x$hyperpar[[t]]))
##     }
##   }

##   ## Random effects
##   if (gf$n.random>0 & include_re) {
##       nms <- unlist(lapply(gf$random.spec, `[[`, "label"))
##       i_lp <- contents$start[contents$tag %in% nms]
##       lp <- ii_lp <- vector('list', length(nms))
##       for (i in 1:length(nms)) {
##           ii_lp[[i]] <- i_lp[[i]]:(i_lp[[i]] + (contents$length[contents$tag %in% nms] - 1)[i])
##           lp[[i]] <- t(sapply(draws, function(x) x$latent[ii_lp[[i]]]))
##           if (new_random_levels) lp[[i]] <- apply(lp[[i]], 2,
##                                                   function(x) rnorm(n = length(x), mean = 0, sd = theta[[1+i]]))
##       }
##       lp <- do.call(`cbind`, lp)
##       ## i_lp <- unlist(lapply(1:length(i_lp),
##       ##                       function(i)
##       ##                           i_lp[i]:(i_lp[i] + (contents$length[contents$tag %in% nms] - 1)[i])))
      
##       ## i_lp <- i_lp:(i_lp + contents$length[contents$tag %in% nms] - 1)
##       ## lp <- t(sapply(draws, function(x) x$latent[i_lp]))
##       ## if (new_random_levels) lp <- apply(lp, 2, function(x) rnorm(n = length(x), mean = 0, sd = theta[[2]]))
##       r <- tcrossprod(Zmat, lp)
##       link <- object$misc$linkfunctions$names
##       out <- my_ilink(b+r, link) 
##   }

##   fam <- object$.args$family
##   N <- nrow(out)
##   if (fam == 'gaussian') {
##     rdist <- function(N, x, sigma) rnorm(N, x, sigma)
##     out <- apply(out, 2, rdist, N=N, sigma=theta[[1]])
##   }
##   if (fam=='poisson') {
##     rdist <- function(N, x) rpois(N, x)
##     out <- apply(out, 2, rdist, N=N)
##   }
##   if (fam=="nbinomial") {
##       wch <- grep(
##           "size for the nbinomial observations",
##           names(draws[[1]]$hyperpar)
##       )
##       size <- sapply(draws, function(x) x$hyperpar[[wch]])
##       ## size <- inla.hyperpar.sample(n = ndraws, result = mod.inla)[, wch]
##       ## rdist <- function(N, x, size) MASS::rnegbin(N, mu = x, size) #rnbinom(N, mu = x, prob)
##       rdist <- function(N, x, size) rnbinom(N, mu = x, size) #rnbinom(N, mu = x, prob)
##       out <- apply(out, 2, rdist, N=N, size=size)
##   }
##   if (fam=="binomial") {
##     ## size <- ifelse(is.null(object$.args$Ntrials), 1, object$.args$Ntrials)
##     if (is.null(object$.args$Ntrials)) {
##         size <- rep(1, N)
##     } else if (length(object$.args$Ntrials) == 1) {
##         size <- rep(object$.args$Ntrials, N)
##     } else {
##         size <- object$.args$Ntrials[1:N]
##     }
##     rdist <- function(N, x, size) rbinom(N, prob = x, size)
##     out <- apply(out, 2, rdist, N=N, size=size)
##   }
##   if (fam=="beta") {
##       wch <- grep(
##           "precision parameter for the beta observations",
##           names(draws[[1]]$hyperpar)
##       )
##       phi <- sapply(draws, function(x) x$hyperpar[[wch]])
##       rdist <- function(N, x, phi) rbeta(N, shape1 = x*phi, shape2 = (1 - x)*phi)
##       ## out <- apply(out, 2, rdist, N=N, phi=phi)
##       out <- sapply(seq_len(ncol(out)), function(i)
##           rdist(N = N, x = out[,i], phi = phi[i])) %>%
##           matrix(ncol = ncol(out), byrow = TRUE)
##   }
##    if (fam=="betabinomial") {
##       wch <- grep(
##           "overdispersion for the betabinomial observations",
##           names(draws[[1]]$hyperpar)
##       )
##       phi <-1/sapply(draws, function(x) x$hyperpar[[wch]])
##       if (is.null(object$.args$Ntrials)) {
##           size <- rep(1, N)
##       } else {
##           size <- object$.args$Ntrials[1:N]
##       }
##       ## size <- ifelse(is.null(object$.args$Ntrials), 1, object$.args$Ntrials)
##       rdist <- function(N, x, phi, size) {
##           rbinom(N,
##                  prob = rbeta(N, shape1 = x*phi, shape2 = (1 - x)*phi),
##                  size = size)
##           }
##       out <- sapply(seq_len(ncol(out)), function(i)
##           rdist(N = N, x = out[,i], phi = phi[i], size = size)) %>%
##           matrix(ncol = ncol(out), byrow = TRUE)
##   }
##   if (fam=="zeroinflatedpoisson1") {
##       wch <- grep(
##           "zero-probability parameter for zero-inflated poisson_1",
##           names(draws[[1]]$hyperpar)
##       )
##       phi <- sapply(draws, function(x) x$hyperpar[[wch]])
##       rdist <- function(N, x, phi) rbinom(N, size = 1, prob = 1 - phi) * rpois(N, lambda = x)
##       out <- apply(out, 2, rdist, N=N, phi=phi)
##   }
  
  
##   t(out)
## }

## my_ilink <- function(x, link) {
##   switch(link, identity = x, log = exp(x), logm1 = expp1(x), 
##          log1p = expm1(x), inverse = 1/x, sqrt = x^2, `1/mu^2` = 1/sqrt(x), 
##          tan_half = 2 * atan(x), logit = inv_logit(x), probit = pnorm(x), 
##          cauchit = pcauchy(x), cloglog = inv_cloglog(x), probit_approx = pnorm(x), 
##         softplus = log1p_exp(x), stop2("Link '", link, "' not supported."))
##   }

## ---- pit_hist_plot function
pit_hist_plot <- function(mod, i.mod) {
    pit <- mod$cpo$pit[i.mod]
    g <- ggplot(data = NULL, aes(x = pit)) +
        ## geom_histogram() +
        geom_density() +
        scale_y_continuous("Density") +
        scale_x_continuous("PIT") +
        theme_bw()
    return(g)
}
## ----end
## ---- pit_resid_plot function
pit_resid_plot <- function(mod, i.mod) {
    pit <- mod$cpo$pit[i.mod]
    g <- ggplot(data = NULL, aes(y = pit, x = mod$summary.fitted.values$mean[i.mod])) +
        geom_point() +
        scale_y_continuous("PIT values") +
        scale_x_continuous("Posterior fitted values") +
        theme_bw()
    return(g)
}
## ----end
## ---- pit_plot function
pit_plot <- function(mod, i.mod) {
    pit <- na.omit(mod$cpo$pit[i.mod])
    data <- na.omit(mod$.args$data[i.mod,])
    n <- nrow(data)
    g <- ggplot(data = NULL, aes(y = pit, x = 1:n)) +
        geom_point() +
        scale_y_continuous("PIT values") +
        scale_x_continuous("1:n") +
        theme_bw()
    return(g)
}
## ----end
## ---- pit_qq_plot function
pit_qq_plot <- function(mod, i.mod, logit_scale = TRUE) {
    logit <- function(x) log(x/(1-x))
    pit <- na.omit(mod$cpo$pit[i.mod])
    data <- na.omit(mod$.args$data[i.mod,])
    n <- nrow(data)
    uniquant <- (1:n)/(n+1)
    g <- ggplot(data = NULL, aes(y = sort(pit), x = uniquant))
    if (logit_scale) 
        g <- ggplot(data = NULL, aes(y = logit(sort(pit)), x = logit(uniquant)))
    g <- g + geom_abline(slope = 1) +
        geom_point() +
        scale_y_continuous("Sorted PIT values") +
        scale_x_continuous("Uniform quantiles") +
        theme_bw()
    return(g)
}
## ----end
## ---- pp_check.inla function
pp_check.foo <- function(object, type = c("multiple", "overlaid"), ...) {
  type <- match.arg(type)
  y <- object[["y"]]
  yrep <- object[["yrep"]]

  if (type == "overlaid") {
    ppc_dens_overlay(y, yrep, ...) 
  } else {
    ppc_hist(y, yrep, ...)
  }
}
## ----end


## Take a data.frame that has values for each Dist and calculate the
## fold difference between Before and each of the others. Example of
## the data frame follows

## # A tibble: 12 × 8
##     Dist.time Dist   Values
##     <fct>     <chr>     <dbl>
##   1 Before    Before 0.00730
##   3 After     b      0.000938
##   5 After     c      0.00209
##   7 After     d      0.00283
##   9 After     s      0.00271
##  11 After     u      0.00607


## ---- before_vs_afters_function
before_vs_afters <- function(.x) {
  .x <- .x |>
    ## mutate(Dist = factor(Dist,
    ##   levels = c("Before", "s", "c", "d", "b", "u")
    ## )) |>
    mutate(Dist = forcats::fct_relevel(Dist, "Before")) |>
    arrange(Dist)
  N <- nrow(.x)
  xmat <- cbind(-1, 1 * contr.treatment(N, base = 1, contrast = TRUE))
  xmat <- xmat[-1, ]
  x <- log(as.vector(as.vector(.x$Values)))
  ## data.frame(
  ##   Dist = .x$Dist[-1],
  ##   Values = exp(as.vector(x %*% t(xmat)))
  ## )
  if (N>2) {  #if there are more than 2 levels
    data.frame(
      Dist = .x$Dist[-1],
      Values = exp(as.vector(x %*% t(xmat)))
    )
  } else {  ## if there is only one contrast
    data.frame(
      Dist = .x$Dist[-1],
      Values = exp(as.vector(x %*% (xmat)))
    )
  }
}
## ----end

## ---- brm_generate_newdata_no_dist
brm_generate_newdata_no_dist <- function(data, GBR = FALSE) {
  newdata <-
    data |>
    mutate(zone_depth = NA) |> 
    dplyr::select(SecShelfYr, zone_depth) |>
    ## filter(SecShelf == "CA M") |> 
    distinct() |> 
    mutate(AIMS_REEF_NAME = NA, Site = NA, Transect = NA, event = NA, total.points = 1)
  if (GBR) newdata <- newdata |> mutate(SecShelf = NA)
  return(newdata)
}
## ----end

## ---- brm_generate_newdata
brm_generate_newdata <- function(data, GBR = FALSE) {
  if (GBR) {
    newdata <-
      data |>
      dplyr::select(Dist, s, c, d, b, u) |>
      ## filter(SecShelf == "CA M") |> 
      distinct() |>
      rowwise() |>
      mutate(Dist2 = paste(c("s", "c", "d", "b", "u")[which(c(s, c, d, b, u) == 1)], collapse = "")) |>
      ungroup()
  } else {
    newdata <-
      data |>
      mutate(zone_depth = NA) |> 
      dplyr::select(SecShelf, zone_depth, Dist, s, c, d, b, u) |>
      ## filter(SecShelf == "CA M") |> 
      distinct() |>
      rowwise() |>
      mutate(Dist2 = paste(c("s", "c", "d", "b", "u")[which(c(s, c, d, b, u) == 1)], collapse = "")) |>
      ungroup()
  }

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
    mutate(AIMS_REEF_NAME = NA, Site = NA, Transect = NA, event = NA, total.points = 1) |>
    filter(!is.na(s)) |>
    mutate(Dist3 = ifelse(cummulative, Dist3, Dist2))
  if (GBR) newdata <- newdata |> mutate(SecShelf = NA)
  return(newdata)
}
## ----end

## ---- brm_generate_cellmeans_no_dist_function
brm_generate_cellmeans_no_dist <- function(newdata, pred) {
  newdata |>
    cbind(pred) |> 
    pivot_longer(cols = matches("^[0-9]*$"), names_to = ".draws") |>
    mutate(Pred = (value)) |>
    separate(SecShelfYr, into = c("A_SECTOR", "SHELF", "REPORT_YEAR"), sep = " ") 
}
## ----end

## ---- brm_generate_cellmeans_function
brm_generate_cellmeans <- function(newdata, pred) {
  newdata |>
    cbind(pred) |> 
    pivot_longer(cols = matches("^[0-9]*$"), names_to = ".draws") |>
    mutate(Pred = (value)) |>
    mutate(
      Dist = ifelse(Dist2 == "" | is.na(Dist2), "Before", Dist2),
      Dist = forcats::fct_relevel(Dist, "Before"),
      Dist3 = forcats::fct_relevel(Dist3, "Before")
    ) |> 
    separate(SecShelf, into = c("A_SECTOR", "SHELF"), sep = " ") 
}
## ----end

## ---- disturbance_palette
dist_palette <- tribble(
  ~lab, ~color,
  "Before", "#000000",
  "b",    "#e41a1c",
  "s",    "#377eb8",
  "sb",   "#8E4C6A",
  "d",    "#c6c6c6",
  "c",    "#4daf4a",
  "cb",   "#996533",
  "sc",   "#429781",
  "scb",  "#786D5F",
  "cd",   "#275825",
  "cu",   "#a69725",
  "u",    "#ff7f00",
  "su",   "pink",
  "sd",   "purple",
  )
## ----end

## ---- brm_partial_plot_no_dist_function
brm_partial_plot_no_dist <- function(cellmeans_summ_brm) {
  p <-
    cellmeans_summ_brm |>
    mutate(
      A_SECTOR = factor(A_SECTOR,
        levels = c("CG", "PC", "CL", "CA", "IN", "TO", "WH", "PO", "SW", "CB"),
        labels = c(
          "Cape Grenville", "Princess Charlotte Bay", "Cooktown Lizard", "Cairns",
          "Innisfail", "Townsville", "Whitsunday", "Pompey", "Swain", "Capricon Bunker"
        )
      ),
      SHELF = factor(SHELF,
        levels = c("I", "M", "O"),
        labels = c("Inshore", "Midshelf", "Offshore")
      )
    ) |>
    ggplot(aes(y = median, x = as.numeric(as.character(REPORT_YEAR)))) +
    geom_ribbon(aes(ymin = lower, ymax = upper), fill = "orange", alpha = 0.3) +
    geom_line() +
    geom_pointrange(aes(ymin = lower, ymax = upper),
      position = position_dodge(width = 0.5),
      show.legend = FALSE) +
    ## scale_y_continuous("Seriatopora cover (%)", trans = scales::log_trans(), labels =  function(x) x*100) +
    theme_bw() +
    theme(
      axis.title.y = element_text(size = rel(1.25),
        margin = margin(r = 1, unit = "char")),
      axis.title.x = element_text(size = rel(1.25),
        margin = margin(t = 1, unit = "char")),
      axis.text.y = element_text(size = rel(1.0)),
      axis.text.x = element_text(angle = 0, hjust = 0.5),
      strip.background = element_rect(fill = "#95b8d1"),
      strip.text = element_text(size = rel(1.25))) +
    facet_grid(A_SECTOR ~ SHELF,
      scales = "free",
      labeller = label_wrap_gen(width = 12, multi_line = TRUE)
    )
  p
}
## ----end
## ---- brm_partial_plot_function
brm_partial_plot <- function(cellmeans_summ_brm) {
  before_underlay_band <- cellmeans_summ_brm |>
    ungroup() |>
    filter(Dist2 == "Before") |>
    dplyr::select(A_SECTOR, SHELF, Dist2, lower, upper) |>
    distinct() |>
    mutate(
      Dist2 = factor(Dist2,
        levels = c("Before", "b", "c", "s", "d", "u"),
        labels = c("Before", "Bleaching", "COTS", "Storms", "Disease", "Unknown")
      ),
      ## A_SECTOR = factor(A_SECTOR,
      ##   levels = c("CG", "CL", "CA", "IN", "TO", "WH", "PO", "SW", "CB"),
      ##   labels = c(
      ##     "Cape Grenville", "Cooktown Lizard", "Cairns",
      ##     "Innisfail", "Townsville", "Whitsunday", "Pompey", "Swain", "Capricon Bunker"
      ##   )
      ## ),
      A_SECTOR = factor(A_SECTOR,
        levels = c("CG", "PC", "CL", "CA", "IN", "TO", "WH", "PO", "SW", "CB"),
        labels = c(
          "Cape Grenville", "Princess Charlotte Bay", "Cooktown Lizard", "Cairns",
          "Innisfail", "Townsville", "Whitsunday", "Pompey", "Swain", "Capricon Bunker"
        )
      ),
      SHELF = factor(SHELF,
        levels = c("I", "M", "O"),
        labels = c("Inshore", "Midshelf", "Offshore")
      )
    )

  p <-
    cellmeans_summ_brm |>
    mutate(Dist2 = factor(Dist2,
      levels = c("Before", "b", "c", "s", "d", "u"),
      labels = c("Before", "Bleaching", "COTS", "Storms", "Disease", "Unknown"))) |> 
    ## mutate(A_SECTOR = factor(A_SECTOR,
    ##   levels = c("CG", "CL", "CA", "IN", "TO", "WH", "PO", "SW", "CB"),
    ##   labels = c(
    ##     "Cape Grenville", "Cooktown Lizard", "Cairns",
    ##     "Innisfail", "Townsville", "Whitsunday", "Pompey", "Swain", "Capricon Bunker"
    ##   )
    ## ),
    mutate(A_SECTOR = factor(A_SECTOR,
      levels = c("CG", "PC", "CL", "CA", "IN", "TO", "WH", "PO", "SW", "CB"),
      labels = c(
        "Cape Grenville", "Princess Charlotte Bay", "Cooktown Lizard", "Cairns",
        "Innisfail", "Townsville", "Whitsunday", "Pompey", "Swain", "Capricon Bunker"
      )
    ),
    SHELF =  factor(SHELF, levels = c("I", "M", "O"),
      labels = c("Inshore", "Midshelf", "Offshore")))|> 
    mutate(Dist.time = ifelse(Dist == "Before", "Before", "After")) |>
    mutate(Dist.time = factor(Dist.time, levels = c("Before", "After"))) |>
    mutate(Dist = forcats::fct_relevel(Dist, "Before")) |>
    ggplot(aes(y = median, x = Dist2, colour = Dist3, shape = cummulative)) +
    geom_rect(
      data = before_underlay_band,
      aes(ymin = lower, xmin = -Inf, ymax = upper, xmax = Inf),
      fill = "#d6d6d6",
      inherit.aes = FALSE
    ) +
    geom_pointrange(aes(ymin = lower, ymax = upper),
      position = position_dodge(width = 0.5),
      show.legend = FALSE) +
    scale_shape_manual(values = c(16, 21)) +
    scale_y_continuous("Seriatopora cover (%)", labels =  function(x) x*100) +
    scale_x_discrete("Disturbance type",
      breaks = c("Before", "Bleaching", "COTS", "Storms", "Disease", "Unknown"),
      expand = c(0, 0.3)) +
    scale_colour_manual("Disurbances",
      breaks = dist_palette$lab,
      values = dist_palette$color) +
    theme_bw() +
    theme(
      axis.title.y = element_text(size = rel(1.25),
        margin = margin(r = 1, unit = "char")),
      axis.title.x = element_text(size = rel(1.25),
        margin = margin(t = 1, unit = "char")),
      axis.text.y = element_text(size = rel(1.0)),
      axis.text.x = element_text(angle = 0, hjust = 0.5),
      strip.background = element_rect(fill = "#95b8d1"),
      strip.text = element_text(size = rel(1.25))) +
    facet_grid(A_SECTOR ~ SHELF,
      scales = "free",
      labeller = label_wrap_gen(width = 12, multi_line = TRUE)
    )
  p
}
## ----end

## ---- brm_calc_effect_function
brm_calc_effect <- function(cellmeans_brm, GBR = FALSE) {
  if (GBR) {
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
  } else {
    cellmeans_brm |>
      mutate(Values = Pred,
        Dist =  Dist3) |>
      nest(.by = c(A_SECTOR, SHELF, .draws)) |>
      mutate(eff = map(
        .x = data,
        .f = ~ before_vs_afters(.x)
      )) |>
      dplyr::select(-data) |>
      unnest(c(eff)) |>
      mutate(Dist2 = Dist) |>
      mutate(number_of_dist = str_length(Dist2)) |> 
      mutate(cummulative = ifelse(str_length(Dist2) > 1, TRUE, FALSE)) |> 
      separate_longer_position(Dist2, width = 1, keep_empty = TRUE) 
  }
}
## ----end


## ---- brm_calc_effect_function
brm_calc_effect_hier <- function(cellmeans_brm) {
  cellmeans_brm |>
    mutate(Values = value) |>
    nest(.by = c(A_SECTOR, SHELF, .draws)) |>
    mutate(eff = map(
      .x = data,
      .f = ~ before_vs_afters(.x)
    )) |>
    dplyr::select(-data) |>
    unnest(c(eff)) 
}
## ----end

## ---- brm_effects_plot
brm_effects_plot <- function(eff_brm) {
  dist_palette <- dist_palette |>
    filter(lab != "Before") |>
    droplevels()

  dist_order <- tribble(
    ~levels, ~labels,
    "b",     "Bleaching",
    "c",     "COTS",
    "s",     "Storms",
    "u",     "Unknown",
    "d",     "Disease"
  )
  dist_order <- dist_order[5:1,]

  dist_full_order <- tribble(
    ~levels, ~labels,
    "b",     "Bleaching",
    "c",     "COTS",
    "cb",    "COTS/Bleaching",
    "s",     "Storms",
    "u",     "Unknown",
    "d",     "Disease",
    "sb",    "Storms/Bleaching",
    "sc",    "Storms/COTS",
    "scb",   "Storms/Bleaching/COTS",
    "sd",    "Storms/Disease",
    "su",    "Storms/Unknown",
    "db",    "Disease/Bleaching",
  )

  dist_palette <- dist_palette |>
    full_join(dist_full_order, by = c("lab" = "levels"))

  eff_brm |>
    mutate(Dist2 = forcats::fct_relevel(Dist2, dist_order$levels)) |>
    ## mutate(Dist = forcats::fct_relevel(Dist, "b", "c", "s", "u", "d")) |>
    mutate(Dist = factor(Dist, levels = dist_full_order$levels, labels = dist_full_order$labels)) |>
    ## mutate(Dist = forcats::fct_relevel(Dist, dist_order$levels)) |>
    ## mutate(A_SECTOR = factor(A_SECTOR,
    ##   levels = c("CG", "CL", "CA", "IN", "TO", "WH", "PO", "SW", "CB"),
    ##   labels = c(
    ##     "Cape Grenville", "Cooktown Lizard", "Cairns",
    ##     "Innisfail", "Townsville", "Whitsunday", "Pompey", "Swain", "Capricon Bunker"
    ##   )
    ## ),
    mutate(A_SECTOR = factor(A_SECTOR,
      levels = c("CG", "PC", "CL", "CA", "IN", "TO", "WH", "PO", "SW", "CB"),
      labels = c(
        "Cape Grenville", "Princess Charlotte Bay", "Cooktown Lizard", "Cairns",
        "Innisfail", "Townsville", "Whitsunday", "Pompey", "Swain", "Capricon Bunker"
      )
    ),
    SHELF =  factor(SHELF, levels = c("I", "M", "O"), labels = c("Inshore", "Midshelf", "Offshore")))|> 
    ## ggplot(aes(x = Values, y = Dist2, group = Dist, shape = cummulative)) +
    ggplot(aes(x = Values, y = Dist2, group = Dist, shape = factor(number_of_dist))) +
    geom_vline(xintercept = 1, linetype = "dashed") +
    ## stat_slabinterval(aes(color = Dist, point_size = number_of_dist),
    ## stat_summary(geom = "point", aes(color = Dist, size = after_stat(number_of_dist)), fun = mean,
    ##   position = position_dodge(width = 0.5, preserve = "total")) +
    ## geom_point(data = eff_summ_brm_summary,
    ##   aes(x = Mean, color = Dist, shape = cummulative, size = number_of_dist),
    ##   position = position_dodge(width = 0.5, preserve = "total")) +
    stat_slabinterval(aes(color = Dist, point_size = number_of_dist),
      position = position_dodge(width = 0.5, preserve = "total"),
      height = 0,
      fill = "white", 
      ## point_size = 0,
      show.legend = c(size = FALSE)) +
    scale_point_size_continuous(range = c(2, 4), guide = NULL) +
      ## scale_size_binned(breaks = c(1, 2, 3)) +
      ## scale_size_binned(breaks = c(1, 3), limits = c(1,3), range = c(0.5, 10)) +
      ## scale_shape_manual(breaks = c(FALSE, TRUE), values = c(16, 21)) +
    scale_shape_manual(breaks = c(1, 2, 3), values = c(16, 21, 23), guide = NULL) +
    scale_x_continuous("Effect size (% change in Seriatopora cover after disturbance)",
      trans = scales::log2_trans(),
      labels = function(x) round(100*(x-1), 2),
      expand = c(0, 0),
      ## breaks = c(1 / 200, 1 / 50, 1 / 20, 1 / 10, 1 / 5, 0.3, 1 / 2, 1, 2, 3, 6),
      breaks = c(1 / 200, 1 / 50, 1/20, 1 / 10, 1 / 4, 1 / 2, 1, 2, 4, 10),
      #limits = c(0.001, 10)
      ) +
    scale_colour_manual("Disurbances",
      guide = guide_legend(),
      breaks = dist_palette$labels,
      values = dist_palette$color) +
    scale_y_discrete("Disturbance type",
      ## breaks = c("b", "c", "s", "u", "d"),
      ## labels = c("Bleaching", "COTS", "Storms", "Unknown", "Disease"),
      breaks = dist_order$levels,
      labels = dist_order$labels) +
      ## expand = c(0, 0.3)) +
    theme_bw() +
    theme(axis.title.y = element_text(size = rel(1.25), margin = margin(r = 1, unit = "char")),
      axis.title.x = element_text(size = rel(1.25), margin = margin(t = 1, unit = "char")),
      axis.text.y = element_text(size = rel(1.0)),
      strip.background = element_rect(fill = "#95b8d1"),
      strip.text = element_text(size = rel(1.25))) +
    facet_grid(A_SECTOR ~ SHELF,
      scales = "free",
      labeller = label_wrap_gen(width = 12, multi_line = TRUE)
    ) +
    coord_cartesian(clip = "on", xlim = c(0.001, 11)) #+
    #guides(.width = FALSE, colour = )

}
## ----end

## ---- brm_effects_plot_gbr
brm_effects_plot_gbr <- function(eff_brm, height = 1.5, p = 0.1, dist_order = NULL, label_offset = 0.4) {

  if (is.null(dist_order)) {
    dist_order <- tribble(
      ~levels, ~labels,
      "b",     "Bleaching",
      "c",     "COTS",
      "s",     "Storms",
      "u",     "Unknown",
      "d",     "Disease"
    )
  }
  ## eff_brm <- eff_brm |>
  ##   mutate(Dist = factor(Dist, levels = rev(levels(Dist))))

  dist_order <- dist_order[5:1,]
  my_density <- function(x, p) {
    D <- density(log(x))
    exp(as.vector(D$x[D$y > p][1]))
  }

  eff_brm_lab <- eff_brm |>
    dplyr::select(-.draws) |>
    group_by(Dist) |>
    summarise_draws(
      median,
      HDInterval::hdi,
      D1 = ~ my_density(., p),
      D = ~ quantile(., prob = 0.05, names = FALSE),
      P = ~ mean(. < 1)
    )

  small_dist_palette <- dist_palette |>
    filter(lab %in% c("b", "c", "d", "s", "u")) |>
    droplevels() |>
    rename(Dist = lab)

  p <-
    eff_brm |>
    left_join(small_dist_palette) |>
    ## mutate(Dist = forcats::fct_relevel(Dist, "b", "c", "s", "u", "d")) |>
    mutate(Dist = forcats::fct_relevel(Dist, dist_order$levels)) |>
    ggplot(aes(x = Values, y = Dist)) +
    geom_vline(xintercept = 1, linetype = "dashed") +
    stat_slabinterval(aes(color = Dist),
      position = position_dodge(width = 0.3, preserve = "single"),
      show.legend = FALSE) +
    ## stat_slab(height = height/2) +
    stat_slab(aes(fill = Dist, slab_alpha = after_stat(x) < 1),
      fill_type = "segments",
      height = height,
      expand = FALSE, trim = TRUE, density = "bounded",
      width = 0.95,
      colour = "grey10", alpha = 0.4, show.legend = FALSE) +
    geom_text(data = eff_brm_lab |>
                ## mutate(Dist = factor(Dist, levels = c("b", "c", "s", "u", "d"))),
                mutate(Dist = factor(Dist, levels = dist_order$levels)),
      aes(y = as.numeric(as.factor(Dist))+0.2,
        x = D1, #1.5,
        label = sprintf("P[decline]==%0.3f", P)
      ),
      nudge_y =  0.2, hjust = 1, parse = TRUE) +
    geom_label(
      data = eff_brm_lab |>
        ## mutate(Dist = factor(Dist, levels = c("b", "c", "s", "u", "d"))),
        mutate(Dist = factor(Dist, levels = dist_order$levels)),
      aes(
        y = as.numeric(as.factor(Dist)) + label_offset,
        x = median,
        label = sprintf("%0.2f%%", 100 * (median - 1))
      )
    ) +
    scale_fill_viridis_d() +
    scale_colour_viridis_d() +
    scale_slab_alpha_discrete(range = c(0.2, 0.7))+
    scale_x_continuous("Effect size (% change in Seriatopora cover after disturbance)",
      trans = scales::log2_trans(),
      labels = function(x) round(100*(x-1), 2),
      expand = c(0, 0),
      breaks = c(1 / 200, 1 / 50, 1/20, 1 / 10, 1 / 4, 1 / 2, 1, 2, 4, 10),
      ) +
    scale_y_discrete("Disturbance type",
      ## breaks = c("b", "c", "s", "u", "d"),
      ## labels = c("Bleaching", "COTS", "Storms", "Unknown", "Disease"),
      breaks = dist_order$levels,
      labels = dist_order$labels,
      expand = c(0, 0.3)
    ) + 
    ## limits = rev(dist_order$levels)) +
    theme_bw() +
    theme(axis.title.y = element_text(size = rel(1.25), margin = margin(r = 1, unit = "char")),
      axis.title.x = element_text(size = rel(1.25), margin = margin(t = 1, unit = "char")),
      axis.text.y = element_text(size = rel(1.25))) +
    coord_cartesian(clip = "on", xlim = c(0.001, 11))
  p
}
## ----end
