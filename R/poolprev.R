# Save this file as `R/poolprev.R`

#' Estimate prevalence from pooled test results
#'
#' `poolprev` is used to estimate prevalence from pooled test results over time.
#' It fits a Bayesian model using the Stan framework.
#'
#'
#'**Models**
#'
#' Two models are implemented, which differ in the way the correlation over time is handled.
#'
#' In the default model (selected by setting `method="GP"`), a hidden Gaussian process (GP) model the true prevalence in order to take the correlation over time into account.
#' The `GP` method uses a squared exponential kernel to characterize the correlation of the GP.
#' This basically assumes that the correlation between the prevalence of two timepoints only depends on the distance between the timepoints (i.e. the time difference).
#' Two parameters are used to model the squared exponential kernel (also known as "exponentiated quadratic kernel"): the lengthscale `lambda` and the output variance `alpha^2`.
#' The lengthscale `lambda` determines the length of the "wiggles".
#' The output variance `alpha^2` is a scale factor that determines the average distance between the realizations of the GP.
#' It characterizes the variance of the prevalence over time.
#'
#' The second (more basic) model (selected by setting `method="timepoint"`) estimates the prevalence at each timepoint without considering any correlation over time.
#' It assumes one parameter for each timepoint in order to estimate the prevalence over time.
#' The `timepoint` model is provided for the sake of comparison but the `GP` model should generally be preferred.
#'
#' **Specificity and sensitivity**
#'
#' The models adjust for imperfect sensitivity and specificity.
#' They assume a default fixed specificity of 99.5% (can be modified through the `spec` argument).
#' The sensitivity is fitted together with the prevalence parameters by using the data reported in Bendavid et al. (2020), which includes three sensitivity studies.
#'
#' **Inference**
#'
#' The models calculate the pool test positivity over time from the prevalence (and the pool sizes `pool.size` given in `data`).
#' It assumes that the number of positive pools `n.pos.pools` follow a binomial distribution with probability being the test positivity.
#'
#' When multiple pooled test results are provided for at least one timepoint (for the same population), the model replace the binomial distribution a beta-binomial distribution.
#' Overdispersion is characterized with the `kappa` parameter.
#'
#' If the `pop` column in the `data` argument is absent, the model assumes only one population (i.e. it models one prevalence over time).
#' If it is provided, the models estimate the prevalence over time for each of the populations present in `data`.
#' This allows for a unique (and simultaneous) estimation of the sensitivity.
#' The user should however consider running the populations separately, if some issues about the convergence of the Stan model are observed.
#'
#' **Priors**
#'
#' The priors can be specified through the `prior` argument, which a named list of vectors containing the hyperparameters of the priors.
#' The names of the elements in the list should match with the names of the parameters.
#'
#' The lengthscale `lambda` has truncated normal prior distribution with default mean of 0 and standard deviation of 2.
#' The output standard deviation `alpha` has a truncated normal prior with default mean of 0 and standard deviation of 2.
#' The prior of the overdispersion `kappa` is obtained by adding `2` to an exponential distribution with default mean of 10.
#' The parameter representing the average prevalence `logit_prev` has logit-normal prior, with default mean of -4 and standard deviation of 2.
#' This corresponds to an average prevalence of 1.8% with 95% confidence intervals (CI) of (0.3%,48%).
#' The sensitivity `sens` has a beta prior distribution whose default hyperparameters are respectively 190 and 40.
#' This corresponds to a mean of approximately 82.6% and with 95% CI of (77.5%, 87.2%).
#'
#'
#'
#' @export
#' @param data Datasets with the following variables:
#'
#' \itemize{
#'   \item time. Time (should be numeric)
#'   \item pop. Population name (if missing, it assumes an unique population to analyse)
#'   \item pool.size. Pool size (number of samples by pool)
#'   \item n.pools. Number of pools
#'   \item n.pos.pools. Number of positive pools.
#' }
#' @param method Method used to analyse data ("timepoint" or "GP")
#' @param spec Test specificity
#' @param prior list of values for prior hyperparameters.
#' @param return.par a logical indicating whether to return GP parameter estimates.
#' @param return.stanfit a logical indicating whether to return the stan model.
#' @param ... Arguments passed to `rstan::sampling` (e.g. iter, chains). If there are divergences,
#'     add `control=list(adapt_delta=0.99)`.
#' @return A list containing the following elements:
#' \itemize{
#'   \item `prev`. A data frame displaying estimated prevalence over time (and populations).
#'   \item `prev_ratio`. Ratio between two consecutive prevalence estimates.
#'   \item `mod_par`. Estimates of some model parameters (returned if `return.par` is `TRUE`).
#'   \item `sampler_par`. Stan diagnosis parameters (returned if `return.par` is `TRUE`).
#'   \item `time`. Stan sampling time (returned if `return.par` is `TRUE`).
#'   \item `stanfit`. `stanfit` object (returned if `return.stanfit` is `TRUE`).
#' }
#' @importFrom dplyr %>%
#' @examples
#'
#' ######################
#' #Example 1
#' #load data
#' data(epipool_data1)
#'
#' #run model with method assuming a prevalence parameter every week
#' out_timepoint = poolprev2(epipool_data1[[1]], method="timepoint", return.par=TRUE,return.stanfit=FALSE) #takes a few seconds
#' #run model with Gaussian process
#' out_GP = poolprev2(epipool_data1[[1]], method="GP", return.par=TRUE,return.stanfit=FALSE) #can take a few minutes
#'
#' #plot
#' library(ggplot2)
#' library(dplyr)
#' dplyr::left_join(epipool_data1[[1]],
#'           rbind(out_timepoint$prev %>% dplyr::mutate(method="timepoint"),
#'                 out_GP$prev %>% dplyr::mutate(method="GP")),
#'           by="time") %>%
#'   ggplot(aes(x=time))+
#'   geom_ribbon(aes(ymin=lwr,ymax=upr),fill="black",alpha=0.1)+
#'   geom_line(aes(y=mean),col="black",size=0.8)+
#'   geom_point(aes(y=mean),col="black",size=2)+
#'   geom_point(aes(y=prev),col="red",size=2) +
#'   facet_grid(.~method) +
#'   theme_bw() +
#'   scale_y_continuous(name="Prevalence", labels = scales::percent)
#'
#' ######################
#' #Example 2: estimate the prevalence of two populations
#' #load data
#' data(epipool_data2)
#' data(epipool_data3)
#' data = rbind(epipool_data2[[1]] %>% dplyr::mutate(pop="pop1"),
#'             epipool_data3[[1]] %>% dplyr::mutate(pop="pop2"))
#'
#' #run model with GP
#' out_timepoint = poolprev2(data, method="timepoint", return.par=TRUE) #takes a few seconds
#' out_GP1 = poolprev2(data, method="GP", return.par=TRUE) #can take a few minutes
#' out_GP2 = poolprev2(data, method="GP",prior=list(lambda=c(0,1),alpha=c(0,1)),
#'                      return.par=TRUE) #can take a few minutes
#'
#' #plot prevalence
#' library(ggplot2)
#' left_join(data,
#'           rbind(out_timepoint$prev %>% dplyr::mutate(method="timepoint"),
#'                 out_GP1$prev %>% dplyr::mutate(method="GP"),
#'                 out_GP2$prev %>% dplyr::mutate(method="GP with smaller lengthscale/sd priors")), by=c("time","pop")) %>%
#' ggplot(aes(x=time))+
#' geom_ribbon(aes(ymin=lwr,ymax=upr),fill="black",alpha=0.1)+
#' geom_line(aes(y=mean),col="black",size=0.8)+
#' geom_point(aes(y=mean),col="black",size=2)+
#' geom_point(aes(y=prev),col="red",size=2) +
#' facet_grid(method~pop) +
#' theme_bw() +
#' scale_y_continuous(name="Prevalence", labels = scales::percent)
#'
#' #compare parameter estimates
#' print(out_timepoint$par)
#' print(out_GP1$par)
#' print(out_GP2$par)


poolprev <- function(data, method="GP",
                      spec=0.995,
                      prior=list(lambda=c(0,2),
                                 alpha=c(0,2),
                                 kappa=10,
                                 logit_prev=c(-4,2),
                                 sens=c(190,40)),
                      return.par=FALSE,
                      return.stanfit=FALSE, ...){
  #check that time is numeric and not a date
  #Prior
  #default prior (same as in argument)
  prior_default=list(lambda=c(0,2), #lengthscale
                      alpha=c(0,2), #GP standard deviation
                      kappa=10, #phi parameter
                      logit_prev=c(-4,2),
                      sens=c(190,40))
  #replace missing element in the user prior by the default prior
  prior = c(prior,prior_default[setdiff(names(prior_default),names(prior))])

  #remove row with na values
  dim.data1 = dim(data)[1]
  if(!any(colnames(data)=="pop")){data$pop="pop1"} #create a pop variable with unique population if missing
  data = data %>%
    dplyr::select(time,pop,n.pools,n.pos.pools,pool.size) %>%
    tidyr::drop_na()
  dim.data2 = dim(data)[1]
  if(dim.data2<dim.data1){print(paste0(dim.data1-dim.data2, " rows dropped due to missing values"))}

  #rank time and calculate number of regions
  N=dim(data)[1] #number of observation
  t.unique = data$time %>% unique() %>% sort() #distinct time points with data
  t.f = t.unique #time points on which GP is calculated (could include some prediction point in addition to t_unique)
  pop_unique = data$pop %>% unique() %>% sort()

  #data
  data = data %>%
    dplyr::mutate(rank.t = dplyr::dense_rank(time),
                  pop = as.numeric(factor(pop,levels=pop_unique)))

  #matrix indication when using overdispersion parameter phi (and if yes, it indicates the position in the phi vector)
  phi_pos = data %>%
    dplyr::group_by(pop,rank.t) %>%
    dplyr::summarise(n.regions=dplyr::n(),.groups="drop_last") %>%
    dplyr::summarise(use_phi = as.numeric(max(n.regions)>1),.groups="drop") %>%
    dplyr::arrange(pop) %>%
    dplyr::mutate(phi_pos  = dplyr::row_number(use_phi)-sum(use_phi==0)) %>% dplyr::pull(phi_pos)
  phi_pos [is.na(phi_pos ) | phi_pos<1] = 0

  #stan data
  standata = list(N = N, #number of observations
                  N_t = length(t.f), #number of distinct time points
                  N_pop = length(pop_unique),
                  rank_t = structure(data$rank.t,dim=N), #rank of the time point related to the observations
                  pop = structure(data$pop,dim=N),
                  t = t.f, #time points where prevalence is calculated (could include prediction for GP)
                  s = structure(data$pool.size,dim=N), #number of samples by pool
                  n = structure(as.integer(data$n.pools),dim=N), #number of pools
                  k = structure(as.integer(data$n.pos.pools),dim=N), #number of positive pools
                  #data sensitivity specificity
                  J_sens = 3,
                  y_sens = c(78,27,25),
                  n_sens = c(85,37,35),
                  spec=spec,
                  #hyperparameters of the priors
                  p_sens = prior$sens,
                  p_intercept = prior$logit_prev,
                  p_lambda = prior$lambda,
                  p_alpha = prior$alpha,
                  #overdispersion
                  N_phi=sum(phi_pos>0),
                  use_phi=structure(phi_pos,dim=length(pop_unique)),
                  p_phi = 1/prior$kappa,
                  #inference
                  inference=1)

  #run the stan model
  if(all(phi_pos==0)){
    if(method=="timepoint"){
      stan <- rstan::sampling(stanmodels$mod_timepoint, data = standata, ...,cores = getOption("mc.cores", 4L))
    }else if(method=="GP"){
      stan <- rstan::sampling(stanmodels$mod_GP, data = standata, ...,cores = getOption("mc.cores", 4L))
    }else{stop("Method should be either timepoint or GP")}
  }else{
    if(method=="timepoint"){
      stan <- rstan::sampling(stanmodels$mod_timepoint_overdisp, data = standata, ...,cores = getOption("mc.cores", 4L))
    }else if(method=="GP"){
      stan <- rstan::sampling(stanmodels$mod_GP_overdisp, data = standata, ...,cores = getOption("mc.cores", 4L))
    }else{stop("Method should be either timepoint or GP")}
  }

  #prevalence estimates
  prev = rstan::summary(stan,par="prev_f")$summary %>%
    tibble::as_tibble() %>%
    dplyr::mutate(time = rep(t.f,length(pop_unique)),
                  pop = rep(pop_unique,each=length(t.f))) %>%
    dplyr::select(pop,time,mean, median=`50%`,lwr=`2.5%`,upr=`97.5%`)

  #sampler parameters
  sampler_par <- rstan::get_sampler_params(stan, inc_warmup = FALSE)
  sampler_par <- do.call(rbind, sampler_par) %>% tibble::as_tibble() %>%
    dplyr::summarise(mean_treepdepth = mean(treedepth__),
                     n_treedepth_10orhigher = sum(treedepth__>=10),
                     n_divergence = sum(divergent__))

  #other parameter estimates
  if(return.par){

    #model par
    mod_par = rstan::summary(stan, par=c("kappa","sens"))$summary  %>%
      tibble::as_tibble() %>%
      dplyr::mutate(par=rep(c("kappa","sens"),each=sum(phi_pos>0)),
                    pop = rep(pop_unique[which(phi_pos>0)],2)) %>%
      dplyr::select(pop,par,mean, median=`50%`,lwr=`2.5%`,upr=`97.5%`)
    if(method=="GP"){
      mod_par = rbind(mod_par,
                  rstan::summary(stan, par=c("lambda","alpha"))$summary %>%
                  tibble::as_tibble() %>%
                  dplyr::mutate(par=rep(c("lambda","alpha"),each=length(pop_unique)),
                                pop = rep(pop_unique,2)) %>%
                  dplyr::select(pop,par,mean, median=`50%`,lwr=`2.5%`,upr=`97.5%`))
    }
    mod_par = mod_par %>% dplyr::left_join(data.frame(par=c("kappa","sens","lambda","alpha"),
                                       desc=c("overdispersion parameter","sensitivity","GP lengthscale","GP sd")),
                            by="par")

    #prevalence ratio
    print("--------------------------------")
    print(rstan::summary(stan)$summary %>% rownames())
    prev_ratio = rstan::summary(stan,par="prev_ratio")$summary %>%
      tibble::as_tibble() %>%
      dplyr::mutate(time = rep(t.f[-1],length(pop_unique)),
                    pop = rep(pop_unique,each=length(t.f)-1)) %>%
      dplyr::select(pop,time,mean, median=`50%`,lwr=`2.5%`,upr=`97.5%`)

    #sampler parameters
    sampler_par <- rstan::get_sampler_params(stan, inc_warmup = FALSE)
    sampler_par <- do.call(rbind, sampler_par) %>% tibble::as_tibble() %>%
      dplyr::summarise(mean_treepdepth = mean(treedepth__),
                       n_treedepth_10orhigher = sum(treedepth__>=10),
                       n_divergence = sum(divergent__))

    #sampling time
    time = rstan::get_elapsed_time(stan) %>% tibble::as_tibble() %>%  dplyr::pull(sample) %>% max()

    out = list(prev=prev, prev_ratio = prev_ratio, mod_par = mod_par, sampler_par = sampler_par, time = time)
  }else{
    out = list(prev=prev)
  }
  if(return.stanfit){
    out = c(out, stanfit=stan)
  }

  return(out)
}
