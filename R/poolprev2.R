# Save this file as `R/poolprev2.R`

#' Estimate prevalence from pooled test results
#'
#' `poolprev` is used to estimate prevalence from pooled test results over time.
#' It fits a Bayesian model using the Stan framework.
#'
#' Two models are implemented, which differ in the way the correlation over time is handled.
#' In the default model (selected by setting `method="GP"`), a hidden Gaussian process (GP) model the true prevalence in order to take the correlation over time into account.
#' Another more basic model is also provided (by setting `method="timepoint`), which estimates the prevalence at each timepoint and does not consider any correlation over time.
#' The default model should be preferred.
#'
#' The model adjusts for imperfect sensitivity and specificity.
#' It assumes a default specificity of 99.5% (can be changed through the `spec` argument).
#' The sensitivity is fitted together with the pooled test data by using the data reported in Bendavid et al. (2020), which includes 3 sensitivity studies.
#'
#' When setting ``
#'
#'  `poolprev`  using `rstan` package.
#' Two models are implemented.
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
#' @param method Method used to analyse data (timepoint or GP)
#' @param spec Test specificity
#' @param prior list of values for prior hyperparameters.
#' @param return.par a logical indicating whether to return GP parameter estimates.
#' @param return.stanfit a logical indicating whether to return the stan model.
#' @param ... Arguments passed to `rstan::sampling` (e.g. iter, chains). If there are divergences,
#'     add `control=list(adapt_delta=0.99)`.
#' @return A list containing the following elements:
#' \itemize{
#'   \item prev. A data frame displaying estimated prevalence over time, regions (and populations).
#'   \item pop. Population name (if missing, it assumes an unique population to analyse)
#'   \item pool.size. Pool size (number of samples by pool)
#'   \item n.pools. Number of pools
#'   \item n.pos.pools. Number of positive pools.
#' }
#' @importFrom dplyr %>%
#' @examples

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
#' out_GP2 = poolprev2(data, method="GP",prior=list(lengthscale=c(0,1),sd_GP=c(0,1)),
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


poolprev2 <- function(data, method="GP",
                      spec=0.995,
                      prior=NULL,
                      return.par=FALSE,
                      return.stanfit=FALSE, ...){
  #check that time is numeric and not a date
  #adapt prior to argument
  prior0=list(lengthscale=c(0,2),
             sd_GP=c(0,2),
             phi=1/10,
             logit_prev=c(-4,2))
  if(length(prior)>0){for(i in 1:length(prior)){
    prior0[[which(names(prior0)==names(prior)[i])]] = prior[[i]]
  }}

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
                  p_sens = c(190,40),
                  p_intercept = prior0$logit_prev,
                  p_lambda = prior0$lengthscale,
                  p_alpha = prior0$sd_GP,
                  #overdispersion
                  N_phi=sum(phi_pos>0),
                  use_phi=structure(phi_pos,dim=length(pop_unique)),
                  p_phi = prior0$phi,
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
