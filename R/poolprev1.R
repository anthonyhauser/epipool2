# Save this file as `R/poolprev1.R`

#' Estimate prevalence from pooled test results assuming a fixed prevalence parameter for each time point
#'
#' @export
#' @param data Datasets with the following variables:
#'
#' \itemize{
#'   \item time. Time (should be numeric)
#'   \item pool.size. Pool size (number of samples by pool)
#'   \item n.pools. Number of pools
#'   \item n.pos.pools. Number of positive pools.
#' }
#' @param method Method used to analyse data (timepoint or GP)
#' @param spec Test specificity
#' @param prior list of values for hyperprior parameters.
#' @param ... Arguments passed to `rstan::sampling` (e.g. iter, chains).
#' @return An object of class `stanfit` returned by `rstan::sampling`
#' @importFrom dplyr %>%
#' @examples
#' #load data
#' data(epipool_data1)
#'
#' #run model with method assuming a prevalence parameter every week
#' prev_timepoint = poolprev1(epipool_data1[[1]],method="timepoint") #takes a few seconds
#' #run model with GP
#' prev_GP = poolprev1(epipool_data1[[1]],method="GP",control=list(adapt_delta=0.99,max_treedepth=12)) #can take a few minutes
#'
#' #plot
#' library(ggplot2)
#' left_join(epipool_data1[[1]],
#' rbind(prev_timepoint %>% dplyr::mutate(method="timepoint"),
#'       prev_GP %>% dplyr::mutate(method="GP")),
#'       by="time") %>%
#' ggplot(aes(x=time))+
#' geom_ribbon(aes(ymin=lwr,ymax=upr),fill="black",alpha=0.1)+
#' geom_line(aes(y=mean),col="black",size=0.8)+
#' geom_point(aes(y=mean),col="black",size=2)+
#' geom_point(aes(y=prev),col="red",size=2) +
#' facet_grid(.~method) +
#' theme_bw()

poolprev1 <- function(data, method="GP",
                      spec=0.995,
                      prior=list(lengthscale=c(0,1),
                                 sd_GP=c(0,1),
                                 logit_prev=c(-4,2)), ...){
  #check that time is numeric and not a date

  #remove row with na values
  dim.data1 = dim(data)[1]
  data = data %>%
    dplyr::select(time,n.pools,n.pos.pools,pool.size) %>%
    tidyr::drop_na()
  dim.data2 = dim(data)[1]
  if(dim.data2<dim.data1){print(paste0(dim.data1-dim.data2, " rows dropped due to missing values"))}

  #rank time and calculate number of regions
  data = data %>%
    dplyr::mutate(rank.t = dplyr::dense_rank(time))
  n.regions = data %>%
    dplyr::group_by(rank.t) %>%
    dplyr::summarise(n=n(),.groups="drop") %>%
    dplyr::pull(n) %>% max()

  #data dimension and unique time points
  N=dim(data)[1] #number of observation
  t.unique = data$time %>% unique() %>% sort() #distinct time points with data
  t.f = t.unique #time points on which GP is calculated (could include some prediction point in addition to t_unique)

  #stan data
  standata = list(N = N, #number of observations
                  N_t = length(t.f), #number of distinct time points
                  rank_t = structure(data$rank.t,dim=N), #rank of the time point related to the observations
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
                  p_intercept = prior$logit_prev,
                  p_phi = 1/100,
                  p_lambda = prior$lengthscale,
                  p_alpha = prior$sd_GP,
                  #inference
                  inference=1)

  #run the stan model
  if(n.regions==1){
    if(method=="timepoint"){
      out <- rstan::sampling(stanmodels$mod1_timepoint, data = standata, ...,cores = getOption("mc.cores", 4L))
    }else if(method=="GP"){
      out <- rstan::sampling(stanmodels$mod1_GP, data = standata, ...,cores = getOption("mc.cores", 4L))
    }else{stop("Method should be either timepoint or GP")}
  }else{
    if(method=="timepoint"){
      out <- rstan::sampling(stanmodels$mod2_timepoint, data = standata, ...,cores = getOption("mc.cores", 4L))
    }else if(method=="GP"){
      out <- rstan::sampling(stanmodels$mod2_GP, data = standata, ...,cores = getOption("mc.cores", 4L))
    }else{stop("Method should be either timepoint or GP")}
  }

  prev = rstan::summary(out,par="prev_f")$summary %>%
    as_tibble() %>%
    dplyr::mutate(time= t.f) %>%
    dplyr::select(time,mean, median=`50%`,lwr=`2.5%`,upr=`97.5%`)

  return(prev)
}
