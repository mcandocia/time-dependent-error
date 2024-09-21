# test multiple sets of data
# smaller sizes can get more trials
# larger sizes fewer trials
library(tidyverse)
library(rstan)
library(TSA)

options(mc.cores=4)

#stan_model_compiled = NA
compiled_stan_models = list()

RETURN_PARS_2 =   c('beta0','beta1','beta2','beta3', 'beta4', 'sigma','beta_t','beta_c1','beta_c2','error_0',
    'state_drift_scale','error_0','var_0','vs_k','vs_a','var_0_sigma0','var_decay_rate')

decay_cumsum <- function(x, rate,v0=0, v0_rate = FALSE){
  N = length(x)
  y = numeric(N)
  y[1] = x[1]
  if (v0 != 0){
    y[1] = y[1] + v0 * ifelse(v0_rate,rate,1)
  }
  for (i in 2:N){
    y[i] = y[i-1] * rate + x[i]
  }
  return(y)
}

split_groups <- function(x){
  g = as.numeric(x)
  d = diff(g)
  w = which(d != 0)+1
  starts = c(1, w)
  ends = c(w-1, length(g))
  return(data.frame(s=starts,e=ends,Value=x[starts]))
}


run_model <- function(
  # padding params
  n_legs = 5,
  padding_control = c(2,1/5),
  leg_size_control = c(25,25,2,0.5,1,1),
  # note, if control is set to 0,0,X,Y,0,Z
  # leg sizes will be leg_size_constant + 1
  leg_size_constant = 0,
  padding_constant = 0,

  # time control
  ar_coef = 0.95,
  ar_alpha = 50,
  ar_beta = 1.2,

  # category
  nc_v1 = 10,
  nc_v2 = 5,

  nc_v1_splits = 2,
  c_v1_randomized = c('no','partial','yes'),

  target_time_sd=25,
  standard_noise_factor=1,

  showPlot=TRUE,

  time_noise_type=c('ar','sine'),
  sine_params = list(
    n_periods = 1,
    phase_shift = 0,
    rectify_magnitude=1
  ),
  single_leg_size = 1000,
  skip_pause=FALSE,
  stan_filename='tsa_multileg_test_j.stan',
  n_chains=4,
  return_pars = c('beta0','beta1','beta2','beta3', 'beta4', 'sigma','sigma_t','beta_t','beta_c1','beta_c2','error_0'),
  extra_data = list(),
  force_recompile=FALSE,
  plotOnly=FALSE,
  ... # extra stan parameters
){

  c_v1_randomized = c_v1_randomized[1]

  if (n_legs > 1){
    arima_padding = 1+round(runif(n_legs-1) * padding_control[1] +
                              rexp(n_legs-1,padding_control[2])) +
      padding_constant
    leg_sizes = 1+round(runif(n_legs)^leg_size_control[3] * leg_size_control[1] +
                          runif(n_legs)^leg_size_control[4] * leg_size_control[2] +
                          leg_size_control[5] * rexp(n_legs,leg_size_control[6])) +
      leg_size_constant
    N=sum(leg_sizes)
  } else {
    leg_sizes = single_leg_size
    arima_padding = 0
    N = leg_sizes
  }
  # calculate whole AR series, then cut out undesired parts
  time_noise_type = time_noise_type[1]

  if (time_noise_type == 'ar'){
    time_error = arima.sim(list(ar=ar_coef),n=N+sum(arima_padding))
    time_error_sd = 1/sqrt(1-ar_coef^2)
  } else {
    N_eff = N+sum(arima_padding)
    time_error = sin(
      (0:(N_eff-1))/N_eff * 2 * pi *sine_params$n_periods +
        sine_params$phase_shift
    )
    time_error = sign(time_error) * pmin(abs(time_error), sine_params$rectify_magnitude)
    # note: rectified standard deviation is theoretically lower, but
    # it's reasonably constrained to a square wave, which has sqrt(2)x more
    # after scaling. Also, rectifying shouldn't be too high, just enough to show a flat
    # "cap" can be accounted for
    time_error_sd = pmin(1,sine_params$rectify_magnitude) / sqrt(2)
  }

  if (n_legs > 1){
    cut_starts = 1+cumsum(leg_sizes %>% head(-1)) +
      c(0,cumsum(arima_padding %>% head(-1)))

    cut_ends = cut_starts+arima_padding-1
    cut_bands = unlist(lapply(1:(n_legs-1), function(x) cut_starts[x]:cut_ends[x]))
    time_error2 = time_error[-cut_bands]
    time_diffs = rep(1,N)
    time_diffs[1+cumsum(leg_sizes %>% head(-1))] = arima_padding
  } else {
    time_diffs = rep(1,N)
    time_error2 = time_error
  }

  v1e = 40*(runif(nc_v1) - 0.5)
  v2e = rexp(nc_v2, 1/10)


  time_error_factor = target_time_sd/time_error_sd

  c_v1_randomized = tolower(c_v1_randomized)
  if (c_v1_randomized == 'no' | c_v1_randomized == FALSE){
    v1c = rep(rep(1:nc_v1,nc_v1_splits),each=round(N/nc_v1/nc_v1_splits+0.5)) %>% head(N)
  } else if (c_v1_randomized=='partial' | c_v1_randomized == TRUE){
    #randomized
    split_positions = rank(rnorm(nc_v1 * nc_v1_splits))
    cats = rep(1:nc_v1,nc_v1_splits)[split_positions]
    v1c = rep(cats, each=round(N/nc_v1/nc_v1_splits+0.5)) %>% head(N)
  } else {
    v1c = (1:nc_v1)[round(runif(N,0.5,nc_v1+0.49999999))]
  }
  v2c = rep(1:nc_v2, round(N/nc_v2)+1) %>% head(N)


  # add a sort() to introduce long-term bias with sampling
  # as well as x4
  tdf = tibble(
    x1 = (rnorm(N)*20) + rnorm(N) * 10,
    x2 = rnorm(N)*30,
    x3 = sqrt(abs(x1*x2)) + rnorm(N) * 10,
    x4 = rbinom(N, 1, (1+1:N/20)/(1+N)),
    c1 = v1c,
    c2 = v2c,
    y0 = v1e[c1] + v2e[c2] +
      x1 + 2 * x2 + x3 -3 * x4 +
      rnorm(N)*standard_noise_factor,
    error_t = time_error2*time_error_factor,
    y = y0 + error_t,
    td = time_diffs
  )

  stan_data = list(
    N=N,
    x1 = tdf$x1,
    x2 = tdf$x2,
    x3 = tdf$x3,
    x4 = tdf$x4,
    k1 = nc_v1,
    k2 = nc_v2,
    c1 = model.matrix(~c1,data=tdf %>% mutate(c1=factor(c1)))[,-1],
    c2 = model.matrix(~c2,data=tdf %>% mutate(c2=factor(c2)))[,-1],
    y = as.numeric(tdf$y),
    td = time_diffs,
    ar_alpha=ar_alpha,
    ar_beta=ar_beta
  )

  if (length(extra_data) > 0){
    dn = names(extra_data)
    for (nm in dn){
      stan_data[[nm]] = extra_data[[nm]]
    }
  }

  # plotting

  csplits = split_groups(tdf$c1)
  y_quants = lapply(1:nrow(csplits), function(i){
    s = csplits[i,'s']
    e = csplits[i,'e']
    tdf %>%
      slice(s:e) %>%
      summarize(
        s = s,
        e=e,
        yl = quantile(y,0.25),
        ymed = quantile(y,0.5),
        ymean = mean(y),
        yu = quantile(y,0.75)
      )
  }
  ) %>% bind_rows()

  group_averages = tdf %>% group_by(c1) %>% summarize(ya = mean(y))

  csplits = csplits %>%
    left_join(
      group_averages %>% rename(Value=c1)
    )

  yiqr = max(y_quants$yu) - min(y_quants$yl)

  sr_adjustment = ifelse(nrow(csplits) < 40, 0, 1.6*(nrow(csplits) - 40)/nrow(csplits))
  print(sprintf('SRA: %0.2f', sr_adjustment))
  print(head(csplits))

  # visualize the effect of groupings being time-dependent
  plt <- ggplot(tdf %>% mutate(x=1:n(),c1=factor(c1),c2=factor(c2))) +
    geom_rect(data=csplits,
              aes(xmin=s,xmax=e,color=as.factor(Value),
                  fill=as.factor(Value),
                  ymin=min(c(y_quants$yl,tdf$error_t))-0.05*yiqr,ymax=max(c(y_quants$yu,tdf$error_t))+0.05*yiqr),alpha=0.6) +
    geom_line(aes(x=x,y=error_t)) +
    geom_line(data=y_quants,aes(x=(s+e)/2,y=yu),linetype='dashed') +
    geom_line(data=y_quants,aes(x=(s+e)/2,y=yl),linetype='dashed') +
    scale_color_manual('Category', values=cetcolor::cet_pal(nc_v1,'r3')) +
    scale_fill_manual('Category', values=cetcolor::cet_pal(nc_v1,'r3')) +
    geom_point(data=csplits, aes(x=(s+e)/2,y=ya),color='red',size=2.5-sr_adjustment) +
    geom_point(data=y_quants, aes(x=(s+e)/2,y=ymean),color='blue',size=2.5-sr_adjustment,shape=15) +
    theme_bw() +
    guides(size='none',color='none') +
    xlab('Measurement #') +
    ylab('y') +
    ggtitle('Time Dependent Error Across Time-Dependent Categories',
            subtitle='Dashed lines: 25th and 75th percentile of Y in range\nRed dots: average y for category\nBlue squares: average y for category within contiguous zone')

  if (showPlot){
    print(plt)
  }

  if (plotOnly){
    return(list(
      plt = plt,
      tdf=tdf,
      sr_adjustment=sr_adjustment,
      csplits=csplits,
      group_averages=group_averages,
      N=N,
      time_diffs=time_diffs,
      v1e=v1e,
      v2e=v2e,
      y_quants=y_quants,
      nc_v1=nc_v1,
      nc_v2=nc_v2,
      time_error2=time_error2,
      time_error_factor=time_error_factor,
      time_error_sd=time_error_sd,
      standard_noise_factor=standard_noise_factor,
      ar_coef=ar_coef
    ))
  }


  if (!skip_pause){
    print('waiting 8 seconds...')
    Sys.sleep(8)
  }

  # this doesn't work, but it doesn't make a difference rn
  # TODO: fix and add "force model remake" parameter
  if (!stan_filename %in% names(compiled_stan_models) & !force_recompile){
    print('need to compile model 1st time')
    compiled_stan_models[[stan_filename]] <<- stan_model(stan_filename)
  }
  t1=Sys.time()
  print(N)
  print(t1)
  stan_model = rstan::sampling(
    compiled_stan_models[[stan_filename]],
    data = stan_data,
    pars = return_pars,
    chains=n_chains,
    ...
  )
  t2 = Sys.time()
  print(N)
  print(t2-t1)

  # comparison of "true" vs. linear model vs. Stan model

  lm_mod = lm(y~x1+x2+x3+x4+c1+c2,
                data=tdf %>% mutate(c1=factor(c1),c2=factor(c2)))

  # include intercept
  lm_coefs = c(0,coef(lm_mod)[paste0('c1',2:nc_v1)]) + coef(lm_mod)[1]
  stan_coefs_all = get_posterior_mean(stan_model)
  stan_coefs = stan_coefs_all[grepl('\beta_c1\\[',rownames(stan_coefs_all)),] %>%
    as.data.frame()

  stan_coefs = bind_rows(
    matrix(0,nrow=1,ncol=ncol(stan_coefs)) %>%
      as.data.frame() %>%
      setNames(colnames(stan_coefs)),
    as.data.frame(stan_coefs)
  )

  c1_values = tibble(
    n = 1:nc_v1,
    actual=v1e,
    linear=lm_coefs
  ) %>%
    bind_cols(
     stan_coefs
    ) %>%
    mutate_all(
      function(x) x - x[1]
    )

  return(
    list(
      stan_model=stan_model,
      linear_model=lm_mod,
      data=tdf,
      comparison=c1_values,
      N=N,
      plt=plt,
      time_elapsed=as.numeric(difftime(t2,t1,units='mins')),
      cats = list(c1=v1e,c2=v2e),
      y_quants=y_quants,
      group_averages=group_averages,
      posterior = get_posterior_mean(stan_model),
      sigma_t_true = time_error_factor,
      standard_noise_factor=standard_noise_factor,
      ar_coef = ar_coef
    )
  )
}

get_par <- function(model, pattern, index=1){
  po = model$posterior

  return(po[grepl(pattern,rownames(po)),index])
}

add_results <- function(model, comparison='average'){
  plots = list()
  plot_data = list()
  # index
  if (comparison == 'average'){
    #print('avg')
    j = ncol(model$posterior)
    jc = 'stan_all'
  } else {
    #print('sng')
    j = comparison
    jc = sprintf('stan_%s',j)
  }

  # linear vs. actual & stan vs. actual
  yhat = get_par(model,'y_hat',j)
  yhat_lm = predict(model$linear_model,newdata=model$data %>%
                      mutate(c1=factor(c1),c2=factor(c2))
  )

  stan_coefs = get_par(model,'beta_c1\\[.*',1:ncol(model$posterior))

  prediction_df = tibble(
    actual = model$data$y,
    stan=yhat,
    linear=yhat_lm
  )
  prediction_df_l = prediction_df %>%
    pivot_longer(c('stan','linear'))

  # line fits
  lml = lm(actual~linear,data=prediction_df)
  lms = lm(actual~stan,data=prediction_df)

  lmlab = coef(lml)
  lmsab = coef(lms)

  plot_data$prediction_df = prediction_df
  plot_data$prediction_df_l = prediction_df_l

  (plt<-ggplot(prediction_df_l) +
    geom_point(
      aes(x=actual,y=value,color=name),
      alpha=0.4
    ) +
      theme_bw() +
      geom_abline(slope=1,intercept=0,linetype='dashed',color='#11FF11',size=1.5,alpha=0.7) +
      xlab('True Value') +
      ylab('Model Value') +
      guides(color=guide_legend('Model')) +
      ggtitle('Predictions of Linear Model vs. Stan Model',
              subtitle='Stan model estimates slow-moving time series') +
      scale_color_manual(values=c('linear'='#FF3333',stan='#2222BB'))
  )

  plots[['yhat_comp']] = plt

  r1 = rep(0,ncol(stan_coefs))
  # parameter comparison
  co = model$comparison[,1:3]
  co2 = co %>%
    bind_cols(
      stan_coefs %>% as.data.frame() %>% add_row(.before=1)
    ) %>%
    rename_all(
      function(x) gsub('mean-chain:(\\1)','stan_\\1',x)
    ) %>%
    rename(
      stan_all = `mean-all chains`
    ) %>%
    mutate_all(function(x)coalesce(x,0))

  co3 = co2 %>%
    pivot_longer(
      starts_with('stan_') | matches('linear')
    )

  plot_data$co2 = co2
  plot_data$co3 = co3

  (plt <-ggplot(
    co3 %>% filter(name %in% c('linear',jc)) %>%
      mutate(
        name=ifelse(name=='linear',name,'stan')
        )
    )  +
      geom_point(
        aes(x=actual,y=value,color=name),
        alpha=0.8
      ) +
      theme_bw() +
      geom_abline(slope=1,intercept=0,linetype='dashed',color='#11FF11',size=1.5,alpha=0.7) +
      xlab('True Value') +
      ylab('Model Value') +
      guides(color=guide_legend('Model')) +
      ggtitle('Parameter Estimates of Linear Model vs. Stan Model',
              subtitle='Stan model estimates slow-moving time series') +
      scale_color_manual(values=c('linear'='#FF3333',stan='#2222BB'))
  )

  plots[['param_comp']] = plt

  # time error estimate vs. actual error
  error_t = model$data$error_t
  error_t_estimate = get_par(
    model,
    'error_t_hat'
  )

  time_error_df = tibble(
    time=1:length(error_t),
    actual = error_t - mean(error_t),
    estimate= error_t_estimate - mean(error_t_estimate),
  ) %>% pivot_longer(
    c('actual','estimate')
  )

  (plt <- ggplot(time_error_df) +
      geom_line(
        aes(x=time,y=value,color=name),
        alpha=0.6,
        size=1.1
      ) +
      scale_color_manual(
        values=c(actual='#111',estimate='#F33')
      ) + theme_bw() +
      xlab('Time') +
      ylab('Time-Drift Error') +
      guides(color=guide_legend('Source')) +
      ggtitle(
        'Time-Drift Error Estimate vs. Actual',
        subtitle='both series re-centered at 0'
      )
  )

  plots[['time_error']] = plt
  plot_data$time_error_df = time_error_df

  # add stuff to model
  results = list()

  results$plots = plots
  results$plot_data = plot_data

  model$results = results
  #

  return(model)
}


# I think there's a bug in these
get_error_estimate_cv <- function(model, index=1, plot=TRUE){
  mod = get_posterior_mean(model$stan_model)
  j = index
  mod_de = mod[grepl('drift_errors\\[',rownames(mod)),]
  mod_dv = mod[grepl('drift_var\\[',rownames(mod)),]
  bt = mod[grepl('beta_t$',rownames(mod)),]
  vd = mod[grepl('var_decay_rate',rownames(mod)),]
  v0 = mod[grepl('var_0_sigma0',rownames(mod)),]

  variance = v0[j]*exp(decay_cumsum(mod_dv[,j],vd[j] ))
  error = decay_cumsum(mod_de[,j] * sqrt(variance) , bt[j])

  if (plot){
    plot(model$data$error_t, col='red',type='l')
    lines(error,
          type='l'
    )
  }
  return(error)
}

get_error_estimate <- function(model, index=1, plot=TRUE){


  mod = get_posterior_mean(model$stan_model)

  mod_de = mod[grepl('drift_errors\\[',rownames(mod)),]
  bt = mod[grepl('beta_t$',rownames(mod)),]
  st = mod[grepl('sigma_t$',rownames(mod)),]
  j = index
  error = decay_cumsum(mod_de[,j] * st[j], bt[j])

  yc1 = c(0,mod[grepl('beta_c1\\[.+$',rownames(mod)),j])
  c1e_est = yc1[model$data$c1]
  c1e = c(0,model$cats$c1)[model$data$c1]

  c1e_diff = c1e - c1e_est
  #print(sum(table(table(c1e_diff))))

  # yhat = t_error_hat + y0_hat
  # y0_hat ~ y0 + (c1e_hat-c1e)

  if (plot){
    plot(model$data$error_t, col='red',type='l')
    lines(error,
      type='l'
    )
    #lines(model$data$y, col='#2222FF33', type='l')
    #lines(c1e_diff+error, col='green',type='l')
  }
  return(error)
}

get_numeric_results <- function(model,j=1){
  N = model$N
  time_elapsed = model$time_elapsed
  if (j == ncol(model$posterior)){
    comparison = 'average'
  } else {
    comparison = j
  }
  if (! 'results' %in% names(model)){
    model = add_results(model,comparison)
  }
  sigma = get_par(model,'sigma$',j)
  if (!'standard_noise_factor' %in% names(model)){
    sigma_actual = 1
  } else {
    sigma_actual = model$standard_noise_factor
  }

  if (!'ar_coef' %in% names(model)){
    ar_coef = 0.9995
  } else {
    ar_coef = model$ar_coef
  }

  ar_coef_est = get_par(model,'beta_t$',j)

  sigma_t_true = model$sigma_t_true
  sigma_t = get_par(model,'sigma_t$',j)

  # prediction
  prediction_mse = model$results$plot_data$prediction_df_l %>%
    group_by(name) %>%
    summarize(
      MSE = mean((value-actual)^2),
      R2 = 1-MSE/var(actual)
    ) %>%
    ungroup()

  mse_stan = prediction_mse %>%
    filter(name=='stan') %>% pull(MSE)

  mse_linear = prediction_mse %>%
    filter(name=='linear') %>% pull(MSE)

  # time error
  time_results = model$results$plot_data$time_error_df %>%
    pivot_wider(id_cols='time') %>%
    summarize(
      MSE = mean((actual-estimate)^2),
      R2 = MSE/var(actual)
    )

  time_mse = time_results$MSE
  time_r2 = time_results$R2

  # param accuracy
  if (comparison=='average'){
    target = 'stan_all'
  } else {
    target = sprintf('stan_%s',j)
  }
  param_acc = model$results$plot_data$co3 %>%
    filter(
      name %in% c(target,'linear')# & n != 0
    ) %>%
    group_by(name) %>%
    summarize(
      MSE = compare_cats(value,actual)$mse,
      R2 = compare_cats(value,actual)$r2
    )

  param_mse_stan = param_acc %>%
    filter(name==target) %>% pull(MSE)

  param_mse_linear = param_acc %>%
    filter(name=='linear') %>% pull(MSE)

  param_r2_stan = param_acc %>%
    filter(name==target) %>% pull(R2)

  param_r2_linear = param_acc %>%
    filter(name=='linear') %>% pull(R2)


  return(
    tibble(
      N=N,
      time_elapsed=time_elapsed,
      sigma=sigma,
      sigma_true=sigma_actual,
      sigma_t=sigma_t,
      sigma_t_true=sigma_t_true,
      yhat_mse_stan = mse_stan,
      yhat_mse_linear = mse_linear,
      time_mse = time_mse,
      time_r2 = time_r2,
      param_mse_stan=param_mse_stan,
      param_mse_linear=param_mse_linear,
      param_r2_stan = param_r2_stan,
      param_r2_linear = param_r2_linear,
      ar_coef = ar_coef,
      ar_coef_est = ar_coef_est
    )
  )

}

# general parameter controls for size:
#n_legs = 5, # larger numbers make sense for more data
#padding_control = c(2,1/5), Uniform and Exponential component
#leg_size_control = c(15,15,2,2,1,1), PF,PF,PE,PE,EF,ES
# PF = power factor; PE = power exponent (of unif.)
# EF = exponential factor; ES = exponential scale

# time control
#ar_coef = 0.95; sd will be 1/(sqrt(1-0.ar_coef^2))
#ar_alpha = 50; this should be higher with higher AR; a/(a+b) should be in ballpark
#ar_beta = 1.2; this can be constant

run_model_with_params<-function(l){
  do.call(run_model, l)
}

run_small_model <- function(...){
  params = list(
    # padding params
    n_legs = 5,
    padding_control = c(2,1/5),
    leg_size_control = c(28,28,2,0.5,3,1/5),
    # time control
    ar_coef = 0.95,
    ar_alpha = 50,
    ar_beta = 1.2,
    nc_v1 = 10,
    nc_v2 = 5,
    ...
  )
  run_model_with_params(params)
}

run_medium_model<- function(...){
  params = list(
    # padding params
    n_legs = 50,
    padding_control = c(50,1/20),
    leg_size_control = c(55,80,2,0.7,1,1),
    # time control
    ar_coef = 0.9975,
    ar_alpha = 80,
    ar_beta = 1.2,
    nc_v1 = 10,
    nc_v2 = 5,
    ...
  )
  run_model_with_params(params)
}

run_large_model<- function(...){
  params = list(
    # padding params
    n_legs = 250,
    padding_control = c(250,1/200),
    leg_size_control = c(300,60,4,2,5,1/10),
    # time control
    ar_coef = 0.998,
    ar_alpha = 120,
    ar_beta = 1.2,
    nc_v1 = 10,
    nc_v2 = 5,
    ...
  )
  run_model_with_params(params)
}

run_medium_sine_model <- function(...){
  params = list(
    leg_size_control = c(0,0,1,1,0,1),
    leg_size_constant = 799,
    n_legs = 5,
    padding_control=c(5,1/10),
    padding_constant = 2,
    nc_v1 = 10,
    nc_v2 = 5,
    time_noise_type='sine',
    sine_params = list(
      n_periods = 1.33,
      phase_shift = 0,
      rectify_magnitude=0.85
    ),
    ar_alpha = 150,
    ar_beta = 1.2,
    ...
  )
  run_model_with_params(params)
}

run_small_vdrift_model <- function(...){
  return_pars = c('beta0','beta1','beta2','beta3', 'beta4', 'sigma','beta_t','beta_c1','beta_c2','error_0',
                  'state_drift_scale','error_0','var_0','vs_k','vs_a','var_0_sigma0','var_decay_rate')
  params = list(
    leg_size_control = c(0,0,1,1,0,1),
    leg_size_constant = 199,
    n_legs = 5,
    padding_control=c(5,1/4),
    padding_constant = 2,
    nc_v1 = 12,
    nc_v2 = 5,
    time_noise_type='sine',
    sine_params = list(
      n_periods = 1.33,
      phase_shift = 0,
      rectify_magnitude=0.85
    ),
    ar_alpha = 120,
    ar_beta = 1.2,
    stan_filename='stan_ar_ml_ak_test.stan',
    n_chains=6,
    extra_data = list(
      sigma_var_lower = 0.1,
      sigma_var_upper = 8
    ),
    return_pars = return_pars,
    ...
  )
  run_model_with_params(params)
}

run_medium_vdrift_model <- function(...){
  return_pars = c('beta0','beta1','beta2','beta3', 'beta4', 'sigma','beta_t','beta_c1','beta_c2','error_0',
    'state_drift_scale','error_0','var_0','vs_k','vs_a','var_0_sigma0','var_decay_rate')
  params = list(
    leg_size_control = c(0,0,1,1,0,1),
    leg_size_constant = 899,
    n_legs = 5,
    padding_control=c(5,1/10),
    padding_constant = 2,
    nc_v1 = 12,
    nc_v2 = 5,
    time_noise_type='sine',
    sine_params = list(
      n_periods = 1.33,
      phase_shift = 0,
      rectify_magnitude=0.85
    ),
    ar_alpha = 150,
    ar_beta = 1.2,
    stan_filename='stan_ar_ml_ak_test.stan',
    n_chains=6,
    extra_data = list(
      sigma_var_lower = 0.1,
      sigma_var_upper = 8
    ),
    return_pars = return_pars,
    ...
  )
  run_model_with_params(params)
}

TRIAL_DIRECTORY = 'stan_ar_trials'
dir.create(TRIAL_DIRECTORY,showWarnings=FALSE)



crossdiff <- function(x){
  N = length(x)
  f = expand.grid(1:N,1:N) %>%
    mutate(f=Var1<Var2) %>%
    pull(f)
  #print(expand.grid(x,x))
  (expand.grid(x,x) %>%
    mutate(
      diff = Var1-Var2
    ))[
      f,
    ] %>%
    pull(diff)
}

compare_cats <- function(x,y){
  # 🐈 😺
  dx = crossdiff(x)
  dy = crossdiff(y)
  mse = mean((dx-dy)^2)
  r2 = 1-mse/var(y)
  data.frame(mse=mse,r2=r2)
}
