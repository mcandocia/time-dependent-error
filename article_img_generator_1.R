source('stan_ar_test_run3.R')
source('stan_ar_test_organizer.R')
library(gridExtra)

PLOT_DIRECTORY=file.path('plots','article1')
dir.create(PLOT_DIRECTORY,showWarnings=FALSE,recursive=TRUE)

fpng <- function(x,height=720,width=720,...){
  png(file.path(PLOT_DIRECTORY,x),...)
}

set.seed(2001)
example_sample = run_model_with_params(list(
  # padding params
  n_legs = 1,
  padding_control = c(20,1/20),
  leg_size_control = c(0,0,2,0.7,0,1),
  leg_size_constant=1000,
  # time control
  ar_coef = 0.9995,
  ar_alpha = 120,
  ar_beta = 2.2,
  nc_v1 = 10,
  nc_v2 = 5,
  c_v1_randomized='no',
  nc_v1_splits = 1,
  seed = seed,
  showPlot=FALSE,
  plotOnly=TRUE
))

tdf = example_sample$tdf  %>%
  mutate(t=1:n())

# show basic drift
fpng('basic_time_series.png',res=120)
print(
  ggplot(tdf %>% mutate(t=1:n())) +
    geom_line(aes(x=t,y=error_t)) +
    theme_bw() +
    ylab('Drift Error') +
    ggtitle('AR(1) Time Series with φ=0.9995')
)
dev.off()


# show with categories (simple model)

v1e = example_sample$v1e
nc_v1 = length(example_sample$v1e)
tdf = tdf %>%
  mutate(
    y_basic = error_t + v1e[c1] + rnorm(n())*1
  )
csplits = split_groups(tdf$c1)
y_quants = lapply(1:nrow(csplits), function(i){
  s = csplits[i,'s']
  e = csplits[i,'e']
  tdf %>%
    slice(s:e) %>%
    summarize(
      s = s,
      e=e,
      yl = quantile(y_basic,0.25),
      ymed = quantile(y_basic,0.5),
      ymean = mean(y_basic),
      yu = quantile(y_basic,0.75)
    )
}
) %>% bind_rows() %>%
  mutate(v1e = v1e, direction=ifelse(v1e>0,'up','down'))
group_averages = tdf %>% group_by(c1) %>% summarize(ya = mean(y))
csplits = csplits %>%
  left_join(
    group_averages %>% rename(Value=c1)
  )
yiqr = max(y_quants$yu) - min(y_quants$yl)
sr_adjustment = ifelse(nrow(csplits) < 40, 0, 1.6*(nrow(csplits) - 40)/nrow(csplits))

#plot
plt_basic2 <- ggplot(tdf %>% mutate(x=1:n(),c1=factor(c1),c2=factor(c2))) +
  geom_rect(data=csplits,
            aes(xmin=s,xmax=e,color=as.factor(Value),
                fill=as.factor(Value),
                ymin=min(c(tdf$y_basic,tdf$error_t))-0.05*yiqr,ymax=max(c(tdf$y_basic,tdf$error_t))+0.05*yiqr),alpha=0.6) +
  #geom_line(data=y_quants,aes(x=(s+e)/2,y=yu),linetype='dashed') +
  #geom_line(data=y_quants,aes(x=(s+e)/2,y=yl),linetype='dashed') +
  scale_color_manual('Category', values=cetcolor::cet_pal(nc_v1,'r3')) +
  scale_fill_manual('Category', values=cetcolor::cet_pal(nc_v1,'r3')) +

  geom_line(aes(x=x,y=error_t)) +
  geom_line(aes(x=t,y=y_basic),color='blue') +
  geom_errorbar(data=y_quants %>% filter(direction=='up'),aes(x=(s+e)/2,
                                                              ymin=ifelse(direction=='down',ymean,ymean-v1e),
                                                              ymax=ifelse(direction=='up',ymean,ymean-v1e)
  ),
  color='white',
  size=1.5) +
  geom_errorbar(data=y_quants %>% filter(direction=='down'),aes(x=(s+e)/2,
                                                                ymin=ifelse(direction=='down',ymean,ymean-v1e),
                                                                ymax=ifelse(direction=='up',ymean,ymean-v1e)
  ),
  color='#444',
  size=1.5) +
  geom_point(data=y_quants, aes(x=(s+e)/2,y=ymean),color='#33E',size=2.5-sr_adjustment,shape=15) +
  theme_bw() +
  guides(size='none',color='none') +
  xlab('t') +
  ylab('y') +
  ggtitle('Time Dependent Error Across Time-Dependent Categories',
          subtitle='Blue lines: y (square is average for category within contiguous zone)\nWhite & gray error bars: +/- gap caused by categorical effect') +
  geom_abline(slope=0,intercept=0,linetype='dotted')


fpng('example_effect.png')
print(plt_basic2)
dev.off()



plt_basic2_nobar <- ggplot(tdf %>% mutate(x=1:n(),c1=factor(c1),c2=factor(c2))) +
  geom_rect(data=csplits,
            aes(xmin=s,xmax=e,color=as.factor(Value),
                fill=as.factor(Value),
                ymin=min(c(tdf$y_basic,tdf$error_t))-0.05*yiqr,ymax=max(c(tdf$y_basic,tdf$error_t))+0.05*yiqr),alpha=0.6) +
  #geom_line(data=y_quants,aes(x=(s+e)/2,y=yu),linetype='dashed') +
  #geom_line(data=y_quants,aes(x=(s+e)/2,y=yl),linetype='dashed') +
  scale_color_manual('Category', values=cetcolor::cet_pal(nc_v1,'r3')) +
  scale_fill_manual('Category', values=cetcolor::cet_pal(nc_v1,'r3')) +
  geom_line(aes(x=x,y=error_t)) +
  geom_line(aes(x=t,y=y_basic),color='blue') +
  geom_point(data=y_quants, aes(x=(s+e)/2,y=ymean),color='#33E',size=2.5-sr_adjustment,shape=15) +
  theme_bw() +
  guides(size='none',color='none') +
  xlab('t') +
  ylab('y') +
  ggtitle('Time Dependent Error Across Time-Dependent Categories',
          subtitle='Blue lines: y (square is average for category within contiguous zone)\n') +
  geom_abline(slope=0,intercept=0,linetype='dotted')

fpng('example_effect_nobar.png')
print(plt_basic2_nobar)
dev.off()


plt_basic2_nobar_fx <- ggplot(tdf %>% mutate(x=1:n(),c1=factor(c1),c2=factor(c2))) +
  geom_rect(data=csplits,
            aes(xmin=s,xmax=e,color=as.factor(Value),
                fill=as.factor(Value),
                ymin=min(c(tdf$y_basic,tdf$error_t))-0.05*yiqr,ymax=max(c(tdf$y_basic,tdf$error_t))+0.05*yiqr),alpha=0.6) +
  scale_color_manual('Category', values=cetcolor::cet_pal(nc_v1,'r3')) +
  scale_fill_manual('Category', values=cetcolor::cet_pal(nc_v1,'r3')) +
  geom_line(aes(x=x,y=error_t)) +
  geom_line(aes(x=t,y=y_basic),color='blue') +
  geom_point(data=y_quants, aes(x=(s+e)/2,y=ymean),color='#33E',size=2.5-sr_adjustment,shape=15) +
  #geom_line(data=y_quants, aes(x=(s+e)/2,y=ymean),color='#88E',size=2.5-sr_adjustment) +
  geom_line(data=y_quants,aes(x=(s+e)/2,y=v1e),linetype='dashed',color='red',linewidth=1) +
  geom_point(data=y_quants,aes(x=(s+e)/2,y=v1e),color='red',size=3) +
  theme_bw() +
  guides(size='none',color='none') +
  xlab('t') +
  ylab('y') +
  ggtitle('Time Dependent Error Across Time-Dependent Categories',
          subtitle='Blue lines: y (square is average for category within contiguous zone)\nRed dashed line = connecting red dots of actual categorical effect') +
  geom_abline(slope=0,intercept=0,linetype='dotted')


fpng('example_effect_nobar_fx.png')
print(plt_basic2_nobar_fx)
dev.off()


plt_basic2_nobar_fx_v2 <- ggplot(tdf %>% mutate(x=1:n(),c1=factor(c1),c2=factor(c2))) +
  geom_rect(data=csplits,
            aes(xmin=s,xmax=e,color=as.factor(Value),
                fill=as.factor(Value),
                ymin=min(c(tdf$y_basic,tdf$error_t))-0.05*yiqr,ymax=max(c(tdf$y_basic,tdf$error_t))+0.05*yiqr),alpha=0.6) +
  scale_color_manual('Category', values=cetcolor::cet_pal(nc_v1,'r3')) +
  scale_fill_manual('Category', values=cetcolor::cet_pal(nc_v1,'r3')) +
  geom_line(aes(x=x,y=error_t)) +
  #geom_line(aes(x=t,y=y_basic),color='blue') +
  geom_point(data=y_quants, aes(x=(s+e)/2,y=ymean),color='#33E',size=2.5-sr_adjustment,shape=15) +
  geom_line(data=y_quants, aes(x=(s+e)/2,y=ymean),color='#33E',linewidth=1,linetype='dashed') +
  geom_line(data=y_quants,aes(x=(s+e)/2,y=v1e),linetype='dashed',color='red',linewidth=1) +
  geom_point(data=y_quants,aes(x=(s+e)/2,y=v1e),color='red',size=3) +
  theme_bw() +
  guides(size='none',color='none') +
  xlab('t') +
  ylab('y') +
  ggtitle('Time Dependent Error Across Time-Dependent Categories',
          subtitle='Blue dashed lines: connecting y average for category within contiguous zone\nRed dashed line = connecting red dots of actual categorical effect') +
  geom_abline(slope=0,intercept=0,linetype='dotted')


fpng('example_effect_nobar_fx_v2.png')
print(plt_basic2_nobar_fx_v2)
dev.off()

# SINE MODEL

fpng('sine_model_basic_time_error.png')
sine_model_list[[4]]$results$plots$time_error +
  ggtitle('Rectified Sine Model Time-Drift Estimate vs. Actual',subtitle='N=1,010')
dev.off()

fpng('sine_model_basic_param_comp.png')
sine_model_list[[4]]$results$plots$param_comp +
  ggtitle('Parameter estimate comparison of Stan model vs. linear model')
dev.off()





# how long do the models take to run

n_models = length(model_list)
mn = names(model_list)
get_k <- function(x) as.numeric(gsub('.._(...)..','\\1',x))
kvals = get_k(mn)

model_time_df = bind_rows(
  lapply(
    1:n_models,
    function(x){
      data.frame(
        time=model_list[[x]]$time_elapsed,
        N=model_list[[x]]$N,
        k=kvals[x]
      )
    }
  )
)

fpng('run_times.png')
ggplot(model_time_df %>% mutate(N=factor(N))) +
  geom_boxplot(
    aes(x=N,y=time,fill=N)
  ) +
  theme_bw() +
  xlab('Sample Size') +
  ylab('Run Time (minutes)') +
  guides(fill='none') +
  ggtitle(
    "Run Time of Time-Drift Linear Models by Sample Size"
  )
dev.off()


# results
numeric_results_df = bind_rows(
  lapply(model_list,get_numeric_results)
) %>% mutate(
  model_name = names(model_list),
  k = get_k(model_name),
  time_error_type='AR'
)

sine_numeric_results_df = bind_rows(
  lapply(sine_model_list,get_numeric_results)
) %>%
  mutate(
    model_name = names(sine_model_list),
    k = get_k(model_name),
    time_error_type='sine'
  )

# plot sigma and sigma_t boxplots

(plt_sigma <- ggplot(numeric_results_df %>% mutate(N=factor(N))) +
  geom_boxplot(
    aes(x=N,y=sigma,fill=N)
  ) +
  theme_bw() +
  geom_abline(slope=0,intercept=1,linetype='dashed',color='red') +
  guides(fill='none') +
  ggtitle(
    'Estimate of standard error for models of different sample sizes',
    subtitle='"True" value: σ=1'
  ) +
  xlab('Sample Size') +
  ylab('Standard Error Estimate')
)

(plt_sigma_t <- ggplot(numeric_results_df %>% mutate(N=factor(N))) +
  geom_boxplot(
    aes(x=N,y=sigma_t,fill=N)
  ) +
  theme_bw() +
  geom_abline(slope=0,intercept=numeric_results_df$sigma_t_true[1],linetype='dashed',color='red') +
  guides(fill='none') +
  ggtitle(
    'Estimate of time drift step error',
    subtitle=sprintf('"True" value: σ=%0.2f',numeric_results_df$sigma_t_true[1])
  ) +
  xlab('Sample Size') +
  ylab('Time Drift step error Estimate')
)

(plt_phi <- ggplot(numeric_results_df %>% mutate(N=factor(N))) +
  geom_boxplot(
    aes(x=N,y=ar_coef_est,fill=N)
  ) +
  theme_bw() +
  geom_abline(slope=0,intercept=numeric_results_df$ar_coef[1],linetype='dashed',color='red') +
  guides(fill='none') +
  ggtitle(
    'Estimate of time drift decay coefficient φ',
    subtitle=sprintf('"True" value:  φ=%0.4f',numeric_results_df$ar_coef[1])
  ) +
  xlab('Sample Size') +
  ylab('φ Estimate')
)

fpng('parameter_estimates_boxplots.png',height=1080)
grid.arrange(plt_sigma,plt_sigma_t,plt_phi)
dev.off()


# examples of time error, parameter error, and yhat error

fpng('single_result_param_estimate.png')
model_list[['10_020_1']]$results$plots$param_comp +
  ggtitle('Parameter Estiamtes of Linear Model vs. Stan Model',subtitle='N=1,010')
dev.off()

fpng('single_result_time_error_estimate.png')
model_list[['10_020_1']]$results$plots$time_error +
  ggtitle('Time-drift error estimate vs. actual',subtitle='both series re-centered at 0, N=1,010')
dev.off()

fpng('single_result_yhat_estimate.png')
model_list[['10_020_1']]$results$plots$yhat_comp +
  ggtitle('y estimates of linear model vs. Stan model',subtitle='N=1,010')
dev.off()


# condensed parameter comparison
fpng('r2_comparison.png')
ggplot(numeric_results_df %>%
         transmute(N=factor(N),Stan=param_r2_stan,linear=param_r2_linear) %>%
         pivot_longer(c('Stan','linear'))) +
  geom_boxplot(aes(x=N,fill=name,y=value,color=name)) +
  theme_bw() +
  scale_y_continuous(limits=c(0.75,1.01)) +
  scale_fill_manual(
    'Model',
    values=c('Stan'='#F2F','linear'='#2FE')
  ) +
  scale_color_manual(
    'Model',
    values=c('Stan'='#B2B','linear'='#2BA')
  ) +
  ggtitle('Cross-R² of parameter estimates of linear and Stan models',
          subtitle='1-Cross-MSE/var(cross_cat_diff(β))\ntwo points for linear model cut off below 0.5') +
  xlab('Sample Size') +
  ylab('Cross-R²')
dev.off()


fpng('r2_comparison_sine.png')
ggplot(sine_numeric_results_df %>%
         transmute(N=factor(N),Stan=param_r2_stan,linear=param_r2_linear) %>%
         pivot_longer(c('Stan','linear'))) +
  geom_boxplot(aes(x=N,fill=name,y=value,color=name)) +
  theme_bw() +
  scale_fill_manual(
    'Model',
    values=c('Stan'='#F2F','linear'='#2FE')
  ) +
  theme_bw() +
  scale_color_manual(
    'Model',
    values=c('Stan'='#B2B','linear'='#2BA')
  ) +
  ggtitle('Cross-R² of parameter estimates of linear and Stan models with rectified sine wave drift',
          subtitle='1-Cross-MSE/var(cross_cat_diff(β))') +
  xlab('Sample Size') +
  ylab('Cross-R²')
dev.off()


