source('stan_ar_test_organizer.R')

SAV_DIRECTORY = file.path(
  'stan_ar_trials',
  'trials3'
)
dir.create(SAV_DIRECTORY,showWarnings=FALSE,recursive=TRUE)

model_list = list()
# k of 20, 50, 100
for (nl in c(5, 10, 20, 40)){
  for (k in c(20,50,100)){
    for (i in 1:3){
      key=sprintf('%02d_%03d_%s',nl,k,i)
      print(sprintf('%s:%s:%s',nl,k,i))
      seed = 2024*9*18 + k*i^2*nl^3
      set.seed(seed)

      filename = sprintf('ar_trial_%02d_%03d_%s.RDS',nl,k,i)

      filepath = file.path(SAV_DIRECTORY, filename)

      if (file.exists(filepath)){
        print(sprintf('%s exists, loading & skipping...',filepath))
        model_list[[key]] = readRDS(filepath)
        #model_list[[key]] = add_results(model_list[[key]])
        #saveRDS(model_list[[key]], file=filepath)
        next
      }

      mm = run_model_with_params(list(
        # padding params
        n_legs = nl,
        padding_control = c(20,1/20),
        leg_size_control = c(0,0,2,0.7,0,1),
        leg_size_constant=100,
        # time control
        ar_coef = 0.9995,
        ar_alpha = 120,
        ar_beta = 2.2,
        nc_v1 = k,
        nc_v2 = 5,
        c_v1_randomized='yes',
        nc_v1_splits = 80,
        seed = seed,
        showPlot=FALSE
      ))
      print('adding results...')
      mmr = add_results(mm)
      mmr$seed = seed
      print(sprintf('writing %s to disk', filepath))
      saveRDS(mmr, file=filepath)
      model_list[[key]] = mmr
    }
  }
}

# sinusoidal


sine_model_list = list()
# k of 20, 50, 100
for (nl in c(5, 10, 20)){
  for (k in c(20,100)){
    for (i in 1:1){
      key=sprintf('%02d_%03d_%s',nl,k,i)
      print(sprintf('%s:%s:%s',nl,k,i))
      seed = 2024*9*18 + k*i^2*nl^3
      set.seed(seed)

      filename = sprintf('sine_ar_trial_%02d_%03d_%s.RDS',nl,k,i)

      filepath = file.path(SAV_DIRECTORY, filename)

      if (file.exists(filepath)){
        print(sprintf('%s exists, loading & skipping...',filepath))
        sine_model_list[[key]] = readRDS(filepath)
        next
      }

      mm = run_model_with_params(list(
        # padding params
        n_legs = nl,
        padding_control = c(20,1/20),
        leg_size_control = c(0,0,2,0.7,0,1),
        leg_size_constant=100,
        # time control
        ar_coef = 0.9995,
        ar_alpha = 120,
        ar_beta = 2.2,
        nc_v1 = k,
        nc_v2 = 5,
        c_v1_randomized='yes',
        nc_v1_splits = 80,
        seed = seed,
        showPlot=FALSE,
        time_noise_type='sine',
        sine_params = list(
          n_periods = 1.33,
          phase_shift = 0,
          rectify_magnitude=0.85
        )
      ))
      print('adding results...')
      mmr = add_results(mm)
      mmr$seed = seed
      print(sprintf('writing %s to disk', filepath))
      saveRDS(mmr, file=filepath)
      sine_model_list[[key]] = mmr
    }
  }
}

# "faulty" models


filepath = file.path(SAV_DIRECTORY, 'faulty_sampling_1.RDS')

if (!file.exists(filepath)){
  nl = 10
  k=20
  seed = 12345
  mm_faulty = run_model_with_params(list(
    # padding params
    n_legs = nl,
    padding_control = c(20,1/20),
    leg_size_control = c(0,0,2,0.7,0,1),
    leg_size_constant=100,
    # time control
    ar_coef = 0.9995,
    ar_alpha = 120,
    ar_beta = 2.2,
    nc_v1 = k,
    nc_v2 = 5,
    c_v1_randomized='no',
    nc_v1_splits = 2,
    seed = seed,
    showPlot=TRUE
  ))

  mm_faulty_r = add_results(mm_faulty)
  saveRDS(mm_faulty_r, filepath)
} else {
  mm_faulty_r = readRDS(filepath)
}
