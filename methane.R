library(pacman)
p_load(tidyverse)

methane_rf <- function(CH4_gton) {
  - ( 1.0569 * ( log(0.8) - log(CH4_gton/2.0) ) +
        1.4 * 0.16019 * ( log(0.8)^2 -  log(CH4_gton/2.0)^2 ) )
}

co2_rf <- function(pCO2) {
  - ( 2.6222 * ( log(280.0) - log(pCO2) ) +
        0.2960 * ( log(280.0)^2 - log(pCO2)^2 ) )
}

methane_dynamics <- function(ch4_gton, pco2, ch4_src, co2_src, age_exp, lifetime_0, ch4_gton_0, time_step) {
  ch4_lifetime = lifetime_0 * (ch4_gton / ch4_gton_0)^age_exp
  ch4_sink = ch4_gton / ch4_lifetime

  ch4_gton = ch4_gton + (ch4_src - ch4_sink) * time_step

  pco2 = pco2 + co2_src * time_step

  ch4_rf = methane_rf(ch4_gton)
  co2_rf = co2_rf(pco2)

  df = tibble(ch4_gton = ch4_gton / 2.0, pco2 = pco2,
              ch4_src = ch4_src, ch4_sink = ch4_sink,
              co2_src = co2_src,
              ch4_forcing = ch4_rf, co2_forcing = co2_rf,
              ch4_lifetime = ch4_lifetime)
  invisible(df)
}

calc_methane <- function(natural_rate = 0.26, anthro_rate = 0.2,
                         spike_gton = 50.0, spike_duration = 20.,
                         spinup_years = 50, pre_spike_years = 50,
                         sim_duration = 100) {
  spinup_years = as.integer(round(spinup_years))
  pre_spike_years = as.integer(round(pre_spike_years))
  spike_duration = as.integer(round(spike_duration))
  sim_duration = as.integer(round(sim_duration))

  start_time = -(spinup_years + pre_spike_years)

  ch4_gton_0 = 1.6
  lifetime_0 = 1.6 / 0.26
  age_exp = 0.26

  ch4_gton = 1.0E-3
  spike_rate = spike_gton / spike_duration

  pco2 = 280.0
  co2_tau = 0.0065

  switch_1_time = start_time + spinup_years
  switch_2_time = 0
  switch_3_time = spike_duration

  timespan = 1 + sim_duration - start_time

  data = tibble()

  row_counter = 0

  for (time in seq(start_time, sim_duration, 1)) {
    row_counter = row_counter + 1

    ch4_src = natural_rate
    if (time > switch_1_time)
      ch4_src = ch4_src + anthro_rate

    if (time > switch_2_time && time <= switch_3_time)
      ch4_src = ch4_src + spike_rate

    if (time > switch_1_time)
      co2_src = pco2 * co2_tau
    else
      co2_src = 0.0

    df = methane_dynamics(ch4_gton, pco2, ch4_src, co2_src, age_exp, lifetime_0, ch4_gton_0, 1) %>%
      mutate(year = time)

    ch4_gton = df$ch4_gton * 2.0
    pco2 = df$pco2

    data = bind_rows(data, df)
    if (time > switch_3_time && abs(df$ch4_src - df$ch4_sink) / df$ch4_gton < 1E-3 & (time %% 50) == 0)
      break
  }

  if (time < sim_duration) {
    for (time in seq(time, sim_duration, 50)) {
      df = methane_dynamics(ch4_gton, pco2, ch4_src, co2_src, age_exp, lifetime_0, ch4_gton_0, 50) %>%
        mutate(year = time)
      data = bind_rows(data, df)
    }
  }

  invisible(data)
}
