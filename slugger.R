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

carbon_dynamics <- function(year, ch4_gton, pco2, ch4_src, co2_src, age_exp,
                            lifetime_0, ch4_gton_0, pco2_0, heating_coef,
                            t_prime,
                            time_step) {
  drawdown = tibble(
    time_const = c(300, 5000, 4E5),
    fraction = c(0.75, 0.15, 0.10)
  )

  ch4_sink_baseline = ch4_gton_0 / lifetime_0

  ch4_lifetime = lifetime_0 * (ch4_gton / ch4_gton_0)^age_exp
  ch4_sink = ch4_gton / ch4_lifetime

  delta_ch4 = (ch4_src - ch4_sink) * time_step
  if (time_step > 1 && abs(delta_ch4) >= 0.2 * ch4_gton) {
    message("*** Year ", year, ": Instability detected: Adjusting... ***")
    new_time_step = time_step
    old_ch4_gton = ch4_gton
    unstable = TRUE
    while(new_time_step >= 1 && unstable) {
      unstable = FALSE
      ch4_gton = old_ch4_gton
      max_delta_ch4 = 0
      if (new_time_step <= 1) {
        stop("*** Unstable with one-year time-step ***")
      } else if (new_time_step <= 10) {
        message("*** reverting to one-year time-steps ***")
        new_time_step = 1
      } else {
        new_time_step = new_time_step / 10
        message("*** Unstable dynamics: reducing time-steps to ", new_time_step, " ***")
      }
      for (i in seq(1, time_step, new_time_step)) {
        ch4_lifetime = lifetime_0 * (ch4_gton / ch4_gton_0)^age_exp
        ch4_sink = ch4_gton / ch4_lifetime

        delta_ch4 = (ch4_src - ch4_sink) * new_time_step
        max_delta_ch4 = max(abs(delta_ch4), max_delta_ch4)
        ch4_gton = ch4_gton + delta_ch4
        if (max_delta_ch4 >= 0.2 * old_ch4_gton)
          unstable = TRUE
          message("*** Unstable time step. Reducing. ***")
          break
      }
    }
  } else {
    ch4_gton = ch4_gton + delta_ch4
  }

  if (ch4_gton <= 0) {
    message("*** ch4 = ", ch4_gton, ", lifetime = ", ch4_lifetime, ", old = ", old_ch4_gton, ", src = ", ch4_src, ", sink = ", ch4_sink, ", tstep = ", time_step, " ***")
  }

  if (year > 0) {
    co2_excess = pco2 - pco2_0
    rate = drawdown %>%
      mutate(rate = fraction / time_const * exp(-year / time_const)) %>%
      summarize(rate = sum(rate)) %>%
      simplify()
    co2_sink = rate * co2_excess
    pco2 = pco2 + (co2_src - co2_sink + ch4_sink - ch4_sink_baseline) * time_step
  } else {
    co2_sink = 0.0
  }


  ch4_rf = methane_rf(ch4_gton)
  co2_rf = co2_rf(pco2)

  t_equil = c("co2" = co2_rf, "ch4" = ch4_rf) * heating_coef
  if (year > 0) {
    relax_frac = c(surf = 1. - exp(-time_step / 10.), deep = 1. - exp(-time_step / 1000.))
    t_prime[,"surf"] = t_prime[,"surf"] + (t_equil - t_prime[,"surf"]) * relax_frac["surf"]
    surf_deep = (t_prime[,"deep"] - t_prime[,"surf"]) * relax_frac["deep"]
    t_prime[,"surf"] = t_prime[,"surf"] + surf_deep * 40.
    t_prime[,"deep"] = t_prime[,"deep"] - surf_deep
  }

  df = tibble(year = year,
              ch4_gton = ch4_gton / 2.0, pco2 = pco2,
              ch4_src = ch4_src, ch4_sink = ch4_sink,
              co2_src = co2_src, co2_sink = co2_sink,
              ch4_forcing = ch4_rf, co2_forcing = co2_rf,
              ch4_lifetime = ch4_lifetime,
              warm_atm_co2 = t_prime["co2","surf"],
              warm_atm_ch4 = t_prime["ch4","surf"],
              warm_ocn_co2 = t_prime["co2","deep"],
              warm_ocn_ch4 = t_prime["ch4","deep"]
  )

  invisible(list(t_prime = t_prime, df = df))
}

calc_slugs <- function(ch4_spike_gton = 50.0, ch4_spike_duration = 20.,
                       co2_spike_gton = 50.0, co2_spike_duration = 20.,
                       climate_sensitivity = 3.0,
                       spinup_years = 100, sim_duration = 1E4) {
  spinup_years = as.integer(round(spinup_years))
  ch4_spike_duration = as.integer(round(ch4_spike_duration))
  co2_spike_duration = as.integer(round(co2_spike_duration))
  sim_duration = as.integer(round(sim_duration))

  gton_per_ppm = 2
  opacity_co2 = 4
  delta_t_per_watt = climate_sensitivity / opacity_co2

  natural_rate = 0.26
  ch4_gton_0 = 1.6
  lifetime_0 = ch4_gton_0 / natural_rate
  age_exp = 0.26

  ch4_gton = 1.0E-3
  ch4_spike_rate = ch4_spike_gton / ch4_spike_duration

  pco2_0 = 280.0
  pco2 = pco2_0
  co2_spike_rate = co2_spike_gton / co2_spike_duration

  start_time = -spinup_years

  timespan = 1 + sim_duration - start_time

  data = tibble()

  t_prime = matrix(c(0,0,0,0), nrow = 2,
                   dimnames = list(gas = c("co2", "ch4"), layer = c("surf", "deep")))

  for (time in seq(start_time, sim_duration, 1)) {
    ch4_src = natural_rate

    if (time > 0 && time <= ch4_spike_duration)
      ch4_src = ch4_src + ch4_spike_rate

    if (time > 0 && time <= co2_spike_duration)
      co2_src = co2_spike_rate
    else
      co2_src = 0.0

    res = carbon_dynamics(time, ch4_gton, pco2, ch4_src, co2_src, age_exp, lifetime_0, ch4_gton_0, pco2_0, delta_t_per_watt, t_prime, 1)
    t_prime = res$t_prime
    df = res$df

    ch4_gton = df$ch4_gton * 2.0
    pco2 = df$pco2

    data = bind_rows(data, df)

    if (time > max(co2_spike_duration, ch4_spike_duration) && abs(df$ch4_src - df$ch4_sink) / df$ch4_gton < 2E-3 & (time %% 50) == 0)
      break
  }

  dt = 1
  while(time < sim_duration) {
    dt = dt * 10
    message("*** Adjusting time step to ", dt, " at ", time, " years. ***")
    for (time in seq(time + dt, sim_duration, dt)) {
      ch4_src = natural_rate
      co2_src = 0

      # message(time, " src = ", df$ch4_src, ", sink = ", df$ch4_sink, ", new src = ", ch4_src, ", delta = ", ch4_src - df$ch4_sink, ", ch4 = ", ch4_gton)
      res = carbon_dynamics(time, ch4_gton, pco2, ch4_src, co2_src, age_exp, lifetime_0, ch4_gton_0, pco2_0, delta_t_per_watt, t_prime, dt)
      t_prime = res$t_prime
      df = res$df

      ch4_gton = df$ch4_gton * 2.0
      pco2 = df$pco2

      data = bind_rows(data, df)

      if (abs(df$co2_src - df$co2_sink) / df$ch4_gton < (2E-3 / dt) & (time %% (dt * 10)) == 0)
        break
    }
  }

  invisible(data)
}
