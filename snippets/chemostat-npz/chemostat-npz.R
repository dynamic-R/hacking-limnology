npz <- function(time, y, parms) {
  with(as.list(c(y, parms)), {

    p_growth <- r_pgrow * N / (km_n + N) * P

    z_grazing    <- r_zgraz  * P / (km_p + P) * Z
    z_growth     <- asseff * z_grazing
    n_recycling  <- (1 - asseff) * c_pn * z_grazing

    n_import <- d * N0
    n_export <- d * N
    p_export <- d * P
    z_export <- d * Z

    dN_dt   <- n_import - c_pn * p_growth - n_export + n_recycling
    dP_dt   <- p_growth - z_grazing - p_export
    dZ_dt   <- z_growth - z_export

    list(c(dN_dt, dP_dt, dZ_dt))
  })
}
parms <- c(
  r_pgrow  = 0.5,   # phytoplankton growth parameter 1/d
  r_zgraz  = 0.4,   # zooplankton grazing parameter 1/d

  km_n     = 0.5,   # Monod constant of phyto growth on nutrient (mmol/m3)
  km_p     = 100,   # Monod constant of zoo growth on phyto (mmol/m3)
  c_pn     = 1/106, #stoichiometric conversion from phosporus P to phyto C (P:C ratio)
  asseff   = 0.3,   # zooplankton assimilation efficiency (-)

  d       = 0.1,    # dilution rate 1/d
  N0      = 5       # P in inflow (mmol/m3)
)

times <- seq(0, 40, 0.1)        # (d)

# Zoo and Phyto as C and Phosphorus P (mol/m3)
y  <- c(N = 5, P=1, Z = 10)

parms["d"] <- 0.1
out0 <- ode(init, times, npz, parms)
plot(out0, mfrow = c(1, 3))
