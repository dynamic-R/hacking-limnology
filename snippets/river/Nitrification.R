library("deSolve")

parms <- c(
  r_nitri = 2.0,   # 1/day   nitrification rate
  kO2     = 0.001, # mol/m3  half sat. constant
  k2      = 2.8,   # 1/day  re-aeration rate
  O2sat   = 0.3    # mol/m3 saturation concentration
)

Nitrification <- function(t, y, parms) {
  with(as.list(c(y, parms)), {

    aeration      <- k2 * (O2sat - O2)
    nitrification <- r_nitri * O2/(O2 + kO2) * NH4

    dNH4 <- -nitrification
    dNO3 <- +nitrification
    dO2  <-  aeration - 2 * nitrification
    list(c(dNH4, dNO3, dO2))
  })
}

y0 <- c(NH4 = 0.5, NO3 = 0, O2 = 0.3) # in mol/m3

times <- seq(from = 0, to = 5, length.out=101)

out <- ode(y=y0, times = times, func=Nitrification, parms = parms)
plot(out)
