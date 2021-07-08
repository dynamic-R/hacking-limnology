library("deSolve")
library("rootSolve")

chemostat <- function(time, init, parms) {
  with(as.list(c(init, parms)), {
    r       <- r_max * P / (kp + P)
    r_z     <- r_max_z * Alg / (k_alg + Alg)

    dZ_dt   <- r_z * Z - D * Z
    dAlg_dt <- r * Alg - D * Alg - r_z * Z * 1/Y_z
    dP_dt   <- D * (P0 - P) - r * Alg * 1/Y
    list(c(dZ_dt, dAlg_dt, dP_dt))
   })
}
parms <- c(
  r_max = 0.5,    # 1/d
  kp    = 0.5,    # half saturation constant, P (mol/m3)
  Y     = 106,    # yield coefficient (stoichiometric C:P ratio)
  r_max_z = 0.2,  # max growth rate of zooplankton
  k_alg = 100,    # half sat zoopl. - algae
  Y_z   = 0.1,    # yield zooplankton
  D     = 0.1,    # 1/d
  P0    = 5       # P in inflow (mol/m3)
)

times <- seq(0, 200, 0.1)        # (d)

# Zoo and Phyto as C and Phosphorus P (mol/m3)
init  <- c(Z=10, Alg = 10, P = 5)

parms["D"] <- 0.0
out0 <- ode(init, times, chemostat, parms)

parms["D"] <- 0.1
out1 <- ode(init, times, chemostat, parms)

parms["D"] <- 0.2
out2 <- ode(init, times, chemostat, parms)

parms["D"] <- 0.5
out4 <- ode(init, times, chemostat, parms)

plot(out0, out1, out2, out4, mfrow=c(1,3))

times <- seq(0, 2000, 0.1)        # (d)
parms["D"] <- 0.16
out3 <- ode(init, times, chemostat, parms)
plot(out4)

# scenario 0: equilibrium, algae extinct
# scenario 1: Lotka-Volterra cycle
# scenario 2: equilibrium, zooplankton extinct
# scenario 3: equilibrium, coexistence
# scenario 4: equilibrium, both extinct

# scneario X: damped Lotka-Volterra cycle
times <- seq(0, 2000, 0.1)        # (d)
parms["D"] <- 0.12
outx <- ode(init, times, chemostat, parms)
plot(outx)

