################################################################################
##
## == Simple Nitrification Model demonstrating the Peterson Matrix approach ==
##
## (based on an exercise of Soetaert and Meysman, Ecological Modelling
##  Nitrification in the Schelde river, 1-D)
## Changes (ThPe): river geometry, parameters (unrealistic, for demo only)
################################################################################

library("ReacTran")
library("deSolve")
library("rootSolve")
library("gifski")
#=============================================================================
# Units used in the program
#=============================================================================

# Space  = m
# Time   = day
# Mass   = mol
# Volume = m3

N    <- 100  # number of grid cells
cell <- 1000 # 1000m = 1km / grid cell
Grid <- setup.grid.1D(N = N, L = N * cell)

#=============================================================================
# Model parameters
#=============================================================================


## river characteristics
v  <- 40000       # m/day
w  <- 15          # m
h0 <- 0.35        # m
A  <- w * h0      # m2

parameters <- c(
  r_nitri       = 2.0,   # 1/day   nitrification rate
  kO2           = 0.001, # mol/m3  half sat. constant
  NH4up         = 0.5,   # mol/m3  upper boundary
  NO3up         = 0.0,   # mol/m3  -- " --
  O2up          = 0.3,   # mol/m3  -- " --
  NH4dwn        = 0.0,   # mol/m3  lower boundary
  NO3dwn        = 0.0,   # mol/m3  -- " --
  O2dwn         = 0.3,   # mol/m3  -- " --
  D             = 100,   # m2/day dispersion coefficient
  k2            = 2.8,   # 1/day  re-aeration rate
  O2sat         = 0.3    # mol/m3 saturation concentration
)

#------------------------------------------------------------------------------
# Model definition:
#------------------------------------------------------------------------------

Nitrification <- function(t, y, parameters) {
  with(as.list(c(parameters)),{
    NH4 <- y[1:N]
    NO3 <- y[(N+1):(2*N)]
    O2  <- y[(2*N+1):(3*N)]

    aeration      <- k2 * (O2sat - O2)
    nitrification <- r_nitri * O2/(O2 + kO2) * NH4

    dNH4 <- -nitrification +
      tran.1D(C = NH4, D = D, v = v, C.up = NH4up, C.down = NH4dwn, A = A, dx = Grid)$dC

    dNO3 <- +nitrification +
      tran.1D(C = NO3, D = D, v = v, C.up = NO3up, C.down = NO3dwn, A = A, dx = Grid)$dC

    dO2  <-  aeration - 2 * nitrification +
      tran.1D(C = O2 , D = D, v = v, C.up = O2up,  C.down = O2dwn,  A = A, dx = Grid)$dC

    list(c(dNH4, dNO3, dO2))
  })
}

#=============================================================================
# Steady-state solution
#=============================================================================

std <- steady.1D(y = runif(3*Grid$N), names=c("NH4", "NO3", "O2"),
   func = Nitrification, parms = parameters, nspec = 3, pos = TRUE)

plot(std, mfrow = c(3, 1), las=1, xlab="km", ylab="mol/m3")


#=============================================================================
# Dynamic solution
#=============================================================================

# State variables and initial conditions
NH4_0  <- 0.0 # mol/m3
NO3_0  <- 0.0 # mol/m3
O2_0   <- 0.3 # mol/m3

state <- c(rep(NH4_0, N), rep(NO3_0, N), rep(O2_0, N))

# Time sequence
time_seq <- seq(from = 0, to = 5, length.out=101)

out1D  <- ode.1D(y=state, times=time_seq, func=Nitrification, parms=parameters, nspec=3)
image(out1D, xlab = "time, days", ylab = "Distance, m", grid= Grid$x.mid,
      main =c("NH4","NO3","O2"), add.contour=TRUE, mfrow = c(1, 3))


plot_poly <- function(data, time) {
  poly <- function(x, y, pcol, ylab, ylim) {
    plot(x, y, type="n", las=1, xlab="", ylab="", ylim=ylim)
    polygon(x, y, col=pcol, lty="blank")
    mtext(side=1, line=1.5, text="km")
    mtext(side=2, line=3.5, text=ylab, las=3)

  }
  y <- data[time,-1]
  NH4 <- y[1:N]
  NO3 <- y[(N+1):(2*N)]
  O2  <- y[(2*N+1):(3*N)]
  x   <- c(1, 1:N, N)

  poly(x, c(0, NH4, 0), ylim=c(0, 0.5), pcol="#a6cee3", ylab="NH4 (mol/m3)")
  poly(x, c(0, NO3, 0), ylim=c(0, 0.5), pcol="#b2df8a", ylab="NO3 (mol/m3)")
  poly(x, c(0,  O2, 0), ylim=c(0, 0.3), pcol="#1f78b4", ylab="O2 (mol/m3)")
}

for (tt in 1:length(time_seq))  {
  png(paste0("river", 1000 + tt, ".png"), width=1600, height=800, pointsize = 18)
  par(mfrow=c(3,1))
  par(mar=c(3,5,.5,0), las=1, cex.axis=1.4, cex.lab=1.4, cex.main=2)
  plot_poly(out1D, tt)
  dev.off()
}
gifski(dir(pattern="^river.*png$"), gif_file = "river.gif", width=1600, height=800, delay=1/25)

