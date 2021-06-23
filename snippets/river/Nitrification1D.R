################################################################################
##
## == Simple Nitrification Model demonstrating the Peterson Matrix approach ==
##
## (based on an exercise of Soetaert and Meysman, Ecological Modelling
##  Nitrification in the Schelde river, 1-D)
################################################################################

library("ReacTran")
library("deSolve")
library("gifski")
#=============================================================================
# Units used in the program:
#=============================================================================

# Space  = m
# Time   = day
# Mass   = mg
# Volume = L

N    <- 100
Grid <- setup.grid.1D(N = N, L = 100000)

#=============================================================================
# Model parameters:
#=============================================================================


## river characteristics
v  <- 1000        # m/day
w  <- 15          # m
h0 <- 0.35        # m
A  <- w * h0      # m^2


Temp <- 5 # deg C
k2 <- 0.12 * 1.07^(Temp - 20) # empirical formula, re-aereation coeff. (1/h)

parameters <- c(rMax          = 0.1,   # 1/day
                kO2           = 5,     # mg/L, half sat. constant
                NH4up         = 5,     # mg/L
                NO3up         = 3.1,   # mg/L
                O2up          = 10,    # mg/L
                NH4dwn        = 0.17,  # mg/L
                NO3dwn        = 1.85,  # mg/L
                O2dwn         = 10,    # mg /L
                D             = 10,    # m2/day, dispersion coefficient
                k2            = k2,    # mg/L
                O2sat         = 10     # mg/L
)

#------------------------------------------------------------------------------
# Model definition:
#------------------------------------------------------------------------------

## Peterson Stoichiometry matrix
stoich <- matrix(c(
  # NH4  NO3  O2
    0,   0,   1,    # reaeration
   -1,  +1,  -4.57  # nitrification
  ),  nrow = 2, byrow = TRUE
)

Nitrification <- function(t, y, parameters) {
  with(as.list(c(parameters)),{
    NH4 <- y[1:N]
    NO3 <- y[(N+1):(2*N)]
    O2  <- y[(2*N+1):(3*N)]

    tran <- cbind(
      tran.1D(C = NH4, D = D, v = v, C.up = NH4up, C.down = NH4dwn, A = A, dx = Grid)$dC,
      tran.1D(C = NO3, D = D, v = v, C.up = NO3up, C.down = NO3dwn, A = A, dx = Grid)$dC,
      tran.1D(C = O2 , D = D, v = v, C.up = O2up,  C.down = O2dwn,  A = A, dx = Grid)$dC
    )

    proc <- cbind(
      k2 * (O2sat - O2),             # re-aeration
      rMax * O2/(O2 + kO2) * NH4     # nitrification
    )

    dY   <- tran + proc %*% stoich

    list(c(dY))#, totalN = NO3 + NH4,  reaer = proc[1,], nitrif = proc[2,])
  })
}

#=============================================================================
# Steady-state solution
#=============================================================================

std <- steady.1D(y = runif(3*Grid$N), names=c("NH4","NO3","O2"),
   func = Nitrification, parms = parameters, nspec = 3, pos=TRUE)

plot(std, grid=Grid$x.mid, type="l")

#=============================================================================
# Dynamic solution
#=============================================================================

# State variables and initial conditions:
NH4_0  <- 0.1   # mg/L
NO3_0  <- 0.1   # mg/L
O2_0   <- 10    # mg/L

state <- c(rep(NH4_0, N), rep(NO3_0, N), rep(O2_0, N))

# Time sequence
time_seq <- seq(from = 0, to = 150, by = 1)

out1D  <- ode.1D(y=state, times=time_seq, func=Nitrification, parms=parameters, nspec=3)
image(out1D, xlab = "time, days", ylab = "Distance, m", grid= Grid$x.mid,
      main =c("NH4","NO3","O2"), add.contour=TRUE)


# for (tt in time_seq) {
#   png(paste0("river", 1000 + tt, ".png"), width=600, height=500)
#   par(mar=c(5,5,3,1), las=1, cex.axis=1.5, cex.main=2)
#   plot.1D(out1D,  ask=FALSE, mfrow=c(3,1), delay=10, type="l", lwd=2, cex.lab=2,
#           ylim=list(c(0, 6), c(0, 10), c(0, 10)),
#           xlab=c("", "", "River km"), ylab=c("NH4 (mg/L)", "NO3", "O2"),
#           subset = (time == tt))
#   dev.off()
# }
#

plot_poly <- function(data, time) {
  poly <- function(x, y, pcol, ...) {
    plot(x, y, type="n", las=1, xlab="", ...)
    polygon(x, y, col=pcol, lty="blank")
    mtext(side=1, line=1.5, text="time")
  }
  y <- data[time,-1]
  NH4 <- y[1:N]
  NO3 <- y[(N+1):(2*N)]
  O2  <- y[(2*N+1):(3*N)]
  x <- c(1, 1:N, N)

  poly(x, c(0, NH4, 0), ylim=c(0, 6), pcol="#a6cee3", ylab="NH4 (mg/L)")
  poly(x, c(0, NO3, 0), ylim=c(0, 10), pcol="#b2df8a", ylab="NO3 (mg/L)")
  poly(x, c(0, O2, 0), ylim=c(0, 10), pcol="#1f78b4", ylab="O2 (mg/L)")
}



for (tt in time_seq)  {
  png(paste0("river", 1000 + tt, ".png"), width=1400, height=800, pointsize = 24)
  par(mfrow=c(3,1))
  par(mar=c(3,5,1,0), las=1, cex.axis=1.4, cex.lab=1.4, cex.main=2)
  plot_poly(out1D, tt)
  dev.off()
}
gifski(dir(pattern="^river.*png$"), gif_file = "river.gif", delay=1/25)

