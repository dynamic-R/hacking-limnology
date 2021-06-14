library("deSolve")
library("ReacTran")  # lateral transport
library("gifski")    # for the animation
set.seed(123)

## helper to extract 2D state matrix from 1D vector
statematrix <- function(y, N, i) {
  matrix(y[(i-1)*(N^2) + 1:N^2], N, N)
}

## The model equations
SIR2D <- function (t, y, parms)  {
  S  <- statematrix(y, N, 1)
  I  <- statematrix(y, N, 2)
  R  <- statematrix(y, N, 3)

  infect <- beta * I * S
  recovr <- nu * I

  dS <- -infect
  dI <- infect - recovr  + tran.2D(I, dx = dx, dy = dy, D.x = D, D.y = D)$dC
  dR <- recovr
  list(c(dS, dI, dR))
}

## grid
N     <- 81
zero <- rep(0, N)
ndx <- 1:N^2
dx <- 10/N
dy <- 10/N

## model parameters
D     <- 1e-3
S0    <- 100
I0    <- 1
R0    <- 0.0
beta  <- 0.001
nu    <- 0.02


# initial conditions
S  <- matrix(nrow = N, ncol = N, data = S0)
I  <- matrix(nrow = N, ncol = N, data = 0)
R  <- matrix(nrow = N, ncol = N, data = R0)


I[sample(1:length(I), 10)] <- I0

y <- c(S = S, I = I, R = R)

times <- seq(0, 400, 2)
  out <- ode.2D (y = y, func = SIR2D, t = times, nspec=3, parms = NULL,
                dim = c(N, N), method = "adams")


for (tt in times) {
  png(paste0("sir2d", 1000 + tt, ".png"), width=500, height=500)
  par(cex=2, cex.lab=2, cex.axis=2, cex.main=2.5, las=1)
  par(mar=c(4,5,4,6) + 0.1)
  image(out, zlim=c(0, 60), select=2, subset=(time==tt), axes=FALSE, legend = TRUE, ask=FALSE,
        main="Infected", xlab=tt, ylab="", mfrow=c(1,1))
  dev.off()
}

gifski(dir(pattern="^sir.*png$"), gif_file="sir2d.gif", delay=0.1)
