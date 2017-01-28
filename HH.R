# A solver in R for the Hodgkin-Huxley equations in discrete time
#  Sean Caron (scaron@umich.edu)

# Time step, mS
dt <- 0.01
# Membrane capacitance, mF/cm^2
cm <- 1.0
# Potassium conductance, mS/cm^2
G_K_bar <- 36.0
# Sodium conductance, mS/cm^2
G_Na_bar <- 120.0
# Leak conductance, mS/cm^2
G_L_bar <- 0.3
# Sodium potential, mV
V_Na <- -115.0
# Potassium potential, mV
V_K <- 12.0
# Leakage potential, mV
V_L <- 10.613

# Clamp at fixed voltage? 1 = true, 0 = false
s <- readline(prompt="Clamp at fixed voltage (1=yes/0=no)? ")
s <- as.integer(s)

# Input special preconditions? 1 = true, 0 = false
r <- readline(prompt="Special initial conditions (1=yes/0=no)? ")
r <- as.integer(r)

if (r == 1) {
    # Membrane voltage before t = 0, mV, INPUT ME
    voi <- readline(prompt="Membrane voltage before t = 0 (mV)= ")
    voi <- as.numeric(voi)

    vo <- voi

    # Calculate alpha and beta values for m, n and h
    bn <- 0.125 * exp(vo/80.0)
    if (vo == -10) {
        an <- 0.1
    } else {
        an <- (0.01 * (vo + 10.0)) / (exp((vo + 10.0)/10.0)-1.0)
    }

    if (vo == -25) {
        am <- 1.0
    } else {
        am <- (0.1 * (vo + 25.0)) / (exp((vo + 25.0)/10.0)-1.0)
    }

    bm <- 4.0 * exp(vo/18.0)
    ah <- 0.07 * exp(vo/20.0)
    bh <- 1.0 / (exp((vo + 30)/10.0)+1.0)

    # Calculate m, n and h values for this initial membrane polarization
    xn <- an / (an + bn)
    xm <- am / (am + bm)
    xh <- ah / (ah + bh)
} else {
    # No, the membrane is just at RNMP
    # We just take the (known) m, n and h values for RNMP
    xn <- 0.31768
    xm <- 0.05293
    xh <- 0.59612
}

# What is the stimulus voltage or clamp voltage?
if (s == 0) {
    v <- readline(prompt="Stimulus voltage (mV)= ")
    v <- as.numeric(v)
} else {
    v <- readline(prompt="Clamp voltage (mV)= ")
    v <- as.numeric(v)
}

# Print a data point after how many steps?
j <- readline(prompt="print a point after this many steps= ")
j <- as.integer(j)
k <- j

# How many total steps of the simulation do we want to run?
nmax <- readline(prompt="total number of steps to run= ")
nmax <- as.integer(nmax)

# Do the simulation
result <- matrix(nrow=nmax, ncol=4)
for (i in 1:nmax) {
    vo <- v

    # Calculate alphas and betas for the membrane voltage at this step

    bn <- 0.125 * exp(vo/80.0)
    if (vo == -10) {
        an <- 0.1
    } else {
        an <- (0.01 * (vo + 10.0)) / (exp((vo+10.0)/10.0)-1.0)
    }

    if (vo == -25) {
        am <- 1.0
    } else {
        am <- (0.1 * (vo + 25.0)) / (exp((vo + 25.0)/10.0)-1.0)
    }

    bm <- 4.0 * exp(vo/18.0)
    ah <- 0.07 * exp(vo/20.0)
    bh <- 1.0 / (exp((vo + 30)/10.0)+1.0)

    # Calculate differential m, n and h values at this step
    dxn <- dt * (an * (1.0 - xn) - bn*xn)
    dxm <- dt * (am * (1.0 - xm) - bm*xm)
    dxh <- dt * (ah * (1.0 - xh) - bh*xh)

    xn <- xn + dxn
    xm <- xm + dxm
    xh <- xh + dxh

    # If voltage is not fixed (clamped) then update it
    if (s == 0) {
        # Calculate differential membrane voltage at this step
        dv <- (-dt/cm) * (G_K_bar*(xn^4)*(v-V_K) + G_Na_bar*(xm^3)*xh*(v-V_Na) + G_L_bar*(v-V_L))
        # Find total new membrane voltage by adding the differential to the existing
        v <- v + dv
    }

    # Save all points for plotting
    result[i, 1] <- v
    result[i, 2] <- xn
    result[i, 3] <- xm
    result[i, 4] <- xh

    # Additionally, print out a data point for the user at the specified interval
    kj <- j - k
    if (kj == 0) {
        print(v)
    }

    k <- k - 1

    if (k == 0) {
        k <- j
    }
}

steps <- seq(0, nmax*dt, length=nmax)
par(mfrow=c(2,2))

# Plot potential
matplot(steps, result[,1], type = "l", col=4, main="Potential", xlab="Time (ms)", ylab="Membrane voltage (mV)")
grid()

# Plot m, n and h
n <- result[,2]
m <- result[,3]
h <- result[,4]
matplot(steps, cbind(n, m, h), pch=1, col=c(4,2,6), type="l", main="n, m, h", xlab="Time (ms)", ylab="n, m, h")
legend("topright", legend=c("n", "m", "h"), pch=1, col=c(4, 2, 6))
grid()

# Plot conductance
G_K <- G_K_bar * n^4
G_Na <- G_Na_bar * ((m^3) * h)
matplot(steps, cbind(G_Na, G_K), pch=1, col=c(4, 2), type="l", main="Ion Channel Conductivity", xlab="Time (ms)", ylab="Conductivity (mmho/cm^2)")
legend("topright", legend=c("G_K", "G_Na"), pch=1, col=c(4,2))
grid()

# Plot current
I_K <- (result[,1] - V_K) * G_K
I_Na <- (result[,1] - V_Na) * G_Na
matplot(steps, cbind(I_Na, I_K), pch=1, col=c(4,2), type="l", main="Ion Channel Currents", xlab="Time (ms)", ylab="Current (uA/cm^2)")
legend("topright", legend=c("I_K", "I_Na"), pch=1, col=c(4,2))
grid()

