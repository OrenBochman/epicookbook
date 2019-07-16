library(odin)
library(reshape2)
library(ggplot2)


sir_ode <- odin::odin({
  ## Derivatives
  deriv(S) <- -b*S*I
  deriv(I) <- b*S*I-g*I
  deriv(R) <- g*I

  ## Initial conditions
  initial(S) <- 0.99
  initial(I) <- 0.01
  initial(R) <- 0.00

  ## parameters
  b <- 0.1
  g <- 0.05
  config(base) <- "sir"
}, ".")

sir_mod <- sir_ode()
times <- seq(0,200,length.out=2001)
sir_out <- sir_mod$run(times)

sir_out_long <- melt(as.data.frame(sir_out),"t")

# Plot
ggplot(sir_out_long, aes(x=t, y=value, colour=variable, group=variable)) +
  geom_line(lwd=2) + xlab("Time")+ylab("Number")
