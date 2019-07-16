library(reshape2)
library(ggplot2)


sir <- function(beta, gamma, N, S0, I0, R0, tf) {
    time <- 0
    S <- S0
    I <- I0
    R <- R0
    ta <- numeric(0)
    Sa <- numeric(0)
    Ia <- numeric(0)
    Ra <- numeric(0)
    while (time < tf) {
        ta <- c(ta, time)
        Sa <- c(Sa, S)
        Ia <- c(Ia, I)
        Ra <- c(Ra, R)
        pf1 <- beta * S * I
        pf2 <- gamma * I
        pf <- pf1 + pf2
        dt <- rexp(1, rate = pf)
        time <- time + dt
        if (time > tf) {
            break
        }
        ru <- runif(1)
        if (ru < (pf1/pf)) {
            S <- S - 1
            I <- I + 1
        } else {
            I <- I - 1
            R <- R + 1
        }
        if (I == 0) {
            break
        }
    }
    results <- data.frame(time = ta, S = Sa, I = Ia, R = Ra)
    return(results)
}

set.seed(42)

sir_out <- sir(0.1/1000,0.05,1000,999,1,0,200)

if(dim(sir_out)[1]==1){
    sir_out <- sir(0.1/1000,0.05,1000,999,1,0,200)
}

head(sir_out)

sir_out_long <- melt(sir_out,"time")

# Plot
ggplot(sir_out_long,aes(x=time,y=value,colour=variable,group=variable)) +
  geom_line(lwd=2) + xlab("Time") + ylab("Number")
