library(deSolve)
library(reshape)
library(ggplot2)
library(plotly)


micro_1 <- function(times,init,parms){
  with(as.list(c(parms,init)), {
  # ODEs
  dS <- r*(1-S/K)*S - beta*S*I
  dI <- beta*S*I-(mu + alpha)*I
  list(c(dS,dI))
  })
}

w <- 1
m <- 10
beta <- 0.0247*m*w**0.44
r <- 0.6*w**-0.27
mu <- 0.4*w**-0.26
K <- 16.2*w**-0.7
alpha <- (m-1)*mu

parms <- c(beta=beta,r=r,mu=mu,K=K,alpha=alpha)
init <- c(S=K,I=1.)
times = seq(0,10,length.out=101)

sir_out <- lsoda(init,times,micro_1,parms)

sir_out_long <- melt(as.data.frame(sir_out),"time")

# Plot
p <- plot_ly(data = sir_out_long, x = ~time[variable=="S"], y = ~value[variable=="S"],type = "scatter", mode = "lines",name = "S(t)") %>%
add_trace(y = ~value[variable=="I"],type = "scatter", mode = "lines",name = "I(t)") %>%
layout(xaxis = list(title="Time"), yaxis = list(title="Number"))
p

m <- c(5,10,20,40)
ws <- 10^seq(-3,3,length.out=601)
betas <- data.frame("ws"=ws, "m5"=0.0247*m[1]*ws^0.44,"m10"=0.0247*m[2]*ws^0.44,"m20"=0.0247*m[3]*ws^0.44,"m40"=0.0247*m[4]*ws^0.44)

# Plot
p <- plot_ly(betas,x=~ws,y=~m5,type = "scatter", mode = "lines",name = "m = 5") %>%
add_trace(y = ~m10,type = "scatter", mode = "lines",name = "m = 10") %>%
add_trace(y = ~m20,type = "scatter", mode = "lines",name = "m = 20") %>%
add_trace(y = ~m40,type = "scatter", mode = "lines",name = "m = 40") %>%
layout(xaxis = list(title="Weight",type = "log"), yaxis = list(title="beta_min",type = "log"))
p
