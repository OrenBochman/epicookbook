library(pomp)
library(reshape2)
library(ggplot2)


sir.step <- "
double lambd = bet * (I+iota) / N;
double ifrac = 1.0 - exp(-lambd *dt);
double rfrac = 1.0 - exp(-gamm*dt);
double infection = rbinom(S, ifrac);
double recovery = rbinom(I, rfrac);
S += -infection;
I += infection - recovery;
R += recovery;
Y += infection;
"

rmeas <- "
Z = Y;
"

pomp(
  data=data.frame(time=seq(0,200,by=0.1),Z=rep(0,2001)),
  times="time",
  t0=0,
  rmeasure=Csnippet(rmeas),
  rprocess=euler.sim(
    step.fun=Csnippet(sir.step),
    delta.t=0.1
  ),
  statenames=c("S","I","R","Y"),
  paramnames=c(
    "bet","gamm","iota", "N"
  ), 
  initializer=function(params, t0, ...) {
    x0 <- c(S=999,I=1,R=0,Y=0)
    x0
  }, 
  params=c(bet=0.1, gamm=0.05, iota=0.01, N=1000.0)
) -> sir

set.seed(42)

sir_sim <- simulate(sir)

sir_sim

sir_out <- data.frame(Time=sir_sim@times,Cases=as.integer(sir_sim@data))

sir_out_long <- melt(as.data.frame(sir_out),"Time")

# Plot
ggplot(sir_out_long,aes(x=Time,y=value,colour=variable,group=variable)) +
  geom_step(lwd=2) + xlab("Time") + ylab("Cases") + theme(legend.position="none")
