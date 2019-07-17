library(GillespieSSA)
library(reshape2)
library(ggplot2)


a <- c("beta*S*I","gamma*I")
nu <- matrix(c(-1,+1,0,0,-1,+1),nrow=3,ncol=2,byrow=FALSE)

# GillespieSSA
parms <- c(beta=0.1/1000,gamma=0.05)
x0 <- c(S=999,I=1,R=0)
tf <- 200

set.seed(42)

sir_out <- ssa(x0,a,nu,parms,tf=tf,simName="SIR")

while(sir_out$stats$nSteps==1){
    sir_out <- ssa(x0,a,nu,parms,tf=tf,simName="SIR")
}

head(sir_out$data)

sir_out_long <- melt(as.data.frame(sir_out$data),"V1")

# Plot
ggplot(sir_out_long,aes(x=V1,y=value,colour=variable,group=variable)) + geom_line(lwd=2) +
  xlab("Time") + ylab("Number")
