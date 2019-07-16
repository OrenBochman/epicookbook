library(deSolve)
library(ggplot2)
library(reshape2)


sis_ode <- function(times,x,parms){
  with(as.list(c(parms,x)),{
    # ODEs
    SH<-nH-x[[1]] #Calculate SH as nH-IH
    SL<-(1-nH)-x[[2]] #Calculate SL as nL-IL
    dIH <-+(betaHH*IH+betaHL*IL)*SH-gamma*IH
    dIL <-+(betaLH*IH+betaLL*IL)*SL-gamma*IL
    
    der<-c(dIH,dIL)
    list(der)
  })
}

parms <- c(betaHH=10,betaHL=0.1,betaLH=0.1,betaLL=1,gamma=1,nH=0.2)
x <- c(IH=0.00001,IL=0.001)
times <-seq(0,15,1)
sis_out <- as.data.frame(lsoda(x,times,sis_ode,parms))

sis_out_long <- melt(sis_out,"time") #Collapse dataset from "wide" to "long" format for plotting
ggplot(sis_out_long,aes(x=time,y=value,colour=variable,group=variable)) + geom_line(lwd=2) +
  labs(x="Time (Years)",y="Proportion of Population",color="Risk Group")
