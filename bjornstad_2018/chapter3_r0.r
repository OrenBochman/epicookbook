library(epimdr)
library(bbmle)
library(statnet)


data(niamey)
head(niamey[, 1:5])

par(mar = c(5,5,2,5))
plot(niamey$absweek, niamey$tot_cases, type = "b", 
     xlab = "Week", ylab = "Incidence")
par(new = TRUE)
plot(niamey$absweek, niamey$cum_cases, type = "l", 
     col = "red", axes = FALSE, xlab = NA, ylab = NA, log = "y")
axis(side = 4)
mtext(side = 4, line = 4, "Cumulative incidence")
legend("topleft", legend = c("Cases", "Cumulative"),
       lty = c(1,1), pch = c(1,NA), col = c("black", "red"))

fit=lm(log(cum_cases)~absweek, subset=absweek<7, 
     data=niamey)
r=fit$coef["absweek"]
V=c(1.5, 1.8)
V*r+1

V = c(1.5, 1.8)
f = (5/7)/V
V * r + 1 + f * (1 - f) * (V * r)^2

llik.cb = function(S0, beta, I) {
    n = length(I)
    S = floor(S0 - cumsum(I[-n]))
    p = 1 - exp(-beta * (I[-n])/S0)
    L = -sum(dbinom(I[-1], S, p, log = TRUE))
    return(L)
}

twoweek = rep(1:15, each = 2)
y = sapply(split(niamey$cases_1[1:30], twoweek), sum)
sum(y)

S0cand = 6500       
betacand = seq(0,10, by = .1)
ll = rep(NA, length(betacand))
for(i in 1:length(betacand)){
     ll[i] = llik.cb(S0 = S0cand, beta = betacand[i],
     I = y)}
plot(ll ~ betacand, ylab = "Neg log-lik", 
     xlab = expression(beta)) 
betacand[which.min(ll)]

betacand = 2.3
S0cand = seq(5920,8000, length = 101)
ll = rep(NA, length = 101)
for(i in 1:101){
     ll[i] = llik.cb(S0 = S0cand[i], beta = betacand, 
     I = y)}
plot(ll ~ S0cand, ylab = "Neg log-lik", 
     xlab = expression(S[0]))
S0cand[which.min(ll)]

fit = mle2(llik.cb, start = list(S0 = 7085, beta = 2.3), 
    method = "Nelder-Mead",data = list(I = y))
summary(fit)
confint(fit)

cov2cor(vcov(fit))

sim.cb=function(S0, beta, I0){
   I=I0
   S=S0
   i=1
   while(!any(I==0)){
        i=i+1
        I[i]=rbinom(1, size=S[i-1], prob=1-
           exp(-beta*I[i-1]/S0))
        S[i]=S[i-1]-I[i]
   }
   out=data.frame(S=S, I=I)
   return(out)
}

plot(y, type="n", xlim=c(1,18), 
   ylab="Predicted/observed", xlab="Week")
for(i in 1:100){
     sim=sim.cb(S0=floor(coef(fit)["S0"]), 
     beta=coef(fit)["beta"], I0=11)
     lines(sim$I, col=grey(.5))
}
points(y, type="b", col=2)

data(flu)
plot(flu$day, flu$cases, type="b", xlab="Day", 
     ylab="In bed", log="y")
tail(flu)

fit=lm(log(cases)~day, subset=day<=5, 
     data=flu)
r=fit$coef["day"]
V=c(2,3)
V*r+1

data(ebola)
par(mar = c(5,5,2,5))
plot(ebola$day, ebola$cases, type="b", xlab="Week", 
     ylab="Incidence")
par(new=T)
plot(ebola$day, ebola$cum_cases, type="l", col="red",
      axes=FALSE, xlab=NA, ylab=NA, log="y")
axis(side = 4)
mtext(side = 4, line = 4, "Cumulative incidence")
legend("right", legend=c("Cases", "Cumulative"),
     lty=c(1,1), pch=c(1,NA), col=c("black", "red"))
tail(ebola)

fit=lm(log(cum_cases)~day, subset=day<100, 
   data=ebola)
r=fit$coef["day"]
V=15
f=.5
V*r+1+f*(1-f)*(V*r)^2

#Data aggregation
cases=sapply(split(ebola$cases, 
   floor((ebola$day-.1)/14)), sum)
sum(cases)

#Removal MLE
fit = mle2(llik.cb, start = list(S0 = 20000, beta = 2), 
   method = "Nelder-Mead",data = list(I = cases))
summary(fit)

confint(fit, std.err = c(100, 0.1))

names(ferrari)
ferrari$Ebolacases95
sum(ferrari$Ebolacases95, na.rm = TRUE)
y = c(na.omit(ferrari$Ebolacases95))

fit = mle2(llik.cb, method = "Nelder-Mead", 
    start = list(S0 = 300, beta = 2), 
    data = list(I = y))
fit
confint(fit, std.err = 2)

data(gonnet)
nwt = network(gonnet, directed = TRUE)
plot(nwt, vertex.col = c(0, rep(1, 17), rep(2, 71)))

F1 = quote(beta * S * I/N)
F2 = 0

Vm1 = quote(mu * E + sigma * E)
Vm2 = quote(mu * I + alpha * I + gamma * I)

Vp1 = 0
Vp2 = quote(sigma * E)

V1 = substitute(a - b, list(a = Vm1, b = Vp1))
V2 = substitute(a - b, list(a = Vm2, b = Vp2))

f11 = D(F1, "E"); f12 = D(F1, "I")
f21 = D(F2, "E"); f22 = D(F2, "I")

v11 = D(V1, "E"); v12 = D(V1, "I")
v21 = D(V2, "E"); v22 = D(V2, "I")

paras = list(S = 1, E = 0, I = 0, R = 0, mu = 0, 
   alpha = 0, beta = 5, gamma = .8, sigma = 1.2, N = 1)
f = with(paras,
   matrix(c(eval(f11),eval(f12),eval(f21),
   eval(f22)), nrow = 2, byrow = TRUE))
v=with(paras,
   matrix(c(eval(v11),eval(v12),eval(v21),
   eval(v22)), nrow=2, byrow=TRUE))

max(eigen(f %*% solve(v))$values)

with(paras,
sigma/(sigma+mu)*beta/(gamma+mu+alpha))

F1 = expression(betai * S * I / N + betah* S * H / N +
   betaf * S * F / N)
F2=0
F3=0
F4=0

Vm1 = quote(sigma * E)
Vm2 = quote(Theta * gammah * I + (1 - Theta) * (1-
   Lambda) * gammar * I + (1 - Theta) * Lambda * 
   gammaf * I)
Vm3 = quote(Lambda * etaf * H + (1 - Lambda) * etar * H)
Vm4 = quote(chi * F)

Vp1 = 0
Vp2 = quote(sigma * E)
Vp3 = quote(Theta * gammah * I)
Vp4 = quote((1 - Theta) * (1 - Lambda) * gammar * I+ 
   Lambda * etaf * H)

V1 = substitute(a - b, list(a = Vm1, b = Vp1))
V2 = substitute(a - b, list(a = Vm2, b = Vp2))
V3 = substitute(a - b, list(a = Vm3, b = Vp3))
V4 = substitute(a - b, list(a = Vm4, b = Vp4))

f11 = D(F1, "E"); f12 = D(F1, "I"); f13 = D(F1, "H") 
     f14 = D(F1, "F")
f21 = D(F2, "E"); f22 = D(F2, "I"); f23 = D(F2, "H") 
     f24 = D(F2, "F")
f31 = D(F3, "E"); f32 = D(F3, "I"); f33 = D(F3, "H") 
     f34 = D(F3, "F")
f41 = D(F4, "E"); f42 = D(F4, "I"); f43 = D(F4, "H") 
     f44 = D(F4, "F")

v11 = D(V1, "E"); v12 = D(V1, "I"); v13 = D(V1, "H")
     v14 = D(V1, "F")
v21 = D(V2, "E"); v22 = D(V2, "I"); v23 = D(V2, "H")
     v24 = D(V2, "F")
v31 = D(V3, "E"); v32 = D(V3, "I"); v33 = D(V3, "H")
     v34 = D(V3, "F")
v41 = D(V4, "E"); v42 = D(V4, "I"); v43 = D(V4, "H")
     v44 = D(V4, "F")

gammah = 1/5 * 7
gammaf = 1/9.6 * 7
gammar = 1/10 * 7
chi = 1/2 * 7
etaf = 1/4.6 * 7
etar = 1/5 * 7
paras = list(S = 1,E = 0, I = 0, H = 0, F = 0,R = 0,
     sigma = 1/7*7, Theta = 0.81, Lambda = 0.81, betai = 0.588, 
     betah = 0.794, betaf = 7.653, N = 1, gammah = gammah,
     gammaf = gammaf, gammar = gammar, etaf = etaf, 
     etar = etar, chi = chi)

f = with(paras, 
matrix(c(eval(f11), eval(f12), eval(f13), eval(f14),
     eval(f21), eval(f22), eval(f23), eval(f24),
     eval(f31), eval(f32), eval(f33), eval(f34),
     eval(f41), eval(f42), eval(f43), eval(f44)),
     nrow = 4, byrow = T))

v = with(paras, 
matrix(c(eval(v11), eval(v12), eval(v13), eval(v14),
     eval(v21), eval(v22), eval(v23), eval(v24),
     eval(v31), eval(v32), eval(v33), eval(v34),
     eval(v41), eval(v42), eval(v43), eval(v44)),
     nrow = 4, byrow = T))

max(eigen(f %*% solve(v))$values)
