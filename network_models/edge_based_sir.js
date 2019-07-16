rk4 = require('https://bundle.run/ode-rk4')
Plotly = require("https://cdn.plot.ly/plotly-latest.min.js")
import {slider} from "@jashkenas/inputs"


sir_config_net = (dydt, y, t) => {
  dydt[0] = -beta*y[0] + beta*dpsi(y[0],k)/dpsi(1,k) + gamma*(1-y[0])
  let S = psi(y[0],k)
  let I = 1-S-y[1]
  dydt[1] = gamma*I
}

psi = (theta,k) => {
  return Math.pow(theta,k)
}

dpsi = (theta,k) => {
  return k*Math.pow(theta,k-1)
}

sir_sol = simulate(sir_config_net,0,[0.999,0],tmax)

copy = (x) => {
  return Object.assign({},x)
}

simulate = (f,t0,y0,tmax) => {
  let step = tmax/200
  let integrator = rk4(y0, f, t0, step)
  let ts = []
  let ys = []
  ts.push(t0)
  ys.push(copy(y0))
  let t = t0
  while (t < tmax) {
    integrator = integrator.step()
    t = t + step
    ts.push(t)
    ys.push(copy(integrator.y))
  }
  return {t:ts, y:ys}
}

// Test
{
  return 'Default parameters: beta = 0.1, gamma = 0.05, k = 5, tmax = 125. theta_final = 0.342509'
}

{
  let start = new Date()
  let sir_sol = simulate(sir_config_net,0,[0.999, 0],tmax)
  let time = new Date() - start
  return `Runs in ${time} milliseconds`
}
