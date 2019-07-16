rk4 = require('https://bundle.run/ode-rk4')
Plotly = require("https://cdn.plot.ly/plotly-latest.min.js")
import {slider} from "@jashkenas/inputs"


sir_forced = (dydt, y, t) => {
  let N = y[0] + y[1] + y[2]
  let force_inf = beta0*(1+beta1*Math.sin(omega*t))*y[0]*y[1]/N
  dydt[0] = mu*N - force_inf - mu*y[0]
  dydt[1] = force_inf - (gamma + mu)*y[1]
  dydt[2] = gamma*y[1] - mu*y[2]
}

sir_sol = simulate(sir_forced,0,[99999, 1, 0],tmax)

copy = (x) => {
  return Object.assign({},x)
}

simulate = (f,t0,y0,tmax) => {
  let step = 1
  let integrator = rk4(y0, f, t0, step)
  let t = t0
  let y = copy(y0)
  let ts = []
  let ys = []
  // Simulate to equilibrium
  while (t < 1000*365) {
    t = t+step
    integrator = integrator.step()
  }
  t = 0
  ts.push(t)
  ys.push(copy(integrator.y))
  while (t < tmax) {
    t = t+step
    integrator = integrator.step()
    ys.push(copy(integrator.y))
    ts.push(t)
  }
  return {t:ts, y:ys}
}

// Test
{
  return 'Default parameters: le=70, beta0 = 1.4, beta0 = 0.05, f = 365, gamma = 0.14, tmax_years = 2. I_final = 14.646455760129378'
}

{
  let start = new Date()
  let sir_sol = simulate(sir_forced,0,[99999, 1, 0],tmax)
  let time = new Date() - start
  return `Runs in ${time} milliseconds`
}
