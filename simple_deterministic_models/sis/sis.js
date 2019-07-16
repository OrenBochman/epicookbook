rk4 = require("https://bundle.run/ode-rk4")
Plotly = require("https://cdn.plot.ly/plotly-latest.min.js")

import {slider} from "@jashkenas/inputs"


function SIS(dydt, y, t) {
  dydt[0] = -beta*y[0]*y[1]/N + mu*y[1];
  dydt[1] = beta*y[0]*y[1]/N - mu*y[1];
}

function simulate(func, step, initial, t0, maxtime){
  var integrator = rk4( initial, func, t0, step )
  var time = []
  var SI = []
  time.push(t0)
  SI.push(initial)
  var t = t0
  while(t<maxtime){
    t += step
    integrator = integrator.step()
    SI.push(copy(integrator.y))
    time.push(t)}

  return {time, SI};
}

function copy(x) {
  return Object.assign({},x)
}

tmax = 100

sis_sol = simulate(SIS, 0.01, [N-I0,I0], 0, 100)

{
 let times = 0;
 let i;
 for (i = 0; i < 1; i++) {
   const start = new Date();
   let sis_sol = simulate(SIS, 0.01, [N-I0,I0], 0, 100);
   const elapsed = new Date() - start;
   times += elapsed;
 }
 return "Average time elapsed over 1 replicates: " + times/1 + " ms";
}

{
 const beta = 1;
 const mu = 0.01
 const test_sol = simulate(SIS, 0.01,[N-I0, I0], 0, tmax)
 return "Final number of infected individuals: " + Math.round(test_sol.SI.map((x)=>{return x[1]}).slice(-1))
}
