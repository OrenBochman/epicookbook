import numpy as np
import math
import pandas as pd
import pythran
from matplotlib import pyplot as plt


#get_ipython().run_cell_magic('writefile', '.pythranrc', '[compiler]\ninclude_dirs=/usr/include/openblas')
#get_ipython().run_line_magic('load_ext', 'pythran.magic')

np.random.seed(123)

### Vanilla Python version
def sir(u,parms,t):
    bet,gamm,iota,N,dt=parms
    S,I,R,Y=u
    lambd = bet*(I+iota)/N
    ifrac = 1.0 - math.exp(-lambd*dt)
    rfrac = 1.0 - math.exp(-gamm*dt)
    infection = np.random.binomial(S,ifrac)
    recovery = np.random.binomial(I,rfrac)
    return [S-infection,I+infection-recovery,R+recovery,Y+infection]

def simulate():
    parms = [0.1, 0.05, 0.01, 1000.0, 0.1]
    tf = 200
    tl = 2001
    t = np.linspace(0,tf,tl)
    S = np.zeros(tl)
    I = np.zeros(tl)
    R = np.zeros(tl)
    Y = np.zeros(tl)
    u = [999,1,0,0]
    S[0],I[0],R[0],Y[0] = u
    for j in range(1,tl):
        u = sir(u,parms,t[j])
        S[j],I[j],R[j],Y[j] = u
    return {'t':t,'S':S,'I':I,'R':R,'Y':Y}

#get_ipython().run_line_magic('timeit', 'simulate()')

sir_out = pd.DataFrame(simulate())
sir_out

### Pythran compiled version

# As the above code only uses simple Python and Numpy times, it is straightforward to obtain compiled versions of the code using Pythran.
#get_ipython().run_cell_magic('pythran', '-DUSE_XSIMD -march=native -O3', "\nimport numpy as np\nimport math\n\n#pythran export sirp(float64 list, float64 list, float64)\ndef sirp(u,parms,t):\n    bet,gamm,iota,N,dt=parms\n    S,I,R,Y=u\n    lambd = bet*(I+iota)/N\n    ifrac = 1.0 - math.exp(-lambd*dt)\n    rfrac = 1.0 - math.exp(-gamm*dt)\n    infection = np.random.binomial(S,ifrac)\n    recovery = np.random.binomial(I,rfrac)\n    return [S-infection,I+infection-recovery,R+recovery,Y+infection]\n\n#pythran export simulatep()\ndef simulatep():\n    parms = [0.1, 0.05, 0.01, 1000.0, 0.1]\n    tf = 200\n    tl = 2001\n    t = np.linspace(0,tf,tl)\n    S = np.zeros(tl)\n    I = np.zeros(tl)\n    R = np.zeros(tl)\n    Y = np.zeros(tl)\n    u = [999,1,0,0]\n    S[0],I[0],R[0],Y[0] = u\n    for j in range(1,tl):\n        u = sirp(u,parms,t[j])\n        S[j],I[j],R[j],Y[j] = u\n    return {'t':t,'S':S,'I':I,'R':R,'Y':Y}")
#get_ipython().run_line_magic('timeit', 'simulatep()')

# This is around two orders of magnitude faster than the vanilla Python code.
sir_outp = pd.DataFrame(simulatep())
sir_outp

### Visualisation
sline = plt.plot("t","S","",data=sir_out,color="red",linewidth=2)
iline = plt.plot("t","I","",data=sir_out,color="green",linewidth=2)
rline = plt.plot("t","R","",data=sir_out,color="blue",linewidth=2)

plt.xlabel("Time",fontweight="bold")
plt.ylabel("Number",fontweight="bold")

legend = plt.legend(title="Population",loc=5,bbox_to_anchor=(1.25,0.5))
frame = legend.get_frame()
frame.set_facecolor("white")
frame.set_linewidth(0)
