__author__ = 'nlpwi4'
import parm as prm
import numpy as np
import library as lib
import math
import pylab as plt

"""
Compute pressures using the radial diffusion equation
This version does a numerical solution that allows us to vary boundary conditions.
The boundary condition will be a tank of gas.

1/r d/dr (r dP / dr ) = phi mu c / k dP / dt = m dP / dt
"""

gsat = 0.2              # Gas saturation
tankrad = 190.0         # radius of stimulated reservoir
tankvol = (tankrad**2) * math.pi * prm.h * prm.phi * gsat     # Initial volume of gas

nday = 80               # days to run simulation
tstep = 3600.0          # time step in seconds
nstep = int(float(nday) / tstep * 24.0 * 3600.0)

cumul = 0.0
x, lr, gp = [], [], []
for i in range(nstep):
    cumul += prm.q * tstep               # Cumulative injected volume.
    liqrad = np.sqrt(cumul / gsat / prm.phi / prm.h / math.pi)  # Fluid radius. Used for boundary condition

    if liqrad < (tankrad-10.0):
        gasvol = tankvol - cumul              # remaining gas volume
        gaspress = prm.pi * tankvol / gasvol  # gas pressure via ideal gas law

        x.append(float(i*tstep/24.0/3600.0))
        lr.append(liqrad)
        gp.append(gaspress)

P = lib.linesolve()

# plt.figure(21)
# plt.plot(x, lr, 'ro')
# plt.ylabel('Liquid Radius(m)')
# plt.figure(22)
# plt.plot(x, gp, 'ro')
# plt.ylabel('Gas Pressure(MPa)')
# plt.show()


