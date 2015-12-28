__author__ = 'nlpwi4'
import parm as prm
import fdparm as fdprm
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



cumul = 0.0
x, lr, gp, tvec = [], [], [], []
for i in range(fdprm.nstep_fine):
    cumul += prm.q * fdprm.dt_fine               # Cumulative injected volume.
    liqrad = np.sqrt(cumul / gsat / prm.phi / prm.h / math.pi)  # Fluid radius. Used for boundary condition

    if liqrad < (tankrad-10.0):
        gasvol = tankvol - cumul              # remaining gas volume
        gaspress = prm.pi * tankvol / gasvol  # gas pressure via ideal gas law

        x.append(float(i*fdprm.dt_fine/3600.0))
        lr.append(liqrad)
        gp.append(gaspress/1000000)
        tvec.append(i*fdprm.dt_fine)

plt.figure(21)
ax=plt.subplot(121)
plt.plot(x, lr, 'r-')
plt.xlabel('Time(hours)')
plt.ylabel('Liquid Radius(m)')
plt.grid()
plt.subplot(122, sharex=ax)
plt.plot(x, gp, 'r-')
plt.ylabel('Gas Pressure(MPa)')
plt.xlabel('Time(hours)')
plt.grid()
plt.show()

P = lib.linesolve(lr, gp, tvec)

