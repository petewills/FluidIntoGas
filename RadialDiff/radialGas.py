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

# lib.testei(fig=999)


tankvol = (prm.tankrad**2) * math.pi * prm.h * prm.phi * prm.gsat     # Initial volume of gas

cumul = 0.0
x, lr, gp, tvec, fr, frp = [], [], [], [], [], []
for i in range(fdprm.nstep_fine):

    time_hours = float(i*fdprm.dt_fine/3600.0)
    if time_hours <= (fdprm.q_shutintime * 24.0):
        rate = prm.q
    else:
        rate = 0.0

    cumul += rate * fdprm.dt_fine               # Cumulative injected volume.dt_fine is time step in seconds
    liqrad = np.sqrt(cumul / prm.gsat / prm.phi / prm.h / math.pi)  # Fluid radius. Used for boundary condition

    if liqrad < (prm.tankrad-10.0):
        gasvol = tankvol - cumul              # remaining gas volume
        gaspress = prm.pi * tankvol / gasvol  # gas pressure via ideal gas law

        x.append(time_hours /24)
        lr.append(liqrad)
        gp.append(gaspress/1000000)
        fr.append(rate)
        frp.append(rate * 3600 * 24)
        tvec.append(i*fdprm.dt_fine)


P = lib.linesolve(lr, gp, tvec, fr)
sp = np.shape(P)

plt.figure(21)
ax=plt.subplot(141)
plt.plot(x, lr, 'r-')
plt.xlabel('Time(days)')
plt.ylabel('Liquid Radius(m)')
plt.title('Liquid Radius')
plt.grid()
plt.subplot(142, sharex=ax)
plt.plot(x, gp, 'r-')
plt.ylabel('Gas Pressure(MPa)')
plt.xlabel('Time(days)')
plt.title('Gas Pressure')
plt.grid()
plt.subplot(143, sharex=ax)
plt.plot(x, frp, 'r-')
plt.ylabel('Fluid Rate(m3/day)')
plt.xlabel('Time(days)')
plt.title('Fluid Rate')
plt.ylim([-10, prm.q * (3600*24) + 10])
plt.grid()
plt.subplot(144, sharex=ax)
plt.plot(fdprm.tvals/24/3600, P[:, sp[1]/2]/1000/1000, 'ro-', label='P(rmax/2)')
plt.plot(x, gp, 'b-', label='Gas Pressure')
plt.ylabel('Pressure(MPa)')
plt.xlabel('Time(days)')
plt.title('Boundary condition versus gas pressure at middle radius')
# plt.ylim([0, 10])
plt.legend(bbox_to_anchor=(1.2, 0.1))
plt.grid()
plt.show()
