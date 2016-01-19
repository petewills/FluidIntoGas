__author__ = 'peterwills'
import numpy as np
import sys as sys
import pylab as plt
import parm as prm

def irreg_cluster(pwr=1.5, deps=1.0):
    """
    create an irregular grid based on power law map function. Using pwr=1 will give you a regular grid.
    :param pwr: how much to expand the grid for larger r. 1.0 is regular
    :param deps: regular grid interval in meters
    :return: r, rreg, dr, rp, rpp
    """
    rseg = [deps, 1000.0]                        # Interval for both grids
    delr = rseg[1] - rseg[0]                    # step for uniform grid
    rreg = np.arange(rseg[0], rseg[1], deps)    # The uniform grid
    nreg = len(rreg)                            # number of nodes

    r, dr = np.zeros(100000), np.ones(100000)*deps
    rp, rpp = np.zeros(10000), np.zeros(10000)
    for i in range(nreg):
        lam = rreg[i] / delr         # (0, 1) to the power
        rp[i] = pwr * lam**(pwr - 1)
        rpp[i] = pwr * (pwr - 1) / delr * (lam**(pwr - 2))
        r[i] = rseg[0] + (lam**pwr) * delr


    r = r[:nreg]
    dr = dr[:nreg]
    rp = rp[:nreg]
    rpp = rpp[:nreg]

    return r, rreg, dr, rp, rpp


# Time stepping for the model
# Time is measured in seconds.
tmax = 3600*24*120                           # total seconds to run simulation
nt = 50
dt = int(tmax / nt)                       # time step in seconds

# Relaxation constant to honour gas pressure outside the liquid sphere
# This constant will multiply the pressure difference in derivs, measured im MPa.
# Should be small as possible whilst maintaining pressure outside liquid(which is plotted)
alpha = 0.0001

# Handle the shutin option
q_shutintime = 80.0                                               # shutin time in days. Make sure its not too big or else gas vanishes
t_init = 3600*24*4; dt_init = 20000
t_bef_shutin = q_shutintime * 3600.0 * 24.0
t_aft_shutin = tmax - t_bef_shutin
if t_aft_shutin <= 3.0*dt or t_bef_shutin <= 3.0*dt:
    print 'need times after and before shutin!', t_bef_shutin,  t_aft_shutin, 3.0*dt
    sys.exit()
shutin_extra_radius = 5.0          # Add this to the liquid radius to get the gas boundary point. Prevents impossible pressureup.

tvals_init = np.arange(dt_init, t_init, dt_init)          # Needed fine at start due to strong gas pressure curvature
tvals_bef = np.arange(t_init, t_bef_shutin, dt)
tvals_bef = np.concatenate([tvals_init, tvals_bef])

tvals_aft = np.arange(0.01, t_aft_shutin + dt, dt)
tvals = np.concatenate([tvals_bef, tvals_aft + t_bef_shutin])

nstep = int(tmax / dt)                    # Number of steps to run before and after shutin
tprev = 0.0
flag = 0

# We have a very fine time grid used to find liquid radius in rhs
dt_fine = 10.0          # Fine grid for search of lr in rhs
nstep_fine = int(tmax / dt_fine)

# r, dr = reg(delr=0.1)
deps = 5.0
r, rreg, dr, rp, rpp = irreg_cluster(pwr=1.5, deps=deps)

plt.figure(61)
ax = plt.subplot(1,3,1)
plt.title("Irregular Grid Map Jacobian rp")
plt.xlabel('r(m)')
plt.ylabel('rp')
plt.plot(r, rp, 'ro')
plt.grid()

plt.subplot(1,3,2, sharex=ax)
plt.title("Irregular Grid Map Jacobian rpp")
plt.xlabel('r(m)')
plt.ylabel('rpp')
plt.plot(r, rpp, 'ro')
plt.grid()


plt.subplot(1,3,3, sharex=ax)
plt.title("Irregular Grid vs Regular grid")
plt.xlabel('r(m)')
plt.ylabel('dr')
plt.plot(r, dr, 'ro', label='irregular')
plt.plot(rreg, np.ones(len(r)), 'bo', label='regular')
plt.ylim([0.0, deps+1.0])
plt.grid()
plt.legend()
# plt.show()
