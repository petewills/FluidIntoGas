__author__ = 'peterwills'
import numpy as np
import sys as sys
import pylab as plt

def irreg_cluster(pwr=1.0, deps=1.0):
    """
    create an irregular grid based on power law map function. Using pwr=1 will give you a regular grid.
    :param pwr: how much to expand the grid for larger r. 1.0 is regular
    :param deps: regular grid interval in meters
    :return: r, rreg, dr, rp, rpp
    """
    rseg = [deps, 500.0]                        # Interval for both grids
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
tmax = 3600*500                           # seconds to run simulation
nt = 20
dt = int(tmax / nt)                       # time step in seconds
nstep = int(tmax / dt)                    # Number of steps to run
tvals = np.arange(0.0, tmax, dt)
tprev = 0.0

# We have a very fine time grid used to find liquid radius in rhs
dt_fine = 10.0          # Fine grid for search of lr in rhs
nstep_fine = int(tmax / dt_fine)

# r, dr = reg(delr=0.1)
r, rreg, dr, rp, rpp = irreg_cluster(pwr=1.5, deps=5.0)

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
plt.grid()
plt.legend()
# plt.show()






"""
def irreg():

    # First node is first physical point in the model
    #
    rseg = [0.1, 2.0, 10.0, 100.0, 1000.0]          # at nodes
    drseg   = [0.5, 1.0,  10.0,   100.0,   1000.0]  # also at nodes
    nseg = 10                                       # single segment n - applies to all segments
                                      # single segment n - applies to all segments
    Nseg = len(rseg) - 1                            # Number of segments
    r, dr = np.zeros(Nseg*(nseg-1)), np.zeros(Nseg*(nseg-1))
    for i in range(Nseg):
        lr, ur = rseg[i], rseg[i+1]
        ldr, udr = drseg[i], drseg[i+1]
        temp, tempr = [], []
        for j in range(nseg-1):
            lamb = (float(j) / float(nseg-2)) **1.6        # 0 to 1
            temp.append((1.0-lamb) * ldr + lamb * udr)
        norm = np.sum(temp)
        tot = lr
        tempr = [lr]
        for j in range(nseg-1):
            temp[j] *= (ur - lr) / norm
        for j in range(nseg-1):
            tot += temp[j]
            tempr.append(tot)
        r[(nseg-1)*i:(nseg-1)*(i+1)] = tempr[:nseg-1]
    for j in range(len(dr)-1):
        dr[j] = r[j+1] - r[j]
    dr[-1] = dr[-2]

    return r, dr
"""
"""
def irreg1(expand=1.0, delr=0.1):
    ""
    create a regular grid
    :param expand: how much to expand the grid for larger r. 1.0 is regular
    :param delr: starting interval
    :return: r, dr
    ""
    rseg = [0.1, 500.0]         # Interval for grid
    rp = rseg[0]
    sz, N = delr, 0
    r, dr = np.zeros(100000), np.zeros(100000)
    while rp < rseg[1]:
        r[N] = rp
        dr[N] = sz
        rp = rp + sz
        sz *= expand
        N += 1

    r = r[:N-1]
    dr = dr[:N-1]

    return r, dr
        """