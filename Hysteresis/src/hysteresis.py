__author__ = 'peterwills'
"""
Simple model to explain cross-plots of well pressure and seismic timeshift.
Convention: slowdown is negative
Required input files for data: ts.dat(timeshifts) and p.dat (pressure)
"""
import sys as sys
py = "/glb/data/eptr_am_8/Wills/Python/PyCharmProjects/"
sys.path.append(py + "FluidIntoGas/RadialDiff")
DATAROOT = py + "FluidIntoGas/Hysteresis/"
import pylab as plt
import numpy as np
import library as lib

alignday=75
correctday = 83
q, pdat, tsname = lib.get_hyst_data(DATAROOT, "ts.dat", "p.dat", alignday=alignday, correctday=correctday)
nwell = len(q)

nday = len(pdat)

# p_init, p_max = 4.0, 11.0            # MPa
t_init, t_max = 60.0, 120.0          # deg C
ts_p_slope =  1.0 / 5.7              # ms/MPa. Positive is speedup. Positive P is negative. Read off soak data slope.
ts_T_slope =  1.0 / 150.00           # ms / deg C. Positive T is negative.

T, P, TS, d = [], [], [], []
for i in range(nday):
    frac = float(i)/float(nday)
    T.append( t_init + (t_max - t_init) * frac)             # T is linear increase
    # x = 2.0 * (0.5 - frac)                                  # -1 to 1
    P.append(pdat[i])
    # P.append(  p_max - x**2 * (p_max - p_init)  )           # Quadratic pressure
    ts = ts_p_slope * (P[i] -  P[0] ) + ts_T_slope * ( T[i] - T[0] )
    TS.append( ts )
    # d.append(i)
TSCAL = TS[alignday]
for i in range(nday):
    TS[i] = TS[i] - TSCAL

lib.plot_data([TS], pdat, tsname, fig=1)          # Consolodated plot of the data
plt.show()

# # Individual attributes
# plt.figure(1)
# ax = plt.subplot(131)
# # plt.plot(d, P, 'ro-')
# plt.scatter( d, P, s=30, c=d, edgecolors='none')
# plt.xlabel('day')
# plt.ylabel('Pressure(MPa)')
# plt.title('Pressure')
# plt.grid()
#
# plt.subplot(132, sharex=ax)
# plt.scatter( d, T, s=30, c=d, edgecolors='none')
# plt.xlabel('day')
# plt.ylabel('Temperature(C)')
# plt.title('Temperature')
# plt.grid()
#
# plt.subplot(133, sharex=ax)
# plt.scatter( d, TS, s=30, c=d, edgecolors='none')
# plt.xlabel('day')
# plt.ylabel('Timeshift(ms)')
# plt.title('Timeshift')
# plt.grid()
#
#
# # Crossplots
# plt.figure(2)
# ax = plt.subplot(111)
# plt.scatter( TS, P, s=30, c=d, edgecolors='none')
# plt.ylabel('Pressure(MPa)')
# plt.xlabel('Timeshift(ms)')
# plt.title('P-TS')
# plt.ylim([0,12])
# plt.xlim([0.2, -1.6])
# # plt.gca().invert_xaxis()
# plt.grid()


plt.show()