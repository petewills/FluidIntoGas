__author__ = 'nlpwi4'
import parm as prm
import numpy as np
import library as lib
import pylab as plt

"""
Compute pressures using the radial diffusion equation, for input into ndi
Used to simulate pressure in 31-08
Creates an ascii version of an int file
"""

# lib.testei(fig=30)

ofile = open("pressure.dat", "w")

# s = phi * mu * cr * r^2 / (4 * k * t)
# P = q * mu / (4*pi*k*h) Ei*-s)
# The classic drawdown shape comes from the Ei at negative argument. My Ei function is really Ei(-x)

mb = prm.k / prm.mu
exfac = 40.0
D = prm.phi * prm.cr * prm.mu / (4.0 * prm.k) * exfac

Pup_days = 74                   # Days of pressure-up
daystep = 1

xwell = 511258
ywell = 6246989
tlim = [1, 74]         # in days
nt = int((tlim[1] - tlim[0]) / float(daystep))

rmax = 250.0            # generate data out to this radius
dx = 5.0
nxmap = int(2.0 * rmax / dx) + 1

print 'D: ', D, D/24/3600

Pdown = np.zeros([nt+1, nxmap, nxmap])
Pup = np.zeros([nt+1, nxmap, nxmap])
Pafter = np.zeros([nt+1, nxmap, nxmap])
Pderiv = np.zeros([nt+1, nxmap, nxmap])
S = np.zeros([nt+1, nxmap, nxmap])
for t in range(nt):
    dt = ((t-1) * daystep + 0.001) * 3600 * 24
    for i in range(nxmap):
        x = (xwell - rmax) + float(i) * dx
        for j in range(nxmap):
            y = (ywell - rmax) + float(j) * dx
            r = np.sqrt((x-xwell)**2 + (y-ywell)**2)

            s = D * r**2 / dt
            Pdown[t, i, j] = prm.p0 * lib.ei(s) / 1000 / 1000        # In MPa

            dt_up = Pup_days * 3600 * 24
            s = D * r**2 / dt_up
            Pup[t, i, j] = prm.p0 * lib.ei(s) / 1000 / 1000#  - Pdown[t, i, j]       # In MPa

            s = D * r**2 / dt
            Pafter[t, i, j] = Pup[t, i, j] - Pdown[t, i, j]  + prm.pi      # In MPa

            # radial derivative
            s = D * r**2 / dt
            Pderiv[t, i, j] = -2.0 * prm.p0 * np.exp(-s) / r / 1000.0 / 1000.0      # In MPa

            s = D * r**2 / dt
            S[t, i, j] = s                                                  # Arg to Ei



plt.figure(1)
nx = int(np.sqrt(nt))
t = tlim[0]
ind = tlim[0]
ip = 1
"""
for i in range(nx):
    for j in range(nx):
        if ip>1:
            plt.subplot(nx, nx, ip, sharex=ax, sharey=ax)
        else:
            ax = plt.subplot(nx, nx, ip)
        plt.imshow(np.transpose(Pdown[ind]), origin='lower', clim=[prm.pi, 11])
        # plt.contour(Pdown[ind])
        plt.colorbar()
        fr = plt.gca()
        fr.axes.get_xaxis().set_visible(False)
        fr.axes.get_yaxis().set_visible(False)
        fr.text(0.5, 0.85, 'Day '+str(t), transform=fr.transAxes, fontsize=12, horizontalalignment='center')
        t += daystep
        ind += 1
        ip += 1
plt.suptitle('Pressure down - absolute')


plt.figure(2)
nx = int(np.sqrt(nt))
t = tlim[0]
ind = tlim[0]
ip = 1
for i in range(nx):
    for j in range(nx):
        if ip>1:
            plt.subplot(nx, nx, ip, sharex=ax, sharey=ax)
        else:
            # ax = plt.subplot(nx, nx, ip)
            plt.subplot(nx, nx, ip, sharex=ax, sharey=ax)
        plt.imshow(np.transpose(Pup[ind]), origin='lower')#, clim=[prm.pi, 11])
        plt.colorbar()
        fr = plt.gca()
        fr.axes.get_xaxis().set_visible(False)
        fr.axes.get_yaxis().set_visible(False)
        fr.text(0.5, 0.85, 'Day '+str(t), transform=fr.transAxes, fontsize=12, horizontalalignment='center')
        t += daystep
        ind += 1
        ip += 1
plt.suptitle('Pressure up - absolute')

plt.figure(3)
nx = int(np.sqrt(nt))
t = tlim[0]
ind = tlim[0]
ip = 1
for i in range(nx):
    for j in range(nx):
        if ip>1:
            plt.subplot(nx, nx, ip, sharex=ax, sharey=ax)
        else:
            # ax = plt.subplot(nx, nx, ip)
            plt.subplot(nx, nx, ip, sharex=ax, sharey=ax)
        plt.imshow(np.transpose(Pafter[ind]), origin='lower')#, clim=[prm.pi, 11])
        plt.colorbar()
        fr = plt.gca()
        fr.axes.get_xaxis().set_visible(False)
        fr.axes.get_yaxis().set_visible(False)
        fr.text(0.5, 0.85, 'Day '+str(t), transform=fr.transAxes, fontsize=12, horizontalalignment='center')
        t += daystep
        ind += 1
        ip += 1
plt.suptitle('Pressure down - subtracted Max Pressure Date')
plt.show()
"""


# Plot radial dependence
plt.figure(5)
rad = [2.0, 4.0, 8.0, 16.0, 32.0, 64.0, 128.0, 264.0, 528.0, 1024.0]
nx = len(rad)
for i in range(nx):
    # if ip>1:
    #     axn = plt.subplot(1, nx, ip, sharex=ax1, sharey=ax1)
    # else:
    #     ax1 = plt.subplot(1, nx, ip)
    Tp, P = [], []
    for t in range(nt):
        dt = (tlim[1]) * 3600 * 24
        s = D * rad[i]**2 / dt
        Pup = prm.pi + prm.p0 * lib.ei(s) / 1000 / 1000

        dt = ((t+1) * daystep) * 3600 * 24
        s = D * rad[i]**2 / dt
        Pdn = prm.p0 * lib.ei(s) / 1000 / 1000
        P.append(Pup - Pdn)
        Tp.append(t+1)
        # if (i==nx-1):
        #     plt.legend()
        # fr = plt.gca()
        # fr.text(0.5, 0.85, 'Radius '+str(rad[i]), transform=fr.transAxes, fontsize=12, horizontalalignment='center')
    plt.plot(Tp, P, 'o-', label='Radius:'+str(rad[i]))
    print 'abs: ', Tp
    print 'abs: ', P
plt.title('Pressure-up relative to background: Time dependence for a given radius')
plt.xlabel('Days of Injection')
plt.ylabel('Pressure-up(MPa)')
plt.legend()
plt.grid()
plt.show()

# Plot radial dependence
plt.figure(4)
nx = int(np.sqrt(nt))
t = tlim[0]
ind = tlim[0]
ip = 1
for i in range(nx):
    for j in range(nx):
        r= []
        yv, yvU, yvD, Sv, yder = [], [], [], [], []
        for k in range(nxmap):
            for l in range(nxmap):
                x = float(k-nxmap/2) * dx
                y = float(l-nxmap/2) * dx
                r.append(np.sqrt(x**2 + y**2))
                yv.append(Pafter[ind, k, l])
                yvD.append(Pdown[ind, k, l])
                yvU.append(Pup[ind, k, l])
                yder.append(Pderiv[ind, k, l])
                Sv.append(S[ind, k, l])
        if ip>1:
            axn=plt.subplot(nx, nx, ip, sharex=ax1, sharey=ax1)
        else:
            ax1 = plt.subplot(nx, nx, ip)
        # plt.plot(Pdown[ind, :, nxmap/2], 'r-', label='Pdown')
        # plt.plot(Pup[ind, :, nxmap/2], 'k-', label='Pup')
        plt.scatter(r, yv, c='b', edgecolors='none', label='Chng post shutin', s=10)
        plt.scatter(r, yvD, c='r', edgecolors='none', label='Pressure-up', s=10)
        plt.scatter(r, yvU, c='k', edgecolors='none', label='Pressure-dwn', s=10)
        if (i==nx-1) and (j==nx-1):
            plt.legend()
        # if ip>1:
        #     axn2 = axn.twinx()
        # else:
        #     axn2 = ax1.twinx()
        # plt.scatter(r, yder, c='g', edgecolors='none', label='P deriv', s=10)

        #     plt.scatter(r, Sv, c='g', edgecolors='none', label='Ei Arg') # , 'ro', , markersize=2)
        fr = plt.gca()
        fr.text(0.5, 0.85, 'Day '+str(t), transform=fr.transAxes, fontsize=12, horizontalalignment='center')
        plt.grid()
        t += daystep
        ind += 1
        ip += 1
plt.suptitle('Radial dependence')
plt.show()


day0 = 1524             # Day when we shutin
ofile.write("%16s%16s\n"% ( Pup, Pdown) )
for it in range(nt):
    for j in range(nxmap):
        for k in range(nxmap):
            x = (xwell - rmax) + float(j) * dx
            y = (ywell - rmax) + float(k) * dx
            ofile.write("%16.2f%16.2f%16d%16.2f%16.2f\n" % (x, y, it+day0, Pup[i, j, k], Pdown[i, j, k], Pafter[i, j, k]))


ofile.close()



