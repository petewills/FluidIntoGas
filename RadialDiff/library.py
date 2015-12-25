__author__ = 'peterwills'
""" Library for the transient well functionality"""
from scipy.special import expn as expn
import pylab as plt
import numpy as np
import fdparm as fdprm
import parm as prm
from scipy.integrate import ode as ode
from scipy.integrate import odeint as odeint
import sys as sys
import math as math

def linesolve():
    """
    Solve the radial diffusivity equation using the method of lines
    :return: Pressure grid in time and space
    """
    def begave(vec):
        """
        Get the average on the first two elements
        :param vec:
        :return: average
        """
        bave = (vec[0] + vec[1]) / 2.0
        #bave = vec[0]

        return bave

    def rhs(y, t):
        """
        Defines a rhs, where lhs is simply the time derivative
        1/r { P' + r P""}
        Note that the arrays have no unphysical nodes
        :param y: a single time step of the radial pressure vector
        :param t: Current time in simulation
        """

        # Use this to limit rhs messages
        if t > fdprm.tprev:
            print 'rhs t', t / 3600, ' hours'
            fdprm.tprev = t

        nr = len(y)             # All points are physically in the model
        y0 = y[1:nr-1]          # Chop two points at either end to form vector differences
        rp = fdprm.rp[1:nr-1]
        rpp = fdprm.rpp[1:nr-1]
        r = fdprm.r                         # irregular grid
        dr = fdprm.dr[1:nr-1]               # always regular grid version
        ym = y[0:nr-2]
        yp = y[2:nr]

        pp_eps = (y0 - ym)/dr
        ppp_eps = ((yp - 2.0*y0 + ym)/dr**2)
        pp = pp_eps / rp
        ppp = ppp_eps / rp**2 - pp_eps * rpp / rp**3

        deriv = np.zeros(nr)
        deriv[1:nr-1] = (pp + r[1:nr-1] * ppp) * prm.nu
        deriv /= r
        deriv[nr-1] = 0.0        # Dirichlet at boundary. There is no change to initial value
        deriv[0] = deriv[1]      # Neuman at zero. [0] and [1] move lockstep having been set in y0

        return deriv

    y0 = np.ones(len(fdprm.r)) * prm.pi                         # y is on the regular grid
    y0[0] = y0[1] + prm.qnorm * begave(fdprm.rp) / begave(fdprm.r) * begave(fdprm.dr)   # Boundary cond set at beginning. derives will preserve it.

    # Loop over times
    h0 = 0.001
    y, output = odeint(rhs, y0, fdprm.tvals, h0=h0, hmax=2000.0, mxstep=2000, full_output=True)


    plot_result(fdprm.r, fdprm.tvals, y, fig=10, compare=True)          # Plot results on IRREGULAR grid

def plot_result(r, tvals, yp, fig=10, compare=False):
    """
    Plot the results and compare with ei if desired
    :param r: radius vector
    :param tvals: time values to compare
    :param yp: result of fd simulation
    :param fig: figure number
    :param compare: True to compare with analytic
    :return:
    """
    plt.figure(fig)
    plt.subplot(1,3,1)
    plt.title('FD Well Pressure(MPa)')
    plt.xlabel('Radial node')
    plt.ylabel('Time step')
    plt.imshow(yp, aspect='auto')
    plt.grid()
    plt.colorbar()

    ax=plt.subplot(1,3,2)
    plt.title('FD Well Pressure')
    plt.xlabel('Radius(m)')
    plt.ylabel('Pressure(MPa)')
    for (i, t) in enumerate(tvals):
        plt.plot(r, yp[i] / 1000 / 1000, label=str(t/3600.0)+' hours')
    plt.grid()
    plt.legend()

    if compare:                     # Plot comparison with exact solution
        plt.subplot(1,3,3, sharex=ax, sharey=ax)
        plt.title('Exact Well Pressure')
        plt.xlabel('Radius(m)')
        plt.ylabel('Pressure(MPa)')
        for (i, t) in enumerate(tvals):
            Pexact = []
            for (j, rp) in enumerate(r):
                D = 1.0 / (4.0 * prm.nu)
                s = D * rp**2 / t
                Pexact.append((prm.pi + prm.p0 * ei(s)) / 1000 / 1000)       # In MPa
            plt.plot(r, Pexact, label=str(t/3600.0)+' hours')
        plt.grid()
        plt.legend()
    plt.show()

def ei(x):
    """
    exponential integral used in transient analysis
    :param x: input arg
    :return: ei
    """

    ei = expn(1, x)
    return ei


def onewell(pi, p0, D, r, t):
    """
    Compute pressure at r from a well at origin
    :param r: distance from well
    :param D: diffusivity
    :param P0: Pressure multiplier
    :param pi: injector pressure
    :param tv: times to calculate
    :return:Pressures at origin
    """
    arg = D * r * r / t
    p = pi + p0 * ei(arg) / 1000000.0       # pressure in MPa

    return p

def onewell_times(x, y, D, p0, pi, tv, fault='NULL'):
    """
    Compute pressure at test point from well at origin
    :param x: x coordinate test point
    :param y: y coordinate of test point
    :param D: diffusivity
    :param P0: Pressure multiplier
    :param pi: injector pressure
    :param tv: times to calculate(seconds)
    :param fault: x coordinate of fault
    :return:Pressures at origin
    """

    n = len(tv)
    pv = np.zeros(n)
    for i in range(n):
        r = np.sqrt(x*x + y*y)
        t = tv[i]

        pv[i] = onewell(pi, p0, D, r, t)
        if fault[0] == 'FAULT':
            x1 = x - 2.0 * fault[1]
            r1 = np.sqrt(x1*x1 + y*y)
            pv[i] += onewell(0.0, p0, D, r1, t)
        elif fault[0] == 'CIRCLE':   # generate a bunch of wells to give circular reservoir limit.
            nc = fault[2]
            c = 2.0 * 3.14159 * (2.0 * fault[1])
            dth = 2.0 * 3.1415 / float(nc)
            th = 0.0
            for j in range(nc-1):
                xc = x - np.cos(th) * 2.0 * fault[1]
                yc = y - np.sin(th) * 2.0 * fault[1]
                rc = np.sqrt(xc*xc+yc*yc)
                pv[i] += onewell(0.0, p0, D, rc, t) / float(nc)
                #print x, xc, rc, pv[i], th
                th += dth

    #sys.exit()
    return pv

def onewell_xy(xv, yv, D, p0, pi, T):
    """
    Compute pressure at origin from a well at this position
    :param x: x coordinate f well vector in meters
    :param y: y coordinate of well vector
    :param D: diffusivity
    :param P0: Pressure multiplier
    :param pi: injector pressure
    :param T: time of calculation in seconds
    :return:Pressures at origin
    """

    n = len(xv)
    pv = np.zeros(n)
    for i in range(n):
        r = np.sqrt(xv[i]*xv[i] + yv[i]*yv[i])

        pv[i] = onewell(pi, p0, D, r, T)

    return pv



def testei(fig=1):
    """
    Test the exponential integral
    :return:
    """
    xv = np.arange(0.001,1.0, 0.001)
    yv = []
    av = []
    for (i, x) in enumerate(xv):
        yv.append(ei(x))
        av.append(-np.log(1.781*x))

    plt.figure(fig)
    plt.subplot(1,2,1)
    plt.loglog(xv, yv, 'r-', label='ei', basex=10)
    plt.loglog(xv, av, 'b-', label='log', basex=10)
    plt.grid(which='both')
    plt.minorticks_on()
    plt.ylim([0.1, 10.0])
    plt.legend()

    plt.subplot(1,2,2)
    plt.plot(xv, yv, 'r-', label='ei')
    plt.plot(xv, av, 'b-', label='log')
    plt.grid()
    plt.minorticks_on()
    # plt.ylim([0.1, 10.0])
    plt.xlim([0,.2])
    plt.legend()
    plt.suptitle('Test of exponential integral')
    plt.show()