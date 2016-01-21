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
import os as os
import math as math


def linesolve(lr, gp, tvec, fr):
    """
    Solve the radial diffusivity equation using the method of lines
    :return: Pressure grid in time and space
    :param lr: radius of liquid front
    :param gp: gas pressure at that front = boundary condition
    :param fr: fluid injection rate.
    """
    def begave(vec):
        """
        Get the average on the first two elements
        :param vec:
        :return: average
        """
        bave = (vec[0] + vec[1]) / 2.0

        return bave

    def rhs(y, t, lr, tv, gp, tbeg_real):
        """
        Defines a rhs, where lhs is simply the time derivative
        1/r { P' + r P""}
        Note that the arrays have no unphysical nodes
        :param y: a single time step of the radial pressure vector
        :param t: Current time in simulation
        :param lr: Fluid front radius
        :param tv: vector of times in the nstep grid(in seconds)
        :param gp: Gas pressure for the boundary condition
        :param tbeg_real: in second odeint, t  starts at zero. Real t starts at tbeg_real
        """

        treal = tbeg_real + t
        # Use this to limit rhs messages
        index2 = int(treal / fdprm.dt_fine)
        index_r = np.searchsorted(fdprm.r, lr[index2] + fdprm.shutin_extra_radius)
        if (t - fdprm.tprev) > fdprm.tmax/10:
            print 'rhs t', treal / 3600 /24, index2, index_r,  ' hours', 'fluid radius: ', lr[index2]+fdprm.shutin_extra_radius, ' searched radius: ', fdprm.r[index_r], 'Gas P: ', gp[index2]
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
        deriv[1:nr-1] = (pp + r[1:nr-1] * ppp) *  prm.k / (prm.phi * prm.mu * prm.cr)
        deriv /= r
        deriv[nr-1] = 0.0        # Dirichlet at boundary. There is no change to initial value
        deriv[0] = deriv[1]      # Neuman at zero. [0] and [1] move lockstep having been set in y0

        # Enforce the lack of change at the liquid radius. Make sure that it does not grow after shutin.
        if treal > fdprm.q_shutintime * 3600 * 24:
            index2 = int(fdprm.q_shutintime * 3600.0 * 24.0 / fdprm.dt_fine)
            index_r = np.searchsorted(fdprm.r, lr[index2] + fdprm.shutin_extra_radius)
            # print 'hit limit', treal/3600/24, fdprm.q_shutintime, '    hours'

        # Handle the boundary condition at large radius, where pressure has to match the gas pressure
        deriv[index_r:] = -fdprm.alpha * (y[index_r:] - 1000*1000*gp[index2])
        deriv[-1] = 0.0


        return deriv

    y0 = np.ones(len(fdprm.r)) * prm.pi                         # y is on the regular grid
    y0[0] = y0[1] + prm.qnorm * begave(fdprm.rp) / begave(fdprm.r) * begave(fdprm.dr)   # Boundary cond set at beginning. derives will preserve it.


    # Loop over times. ODEINT is run twice as we have a shutin period at the end
    h0 = 0.001
    y_bef, output = odeint(rhs, y0, fdprm.tvals_bef, h0=h0, hmax=2000.0, mxstep=2000, full_output=True, args=(lr, tvec, gp, 0.0,))
    print 'Done first ODEINT', np.shape(y_bef), len(fdprm.tvals_bef)

    y0 = y_bef[len(y_bef) - 1, :]              # After the shutin, we have no more fluid flow in boundary condition
    y0[1] = y0[0]

    fdprm.prev = 0.0
    tbeg_real = fdprm.tvals_bef[-1]
    y_aft, output = odeint(rhs, y0, fdprm.tvals_aft, h0=h0, hmax=2000.0, mxstep=2000, full_output=True, args=(lr, tvec, gp, tbeg_real))
    print 'Done second ODEINT', np.shape(y_aft), len(fdprm.tvals_aft)

    y_all = np.concatenate([y_bef, y_aft[0:, :]], axis=0)
    plot_result(fdprm.r, fdprm.tvals, y_all, fig=10, compare=True)          # Plot results on IRREGULAR grid
    plot_result_points(fdprm.r, fdprm.tvals, y_all, fig=11, compare=True, plot_r=[5.0, 10.0, 20.0, 40.0, 80.0, 160.0, 320.0])    # Plot results at spatial points
    return y_all

def plot_result_points(r, tvals, yp, fig=10, compare=False, plot_r=[100.0]):
    """
    Plot the results and compare with ei if desired
    :param r: radius vector
    :param tvals: time values to compare
    :param yp: result of fd simulation
    :param fig: figure number
    :param compare: True to compare with analytic
    :param plot_r: list of spatial points to plot
    :return:
    """
    plt.figure(fig)
    ax = plt.subplot(1,2,1)
    plt.title('FD Well Pressure(MPa)')
    plt.xlabel('Time(days)')
    plt.ylabel('Pressure(MPa)')
    for pr in plot_r:
        index3 = np.searchsorted(fdprm.r, pr)
        y = np.zeros(len(tvals))
        for (i, t) in enumerate(tvals):
            y[i] = yp[i, index3]
        plt.plot(tvals/3600.0/24.0, y / 1000 / 1000, 'o-', label=str(pr)+' meters')
    plt.grid()
    legend = plt.legend(bbox_to_anchor=(1.1, 1.0))
    for label in legend.get_texts():
        label.set_fontsize('x-small')

    if compare:
        plt.subplot(1,2,2, sharex=ax, sharey=ax)
        plt.title('Exact FD Well Pressure(MPa)')
        plt.xlabel('Time(days)')
        plt.ylabel('Pressure(MPa)')
        for pr in plot_r:
            Pexact, x = [], []
            for (i, t) in enumerate(tvals):
                t_zerorate = t - fdprm.q_shutintime* 3600.0 * 24.0
                s = prm.D * pr**2 / t
                s_zr = prm.D * pr**2 / t_zerorate
                x.append(t / 3600.0 / 24.0)
                if (t <= fdprm.q_shutintime * 3600.0 * 24.0):
                    Pexact.append((prm.pi + prm.p0 * ei(s))/1000.0/1000.0)
                else:
                    Pexact.append((prm.pi + prm.p0 * ei(s) - prm.p0 * ei(s_zr)) / 1000 / 1000)       # In MPa
            plt.plot(x, Pexact, label=str(pr)+' meters')
        plt.grid()
        legend = plt.legend(bbox_to_anchor=(1.1, 1.0))
        for label in legend.get_texts():
            label.set_fontsize('x-small')
    tit = make_title( {'cr': prm.cr, 'ex_radius': fdprm.shutin_extra_radius, 'h': prm.h})
    plt.suptitle(tit)

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
        plt.plot(r, yp[i] / 1000 / 1000, label=str(round(t/3600/24, ndigits=4)) +' days')
    plt.grid()
    legend = plt.legend(bbox_to_anchor=(1.1, 1.1))
    for label in legend.get_texts():
        label.set_fontsize('x-small')


    if compare:                     # Plot comparison with exact solution
        plt.subplot(1,3,3, sharex=ax, sharey=ax)
        plt.title('Exact Well Pressure')
        plt.xlabel('Radius(m)')
        plt.ylabel('Pressure(MPa)')
        for (i, t) in enumerate(tvals):
            Pexact = []
            t_zerorate = t - fdprm.q_shutintime* 3600.0 * 24.0
            for (j, rp) in enumerate(r):
                s = prm.D * rp**2 / t
                s_zr = prm.D * rp**2 / t_zerorate
                #print 'times: ', t/(3600.0 * 24.0), fdprm.q_shutintime, s, s_zr
                if (t <= fdprm.q_shutintime * 3600.0 * 24.0):
                    Pexact.append((prm.pi + prm.p0 * ei(s))/1000.0/1000.0)
                else:
                    Pexact.append((prm.pi + prm.p0 * ei(s) - prm.p0 * ei(s_zr)) / 1000 / 1000)       # In MPa
            plt.plot(r, Pexact, label=str(round(t/3600/24, ndigits=4))+' days')
        plt.grid()
        legend = plt.legend(bbox_to_anchor=(1.1, 1.1))
        for label in legend.get_texts():
            label.set_fontsize('x-small')
    tit = make_title( {'cr': prm.cr, 'ex_radius': fdprm.shutin_extra_radius, 'h': prm.h})
    plt.suptitle(tit)

def make_title(titdict):
    """
    Make a super title from the input dictionary
    :param titdict:{label, value} entries
    :return: string title
    """

    label = 'Parameters:  '
    for k, v in titdict.iteritems():
        label += '  ' + k + '=' + str((v)) + ';'
    return label


def ei(x):
    """
    exponential integral used in transient analysis
    :param x: input arg
    :return: ei
    """
    # if x > 1000.0 or x < 1.0/1000.0:
    #     print 'anomalous x in ei', x, expn(1, x)

    ei = expn(1, x)
    return ei


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
    plt.loglog(xv, av, 'b-', label='log approximation', basex=10)
    plt.grid(which='both')
    plt.minorticks_on()
    plt.ylim([0.1, 10.0])
    plt.xlabel('s')
    plt.title('log-log plot')
    plt.legend()

    plt.subplot(1,2,2)
    plt.plot(xv, yv, 'r-', label='ei')
    plt.plot(xv, av, 'b-', label='log approximation')
    plt.grid()
    plt.minorticks_on()
    # plt.ylim([0.1, 10.0])
    plt.xlim([0,.2])
    plt.legend()
    plt.xlabel('s')
    plt.title('linear plot')
    plt.suptitle('Test of exponential integral')
    plt.show()

def get_press_pred(timeshift, best_slope, pbase):
    """
    Get a pressure prediction from timeshift using our slope model
    :param timeshift: a vector of timeshifts
    :param best_slope: slope in MPa/ms
    :param pbase: Where the pressure is at the beginning
    :return: pressure prediction
    """

    usetime = timeshift / best_slope + pbase

    return usetime
# =====================================================================================================================
# STUFF suited to the hysteresis modeling
def predict_shots(best_slope, timeshift, tsname, usename, tsbase=-0.2, pbase=4.0, figno=10, dataplot='on'):
    """
    Use the best slope to predict the fluid shot data and compare
    :param best_slope  Whatever the best slope is
    :param timeshift: array of timeshifts
    :param tsname: list of possible names
    :param usename: the name we use
    :param tsbase: baseline timeshift
    dataplot: plot actual fluid shots or not. If not, we still get pressure conversion
    :return:
    """

    """
    fluid shots:        day     offset      numberDay       Pressure
                        8/15    1457        7               4294
                        9/10    1483        33              4473
                        10/18   1521        71              4921
    """

    fl = {'day': [7, 33, 71], 'P': [4.294, 4.473, 4.921]}
    s = np.shape(timeshift)
    for (i, name) in enumerate(tsname):
        if name == usename:
            usetime = get_press_pred(timeshift[i] - timeshift[i][0], best_slope, pbase)

    x = np.arange(0, s[1])
    plt.figure(figno)
    plt.plot(x, usetime, 'ro-')
    if dataplot=='on':
        plt.plot(fl['day'], fl['P'], 'go', markersize=20)
    plt.xlabel('day number')
    plt.ylabel('computed pressure')
    plt.ylim([4, 11])
    plt.title('Fluid shots versus computed timeshift for: ' + usename + ' using slope: ' + str(1/best_slope) )
    plt.grid()
    plt.show()

    sys.exit()

def get_hyst_data(droot, fn_p, slope = 1/5.3,  alignday=75, correctday=83, f=''):
    """
    Get the data for selected areal points for the hysteresis analysis
    Files must have one header line: fn_p has well designation where the pressure is measured and start date
                                     fn_ts has a label for each areal point (column) and start date
                                     no date column. assume input data is agreeing!
    :param droot: root for data files
    :param fn_p: filename pressure
    :param alignday: puts in an extra timeshift for aligning curves at 20 days (fluid effect)
    :param correctday: corrects the fluid column effects in pressure data
    :return: pressure and timeshift vectors
    """

    f1 = open(droot+'ts' + f + '.dat', 'r')
    f2 = open(droot+fn_p, 'r')
    h1 = f1.readline().split()
    h2 = f2.readline().split()
    if len(h1) < 2 or len(h2) != 2:
        print "bad data file headers: ", h1, h2
        sys.exit()
    print h1, h2
    tsname = h1[1:]
    pname = h2[1]
    print pname, tsname
    nwell = len(tsname)
                           #
    TSDAT = []
    for line in f1.readlines():
        a0 = line.split()
        a1 = []
        for i in range(nwell):
            a1.append(float(a0[i]))
        TSDAT.append(a1)

    PDAT = []
    for line in f2.readlines():
        PDAT.append(float(line))

    # if len(PDAT) != len(TSDAT):
    #     print 'Different data file sizes for pressure and timeshift: ', f, len(PDAT), len(TSDAT)
    #     sys.exit()


    nday = len(TSDAT)
    a = np.array(TSDAT)
    q = []
    # align and correct the data
    for i in range(nwell):
        q.append(np.transpose(TSDAT)[i])
        q[i] -= q[i][alignday]
    for i in range(nday):
        if i >= correctday and i < len(PDAT):
            PDAT[i] = PDAT[i] + 0.7

    plot_data(q, PDAT, tsname, alignday, slope, fig=0)          # Consolidated plot of the data

    return q, PDAT, tsname

def plot_data( qD, PDATD, tsnameD, alignday, slope, fig=0, tit=''):
    """
    Simple consolidated plot of the hysteresis data
    :param qD: timeshift lists
    :param PDATD: pressure data
    :param tsnameD: names of the data in qD
    :return:
    """
    nwell = len(qD)
    nday = len(qD[0])
    npress = len(PDATD)
    nqD = np.zeros([nwell, nday])
    for i in range(nwell):
        for j in range(nday):
            nqD[i,j] = -qD[i][j]

    nyp = 3
    nxp = 3        # Three rows. Upper is basic plots and lower is the indiv well crossplots Last is verification
    # Individual attributes
    plt.figure(fig+1, figsize=(15,10))
    ax = plt.subplot(3,3,1)
    d = np.arange(1, nday + 1)
    plt.scatter(d[:npress], PDATD, s=30, c=d[:npress], edgecolors='none')
    plt.ylim([0, 12])
    plt.xlabel('day')
    plt.ylabel('Pressure(MPa)')
    plt.title('Pressure')
    plt.grid()

    plt.subplot(3, 3, 2, sharex=ax)
    for i in range(nwell):
        plt.plot( d, qD[i], 'o-', label=tsnameD[i])
    plt.xlabel('day')
    plt.ylabel('Timeshift(ms)')
    plt.title('Timeshifts')
    plt.ylim([-1.0, 0.6])
    plt.xlim([0, 160])
    lg=plt.legend()
    for label in lg.get_texts():
        label.set_fontsize('x-small')
    plt.grid()


    # Crossplots
    ax1 = plt.subplot(3, 3, 3)
    for i in range(nwell):
        plt.plot(nqD[i][:npress], PDATD, 'o-', label=tsnameD[i], markersize=5)
    plt.ylabel('Pressure(MPa)')
    plt.xlabel('Timeshift(ms)')
    plt.title("Crossplot")
    plt.xlim([ -0.6, 1.2])
    plt.ylim([0, 12])
    lg=plt.legend()
    for label in lg.get_texts():
        label.set_fontsize('x-small')
    plt.grid()

    plt.subplot(3, 3, 4, sharex=ax1, sharey=ax1)
    ax2 = []
    for i in range(nwell):
        if i > 0:
            ax2.append(plt.subplot( 3, 3, 4 + i, sharex=ax1, sharey=ax1))
        plt.scatter( nqD[i][:npress], PDATD, s=30, c=d[:npress])# , edgecolors='none')
        plt.annotate(tsnameD[i], xy=(0.5, 0.1), xycoords='axes fraction', fontsize=16)
        plt.ylabel('Pressure(MPa)')
        plt.xlabel('Timeshift(ms)')
        plt.title("Crossplot")
        plt.xlim([-0.6, 1.2])
        plt.ylim([0, 12])
        plt.grid()

    plt.subplot(3, 3, 7, sharex=ax, sharey=ax)
    ax3 = []
    for i in range(nwell):
        if i > 0:
            ax3.append(plt.subplot( 3, 3, 7 + i, sharex=ax, sharey=ax))
       #  plt.scatter(d, PDATD, s=30, c=d, edgecolors='none', label='3108 pressure')
        Pdata = get_press_pred(qD[i] - qD[i][0], slope, 4.0)
        Pdata = Pdata + PDATD[alignday] - Pdata[alignday]
        plt.plot(d, Pdata, 'ro-', label='ts predict: ' + tsnameD[i])
        plt.plot(d[:npress], PDATD, 'bo-', label='3108 pressure')
        plt.ylabel('Pressure(MPa)')
        plt.xlabel('day number')
        plt.title("Validate model")
        plt.ylim([0, 12])
        lg=plt.legend()
        for label in lg.get_texts():
            label.set_fontsize('x-small')
        plt.grid()

    plt.suptitle(tit)
    plt.show()

    # output the predicted pressure data to a file
    for i in range(nwell):
        fname = '../pressPred_' + tsnameD[i] + '.dat'
        f = open(fname, 'w')
        f.write('day      timeshift    predictedPressure\n')
        Pdata = get_press_pred(qD[i] - qD[i][0], slope, 4.0)
        Pdata += PDATD[alignday] - Pdata[alignday]
        for j in range(nday):
            f.write("%f %f %f\n" % (j, qD[i][j], Pdata[j]))
        f.close()
        os.system('ls')
        os.system('pwd')

# STUFF suited to the old exact solution.
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
