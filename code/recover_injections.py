# recover rotation periods from light curves in the simulations directory

import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
import os
from simple_acf import simple_acf

SIM_DIR = "data/simulations"
DATA_DIR = "data"
PLOT_PATH = "data/simulations/plots"


def measure_prot(ids, plot_acf=False):
    """
    Measure the rotation periods of a set of light curves.
    :param plot: (optional)
        boolean: create plot of the acf if True.
    """

    # remove previous results file
    if os.path.exists("periods.txt"):
        os.remove("periods.txt")

    # perform acf on a set light curves
    acf_results = np.zeros((len(ids), 5))  # id, p, height, rvar, localph
    for i, id in enumerate(ids):
        print(i, "of", len(ids))
        x, y = np.genfromtxt(os.path.join(SIM_DIR,
                             "lightcurve_{0}.txt".format(str(id).zfill(4)))).T
        period, acf, lags, rvar, height, localph, lppos, rppos = \
            simple_acf(x, y)
        acf_results[i, :] = np.array([id, period, height, rvar, localph])

        if plot_acf:
            plt.clf()
            plt.subplot(2, 1, 1)
            plt.plot(x, y, "k.")
            plt.xlabel("Time (days)")
            plt.ylabel("Normalised flux")
            plt.subplot(2, 1, 2)
            plt.plot(lags, acf)
            plt.xlabel("lags")
            plt.ylabel("acf")
            plt.axvline(period, color="r")
            plt.axvline(lppos, color=".5")
            plt.axvline(rppos, color=".5")
            plt.savefig(os.path.join(PLOT_PATH, "{0}_acf".format(id)))

        if localph < .1:  # remove stars with low localph
            period = 0

        # append results to file
        with open("periods.txt", "a") as f:
            f.write("{0} {1} {2} {3} {4} {5} {6} \n".format(id, period,
                                                            pmin[i], pmax[i],
                                                            height, rvar,
                                                            localph))


def sigma_clipping(x, y, nsigmas, iterations=10):
    for i in range(iterations):
        b, c = fit_line(x, y)
        model = b * x + c
        rms = np.mean((y - model)**2)**.5
        m = np.abs(y - model) < nsigmas * rms
        x, y = x[m], y[m]
    return x, y, b, c


def fit_line(x, y):
    AT = np.vstack((true_periods, np.ones_like(true_periods)))
    ATA = np.dot(AT, AT.T)
    return np.linalg.solve(ATA, np.dot(AT, periods))


if __name__ == "__main__":

    # ids = range(1004)
    # measure_prot(ids, plot_acf=True)

    # load the truth file
    truth = pd.read_csv(os.path.join(DATA_DIR, "final_table.csv"))
    pmin, pmax = truth["P_MIN"], truth["P_MAX"]
    true_periods = np.array(.5*(pmin + pmax))

    # load the results file
    ids, periods, pmins, pmaxs, height, rvar, localph = \
        np.genfromtxt("periods.txt").T

    # remove unreliable rotation periods
    m = (periods > 0) * (localph > .1)  # remove failed period measurements
    true_periods, periods, localph = true_periods[m], periods[m], localph[m]

    # subtract polynomial
    x, y, b, c = sigma_clipping(true_periods, periods, 2.5, 10)
    p = np.polyfit(x, y-(b*x+c), 1)
    # periods -= np.polyval(p, true_periods)

    # plot measured vs true
    lim = .3  # success limit
    plt.clf()
    xs = np.linspace(0, max(true_periods), 100)
    print(b, c)
    plt.plot(xs, xs, "k--")
    plt.plot(xs, xs+xs*lim, "k--")
    plt.plot(xs, xs-xs*lim, "k--")
    plt.plot(true_periods, periods, "k.")
    plt.plot(true_periods, periods - np.polyval(p, periods), "b.")
    # plt.plot(xs, b*xs + c + np.polyval(p, xs), "r")
    plt.plot(xs, b*xs + c + np.polyval(p, xs), "r")
    # plt.scatter(true_periods, periods, marker=".", c=localph,
    #             cmap="GnBu_r", s=50, edgecolor="")
    # plt.colorbar()
    plt.ylim(0, 100)
    plt.xlim(0, 100)
    plt.savefig("test_lph")

    # calculate success rate
    fractional_diff = np.abs(true_periods - periods) / true_periods
    nsuccess = len(fractional_diff[fractional_diff < lim])
    percent = float(nsuccess) / float(len(true_periods)) * 100
    print(percent, "% success")
