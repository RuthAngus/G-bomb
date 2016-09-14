# recover rotation periods from kepler light curves

import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
import os
from simple_acf import simple_acf
from kepler_data import load_kepler_data
import glob

SIM_DIR = "data/simulations"
DATA_DIR = "data"
PLOT_PATH = "data/simulations/plots"


def measure_prot(ids, DATA_DIR, PLOT_DIR, plot_acf=False):
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
        str_id = str(int(id)).zfill(9)
        fnames = glob.glob(os.path.join(DATA_DIR,
                           "{0}/kplr*{0}*llc.fits".format(str_id)))
        print(fnames)
        x, y, yerr = load_kepler_data(fnames)
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
        with open("kplr_periods.txt", "a") as f:
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

    kepler_tgas = pd.read_csv("data/ruth_matched.txt")
