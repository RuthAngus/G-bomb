# recover rotation periods from light curves in the simulations directory

import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
import os
from simple_acf import simple_acf

SIM_DIR = "data/simulations"
DATA_DIR = "data"
PLOT_PATH = "data/simulations/plots"


def recover(plot_acf=False):
    """
    Measure the rotation periods of a set of light curves.
    :param plot: (optional)
        boolean: create plot of the acf if True.
    """

    # remove previous results file
    if os.path.exists("periods.txt"):
        os.remove("periods.txt")

    # load the truth file
    truth = pd.read_csv(os.path.join(DATA_DIR, "final_table.csv"))
    pmin, pmax = truth["P_MIN"], truth["P_MAX"]

    # perform acf on simulated light curves
    ids = range(1004)
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

        # append results to file
        with open("periods.txt", "a") as f:
            f.write("{0} {1} {2} {3} {4} {5} {6} \n".format(id, period,
                                                            pmin[i], pmax[i],
                                                            height, rvar,
                                                            localph))


if __name__ == "__main__":
    recover(plot_acf=False)

#     ids, periods, pmins, pmaxs, height, rvar = np.genfromtxt("periods.txt").T

#     plt.clf()
#     plt.scatter(.5*(pmins + pmaxs), periods, marker=".", c=np.log(rvar), s=50,
#                 edgecolor="")
#     plt.colorbar()
#     plt.ylim(0, 100)
#     plt.xlim(0, 100)
#     plt.savefig("test_rvar")
