# recover rotation periods from light curves in the simulations directory

import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
import os
from simple_acf import simple_acf

SIM_DIR = "data/simulations"
DATA_DIR = "data"


def recover():

    # remove previous results file
    if os.path.exists("periods.txt"):
        os.remove("periods.txt")

    # load the truth file
    truth = pd.read_csv(os.path.join(DATA_DIR, "final_table.csv"))
    N, pmin, pmax = truth["N"], truth["P_MIN"], truth["P_MAX"]

    # perform acf on simulated light curves
    ids = range(1004)
    periods, height, rvar = [], [], []
    for i, id in enumerate(ids):
        print(i, "of", len(ids))
        x, y = np.genfromtxt(os.path.join(SIM_DIR,
                            "lightcurve_{0}.txt".format(str(id).zfill(4)))).T
        periods.append(simple_acf(x, y)[0])
        height.append(simple_acf(x, y)[-1])
        rvar.append(simple_acf(x, y)[3])

        # append results to file
        with open("periods.txt", "a") as f:
            f.write("{0} {1} {2} {3} {4} {5} \n".format(id, periods[i],
                                                        pmin[i], pmax[i],
                                                        height[i], rvar[i]))


if __name__ == "__main__":
    # recover()

    ids, periods, pmins, pmaxs, height, rvar = np.genfromtxt("periods.txt").T

    plt.clf()
    plt.scatter(.5*(pmins + pmaxs), periods, marker=".", c=np.log(rvar), s=50,
                edgecolor="")
    plt.colorbar()
    plt.ylim(0, 100)
    plt.xlim(0, 100)
    plt.savefig("test_rvar")
