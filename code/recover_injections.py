# recover rotation periods from light curves in the simulations directory

import numpy as np
import matplotlib.pyplot as plt

SIM_DIR = "data/simulations"

ids = range(1004)
for id in ids:
    x, y = np.genfromtxt(os.path.join(SIM_DIR,
                         "lightcurve_{0}".format(str(id).zfill(4)))
