# find rotation periods of gaia stars

import numpy as np
import matplotlib.pyplot as plt
import pandas as pd


def match(id1, id2):
    """
    id1 = small
    id2 = big
    """
    matched = []
    inds1, inds2 = [], []
    for i, id in enumerate(id1):
        m = id2 == id
        if len(id2[m]):
            matched.append(id)
            inds2.append(int(np.where(m)[0]))
            inds1.append(i)
    return matched, inds1, inds2


# load matched tgas catalogue
tgas = pd.read_csv("ruth_matched.csv")

# load Amy's catalogue
rot = pd.read_csv("data/Table_1_Periodic.txt")

# load my catalogue
kid, periods = np.genfromtxt("kplr_periods.txt").T

# matched, tgas_inds, rot_inds = match(tgas["kepid"], rot["KID"])
matched, tgas_inds, rot_inds = match(tgas["kepid"], kid)

periods = np.array(rot["Prot"])[rot_inds]
teff = np.array(tgas["teff"])[tgas_inds]
logg = np.array(tgas["logg"])[tgas_inds]
pmra = np.array(tgas["pmra"])[tgas_inds]
pmdec = np.array(tgas["pmdec"])[tgas_inds]
parallax = np.array(tgas["parallax"])[tgas_inds]

m = (logg > 4.2) * (teff < 6250)
periods, teff, logg, pmra, pmdec, parallax = periods[m], teff[m], logg[m], \
        pmra[m], pmdec[m], parallax[m]

plt.clf()
plt.plot(pmra, periods, "k.")
plt.xlabel("pmra")
plt.ylabel("period (days)")
plt.savefig("period_pmra")

plt.clf()
plt.plot(pmdec, periods, "k.")
plt.xlabel("pmdec")
plt.ylabel("period (days)")
plt.savefig("period_pmdec")

plt.clf()
plt.scatter(1./parallax, periods, s=50, c=teff, cmap="GnBu", label="Teff")
plt.colorbar()
plt.xlabel("distance (Kpc)")
plt.ylabel("Rotation period (days)")
plt.xlim(0, 2)
plt.savefig("period_parallax")

plt.clf()
plt.hist(pmra)
plt.savefig("pmra_hist")
