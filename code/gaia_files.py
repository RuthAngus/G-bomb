from __future__ import print_function
import numpy as np
import pandas as pd
import glob
import matplotlib.pyplot as plt

# fnames = glob.glob("*csv")
# db = pd.read_csv(fnames[0])
# print(np.shape(db))
# for f in fnames[1:]:
#     db = pd.concat([db, pd.read_csv(f)])
#     print(np.shape(db))
#
# pd.DataFrame.to_csv(db, "gaia.csv")

tgas = pd.read_csv("../../gaia/Tgas.csv")
tid = tgas["tycho2_id"]
kic = pd.read_csv("../../gaia/gaia-kepler/matched.csv")
ktid = kic["tycho2_match_id"]

# m = kic["tycho2_match_dist_deg"] < 10
print(tid[:10])
print(ktid[:10])
assert 0
print(np.shape(ktid))
print(np.shape(ktid)[m])

plt.clf()
plt.hist(kic["radius"], 50)
plt.savefig("hist")

proper_motions = []
for i, id in enumerate(ktid[m]):
    m = id == tid
    proper_motions.append
