# Crossmatching catalogues: produce ruth_matched.csv

import pandas as pd
from cross_match import convert_to_cartesian
from scipy.spatial import cKDTree as KDTree
import numpy as np


def matching(id1, id2):
    matched = []
    for id in id1:
        m = id2 == id
        if len(id2[m]):
            matched.append(id)
    return matched


def load_kepler():
    return pd.read_csv("data/matched.csv")


def load_tgas():
    return pd.read_csv("../../gaia/Tgas.csv")

if __name__ == "__main__":
    tol = np.radians(.5 / 3600)  # arcsec -> deg -> rad

    print("Loading the Tycho-2 catalog...")
    tycho2 = load_tgas()
    tycho2_xyz = convert_to_cartesian(np.array(tycho2[["ra", "dec"]],
                                               dtype=np.float64))
    tid = np.array(tycho2["tycho2_id"])

    print("Loading the Kepler catalog...")
    kepler = load_kepler()
    kepler_xyz = convert_to_cartesian(np.array(kepler[["ra", "dec"]],
                                               dtype=np.float64))

    # Building KD-tree.
    print("Building KD trees...")
    tycho2_tree = KDTree(tycho2_xyz)
    kepler_tree = KDTree(kepler_xyz)

    # Cross match.
    print("Cross matching trees...")
    match = kepler_tree.query_ball_tree(tycho2_tree, np.sqrt(2-2*np.cos(tol)))
    match_flag = np.zeros(len(kepler), dtype=bool)
    match_id = np.zeros(len(kepler), dtype=np.uint64)
    distances = np.nan + np.zeros(len(kepler))
    tycho_id = []
    for i, m in enumerate(match):
        if len(m):
            # Compute the angular distance.
            d = np.arccos(np.dot(kepler_xyz[i], tycho2_xyz[m].T))
            distances[i] = np.min(d)
            if distances[i] <= tol:
                match_id[i] = m[np.argmin(d)]
                match_flag[i] = True

    print(len(match_id))
    TGAS_PATH = "/export/bbq2/angusr/gaia/Tgas.csv"
    kepler["tycho2_match_id"] = \
        np.array(pd.read_csv(TGAS_PATH)["tycho2_id"].iloc[match_id])
    print(len(np.array(pd.read_csv(TGAS_PATH)["tycho2_id"].iloc[match_id])))
    kepler["pmra"] = \
        np.array(pd.read_csv(TGAS_PATH)["pmra"].iloc[match_id])
    kepler["pmdec"] = \
        np.array(pd.read_csv(TGAS_PATH)["pmdec"].iloc[match_id])
    kepler["parallax"] = \
        np.array(pd.read_csv(TGAS_PATH)["parallax"].iloc[match_id])
    kepler["parallax_error"] = \
        np.array(pd.read_csv(TGAS_PATH)["parallax_error"].iloc[match_id])
    kepler["tycho2_match_dist_deg"] = np.degrees(distances)
    kepler["tycho2_match_ra"] = np.array(tycho2.ra.iloc[match_id])
    kepler["tycho2_match_dec"] = np.array(tycho2.dec.iloc[match_id])
    print(np.shape(kepler))
    kepler_matched = kepler.iloc[match_flag]
    print(np.shape(kepler_matched))
    kepler_matched.to_csv("ruth_matched.csv", index=False)
