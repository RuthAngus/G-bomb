# Crossmatching catalogues

import pandas as pd
from astero import astero
a = astero()

def match(id1, id2):
    matched = []
    for id in id1:
        m = id2 == id
        if len(id2[m]):
            matched.append(id)
    return matched

