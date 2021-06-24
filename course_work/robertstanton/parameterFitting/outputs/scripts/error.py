#!/usr/bin/env python

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from FOX import from_hdf5

hdf5_file: str = 'armc.hdf5'

#Extract ONLY the errors of the accepted iterations
err: pd.DataFrame = from_hdf5(hdf5_file, 'aux_error')
acceptance: np.ndarray = from_hdf5(hdf5_file, 'acceptance')  # Boolean array
accerr = err[acceptance.values].loc[:, "rdf.0"]
print(accerr)

#Line over the iteration with the lowest error
miny = accerr.min()
maxy = accerr.max()
minx = accerr.idxmin()
plt.vlines(minx, maxy+0.5, miny-0.5, colors='k', linestyles= 'dashed')

#Plot
plt.plot(accerr, label="rdf.0")
plt.legend()
plt.xlabel('ARMC iteration')
plt.ylabel('Auxiliary Error')
plt.show()

