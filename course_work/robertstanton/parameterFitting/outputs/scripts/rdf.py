#!/usr/bin/env python

import pandas as pd
from FOX.recipes import get_best, overlay_descriptor, plot_descriptor

hdf5_file: str = 'armc.hdf5'

param: pd.Series = get_best(hdf5_file, name='param')  # Extract the best parameters
rdf: pd.DataFrame = get_best(hdf5_file, name='rdf')  # Extract the matching RDF

#_dct = overlay_descriptor(hdf5_file, name='rdf')
#desired_keys = ["Cs Cs"]
#dct = {k: _dct[k] for k in desired_keys}
#plot(dct)

# Compare the RDF to its reference RDF and plot
rdf_dict = overlay_descriptor(hdf5_file, name='rdf',  i=0)
plot_descriptor(rdf_dict)
