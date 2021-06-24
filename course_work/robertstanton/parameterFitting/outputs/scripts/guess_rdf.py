#!/usr/bin/env python

import pandas as pd
from FOX import MultiMolecule, example_xyz, estimate_lj
xyz_file: str = 'CdSQD.xyz'
atom_subset = ['Cd', 'S']
mol = MultiMolecule.from_xyz(xyz_file)
rdf: pd.DataFrame = mol.init_rdf(atom_subset=atom_subset)
param: pd.DataFrame = estimate_lj(rdf)
print(param)
