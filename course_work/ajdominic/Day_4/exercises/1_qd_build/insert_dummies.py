#!/usr/bin/env python

from scm.plams import Molecule
from CAT.recipes import replace_surface
mol = Molecule('cspbbr3_4.2nm.xyz')
mol_new = replace_surface(mol, symbol='Br', symbol_new='Cl', f=0.8, mode='uniform', displacement_factor=0.7)
mol_new.write('cspbbr3_4.2nm_80Cl.xyz')
