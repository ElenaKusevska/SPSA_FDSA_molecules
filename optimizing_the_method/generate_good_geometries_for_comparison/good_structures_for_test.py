#!/usr/bin/python

#----------------------------------------------------------------
# Program to generate a series of good/optimized results
# from the starting geometry for comparison - a force
# field steepest descent, a force field conjugate
# gradients. All these structures
# and respective geometries can be used as target
# geometries and energies for a genetic algorithm, for
# optimizing the optimization method, for comparison with the
# geometries optimized using the optimization method. 
#
# You may have the idea of comparing against DFT energies/geometries
# again. Please remember that this is not a good idea. DFT 
# optimization is carried out using the DFT model. You are using
# a force field model here. The results are not comparable!
#
# Rotor search results are also not comparable.
#----------------------------------------------------------------

import sys
import os

import math
import random
import numpy

import openbabel
import pybel
from pybel import ob

#----------------------------------------------------
# function to write the coordinates to a .xyz file:
#----------------------------------------------------

def writesxyz(molecule, fname, energy):
    fname = fname + '.xyz'
    outputfile = open(fname, 'a')
    outputfile.write('{:d}'.format(len(molecule.atoms)))
    outputfile.write('\n')
    outputfile.write('{:.12f}'.format(energy)) #there was a string on energy
    outputfile.write('\n') 
    for atom in molecule:
        outputfile.write('{:8} {:20.9f} {:20.9f} {:20.9f}'.format\
                ( str(ob.etab.GetSymbol(int(atom.atomicnum))), \
                atom.coords[0], atom.coords[1], atom.coords[2] ))
        outputfile.write('\n')
    outputfile.close()

#---------------------------------------
# Main program:
#---------------------------------------

# global variable for the MMFF94 force field
force_field = pybel._forcefields["mmff94"] 

# Check if the the correct number of commandline arguments has been used,
# as in:    dihedrals.py  filename.sdf :
if (sys.argv) != 2:
    exit

# read the molecule from the supplied file:
# (several separate optimizations, so several
# separate mol1/mol2/mol3 objects)
# (It is not possible to use ne molecule object,
# but shuffle around several coordinates arrays.
# and reset the coordinates that way. Every time
# a new geometry optimization is started, if
# you want it to be starting from the geometry
# of the input file, you must use a new molecule
# object, unless you have not updated the
# coordinates after the previous optimization. But
# then you can't print the geometries as you go
# along, so that's no good...
filename = sys.argv[1]
extension = os.path.splitext(filename)[1]
mol1 = next(pybel.readfile(extension[1:], filename))
mol2 = next(pybel.readfile(extension[1:], filename))
print("Molecule read. Number of atoms: ", len(mol1.atoms))

# Check if the force field can be set up for this molecule before
# continuing further:
if force_field.Setup(mol1.OBMol) is False:
    print("Cannot set up MMFF94 force field, exiting")
    exit()

outputname = 'good_geometries'

# Sarting geometry:
force_field.SetCoordinates(mol1.OBMol)
E = force_field.Energy(False)
print('starting geometry:', E)
writesxyz(mol1, outputname, E)

# Steepest descent:
force_field.SteepestDescent(1000)
force_field.UpdateCoordinates(mol1.OBMol)
E_conv = force_field.Energy(False)
print('steepest descent:', E_conv)
writesxyz(mol1, outputname, E_conv)

# Starting geometry:
force_field.SetCoordinates(mol2.OBMol)
E = force_field.Energy(False)
print('starting geometry', E)
writesxyz(mol2, outputname, E)

# Conjugate gradients:
force_field.ConjugateGradients(1000)
force_field.UpdateCoordinates(mol2.OBMol)
E_conv = force_field.Energy(False)
print('conjugate gradients', E_conv)
writesxyz(mol2, outputname, E_conv)

