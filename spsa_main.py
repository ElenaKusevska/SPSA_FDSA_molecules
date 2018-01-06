#!/usr/bin/python

import spsa_routines

import sys
import os

import math
import random
import numpy

import openbabel
import pybel
from pybel import ob

random.seed(0)

#-----------------------------------------------------
# Python implemenetation of the SPSA algorithm
# for optimizing the dihedral angles of conformers 
# using MMFF94 energies:
#-----------------------------------------------------

#-----------------------------------------------------
# set up openbabel and read molecule from input file:
#-----------------------------------------------------

# Check if the the correct number of commandline arguments has been used,
# as in:    python3 spsa_main.py  filename.sdf :
if (sys.argv) != 2:
    exit

# read the molecule from the supplied file:
filename = sys.argv[1]
extension = os.path.splitext(filename)[1]
mol = next(pybel.readfile(extension[1:], filename))
print("Molecule read. Number of atoms: ", len(mol.atoms))

# global variable for the MMFF94 force field
force_field = pybel._forcefields["mmff94"]

# Check if the force field can be set up for this molecule:
if force_field.Setup(mol.OBMol) is False:
    print("Cannot set up MMFF94 force field, exiting")
    exit()

# set up the pool of rotatable bonds:
rl = ob.OBRotorList()
rl.Setup(mol.OBMol) # can have values True or False
p = int(rl.Size()) # number of independent variables to optimize
print("Number of rotatable bonds: ", p)

# cycle through the rotatable bonds to create the working array:
rotIterator = rl.BeginRotors()
rotor = rl.BeginRotor(rotIterator) # first rotatable bond
rotors = [] # contains the rotors as c++ objects
angles = numpy.empty(p, dtype=numpy.float64) # contains the 
#                              angles in radians - working array for SPSA
coords_current = mol.OBMol.GetCoordinates()
i = 0
while rotor is not None:
    angle = rotor.CalcTorsion(coords_current)
    angles[i] = angle
    rotors.append(rotor)
    rotors[i].SetToAngle(coords_current, angle)
    rotor = rl.NextRotor(rotIterator)
    i = i + 1

# update the coordinates and compute the energy
mol.OBMol.SetCoordinates(coords_current)
force_field.SetCoordinates(mol.OBMol)
E_old = force_field.Energy(False) # false = don't compute gradients
 
# Write first set of coordinates to file:
outputname1 = 'optimization_steps'
if os.path.isfile('./' + outputname1 + '.xyz'):
    os.remove('./' + outputname1 + '.xyz')
    print( outputname1 + '.xyz', ' found and deleted' )
spsa_routines.writesxyz(mol, outputname1, E_old)

# Write first set of dihedrals and energy to a file:
outputname2 = 'angles_energies.txt'
if os.path.isfile('./' + outputname2):
    os.remove('./' + outputname2)
    print( outputname2, ' found and deleted' )
spsa_routines.writesangles(angles, p, 1, E_old, outputname2)

# Get good converged energies for comparison
targets = numpy.empty(4, dtype=numpy.float64)

#-----------------------------------------------------
# Set up SPSA:
#-----------------------------------------------------

# assign the gain sequence coefficients:
alpha = 0.602
gama = 0.101

#assign the other variables:
c = 0.000005
smallest_initial_change = 0.0003 # smallest desired magnitude change
#                                during the early iterations
                                    
niter = 20.0 # maximum expected/allowed number of iterations
B = 0.09*float(niter)
a = smallest_initial_change * ((B + 1.0)**alpha)

#----------------------------
# the SPSA procedure:
#----------------------------

k = 1 # counter
test = 1.0

outputname3 = 'working_information.txt'
if os.path.isfile('./' + outputname3):
    print( outputname3, ' found and deleted' )
    os.remove('./' + outputname3)

#while test > 0.000001:
while k <= 30:

    spsa_routines.writesarray(angles, p, outputname3, 'old angles')
    
    # compute the k-th coefficients
    ak = a / ((B + k)**alpha)
    ck = c / ((k)**gama)

    # generate the random simultaneous perturbation vector:
    del_k = numpy.empty(p, dtype=numpy.float64)
    del_k = spsa_routines.Bernoulli (del_k,p)
    spsa_routines.writesarray(del_k, p, outputname3, 'del_k')

    # perturb the coordinates and evaluate the function:
    coords_plus = mol.OBMol.GetCoordinates()
    coords_minus = mol.OBMol.GetCoordinates()

    new_angles = numpy.empty(p, dtype=numpy.float64)

    for i in range (0, p, 1):
        new_angles[i] = angles[i] + ck*del_k[i]
        rotors[i].SetToAngle(coords_plus, new_angles[i])
    mol.OBMol.SetCoordinates(coords_plus)
    force_field.SetCoordinates(mol.OBMol)
    E_plus = force_field.Energy(False) 
    spsa_routines.writesarray(new_angles, p, outputname3, 'angles_plus')
        
    for i in range (0, p , 1):
        new_angles[i] = angles[i] - ck*del_k[i]
        rotors[i].SetToAngle(coords_minus, new_angles[i])
    mol.OBMol.SetCoordinates(coords_minus)
    force_field.SetCoordinates(mol.OBMol)
    E_minus = force_field.Energy(False)
    spsa_routines.writesarray(new_angles, p, outputname3, 'angles_minus')

    outputfile = open(outputname3, 'a')
    outputfile.write('{:10} {:10.7f} {:10} {:10.7f}'.format\
            ('E_plus', E_plus, 'E_minus', E_minus))
    outputfile.write('\n')
    outputfile.write('\n')
    outputfile.close()

    # determine the gradient:
    gradient = [0.0]*p
    dE = (E_plus - E_minus) / (2*ck) # multiplicative factor
    for i in range (0, p, 1):
        gradient[i] = dE * (1/del_k[i])
    spsa_routines.writesarray(gradient, p, outputname3, 'gradient')
    
        
    # update the theta estimate, evaluate the function (energy) for the
    # new theta estimate (angles), and write the new coordinates to
    # a file:
    for i in range (0, p, 1):
        angles[i] = angles[i] - ak*gradient[i]
        rotors[i].SetToAngle(coords_current, angles[i])
    mol.OBMol.SetCoordinates(coords_current)
    force_field.SetCoordinates(mol.OBMol)
    E_new = force_field.Energy(False)
    spsa_routines.writesxyz(mol, outputname1, E_new)
    spsa_routines.writesangles(angles, p, k, E_new, outputname2)
    spsa_routines.writesarray(angles, p, outputname3, 'new angles')


    # test for convergence or termination:
    #if (k .ge. 600) break loop
    test = abs(E_new - E_old)

    outputfile = open(outputname3, 'a')
    outputfile.write('{:15} {:15} {:15}'.format('k', 'E_new', 'test'))
    outputfile.write('\n')
    outputfile.write('{:15d} {:15.9f} {:15.9f}'.format(k, E_new, test))
    outputfile.write('\n')
    outputfile.write('\n')
    outputfile.write('---------------------------------------- \n')
    outputfile.write('\n')
    outputfile.close()

    E_old = E_new
    k = k + 1

