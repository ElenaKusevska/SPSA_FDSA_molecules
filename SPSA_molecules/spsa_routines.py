#!/usr/bin/python

import random
import numpy

import openbabel
import pybel
from pybel import ob

import os

#random.seed(0)
#numpy.random.seed(0)

#----------------------------------------------------------------------
# function for determining the lowest energy conformer using different
# methods that are available in openbabe. These conformers/energies can
# used as a test for the SPSA optimization.
def bestconformers(testmolecule, bestenergies):
    force_field = pybel._forcefields["mmff94"]
    force_field.SetCoordinates(testmolecule.OBMol)

    outputxyz = open(OpenBabelEnergies.xyz, 'a')
    outputtxt = open(MethodsUsed.txt, 'a')
    #determine energy
    #write

# Basically, need to have a series of minima / structures and energies,
# for a given initial structure, so that I can compare the results, and
# identify good calculations / combinations of input parameters.
#-----------------------------------------------------------------------

# function to write the coordinates for each step in the optimization
# to an .xyz file:
def writesxyz(molecule, fname, energy, overwrite):
    fname = fname + '.xyz'
    
    # If overwriting, delete file if it already exists:
    if overwrite == True:
        if os.path.isfile('./' + fname):
            os.remove('./' + fname)
            print( fname, ' found and deleted' )

    # edit .xyz file:
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

#function to write the values of the dihedrals and the energy at each
# step in the optimization:
def writesangles(angles, dim, step, energy, fname):
    outputfile = open(fname, 'a')
    outputfile.write('{:10} {:d}'.format('step:', step))
    outputfile.write('\n')
    outputfile.write('{:10} {:.14f}'.format('energy:', energy))
    outputfile.write('\n')
    for i in range (0, dim, 1):
        outputfile.write('{}'.format(angles[i]))
        outputfile.write('\n')
    outputfile.write('{}'.format('---------------------------------'))
    outputfile.write('\n')
    outputfile.write('\n')

def writesarray(array, dim, fname, arrayname):
    outputfile = open(fname, 'a')
    outputfile.write('{:10}'.format(arrayname))
    for i in range (0, dim, 1):
        outputfile.write('{:15.10f}'.format(array[i]))
    outputfile.write('\n')
    outputfile.write('\n')


# function to generate the random perturbation vector:
# (by assigning a vector and populating it with random numbers)
def Bernoulli(a, p):
    for i in range (0, p, 1):
        a[i] = random.random()
        if a[i] > 0.5:
            a[i] = 1.0
        elif a[i] <= 0.5:
            a[i] = -1.0
    
    return a

# using Bernoulli(a,p)
 #b = [0.0]*10
 #b = Bernoulli(b, 10)
 #print (b)
  
 #c = numpy.empty(10, dtype=numpy.float64)
 #c = Bernoulli(c, 10)
 #print(c)

# function to generate the random perturbation vector 
#(by creating an array of random numbers):
#def Bernoulli_numpy_array (a, p):
#    for i in range (0, p, 1):
#        print ('{0:.4f}'.format(a[i]))
#        if a[i] > 0.5:
#            a[i] = 1.0
#        elif a[i] <= 0.5:
#            a[i] = -1.0
#        print ('{0:.4f}'.format(a[i]))
#
#    return a

# using Bernoulli_numpy_array(a,p):
#d = numpy.random.random(10)
#d = Bernoulli_numpy_array(d,10)
#print (d)

def split_uniform(a, p):
    for i in range (0, p, 1):
        a[i] = random.random()
        #if a[i] > 0.5:
        #    a[i] = 1.0
        #elif a[i] <= 0.5:
        #     a[i] = -1.0
        
    return a

