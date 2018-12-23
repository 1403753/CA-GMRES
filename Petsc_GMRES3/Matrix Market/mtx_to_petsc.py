#! /usr/bin/env python

import os
import sys

petsc_dir = '/usr/lib/petsc'
os.environ['PETSC_DIR'] = petsc_dir
sys.path.append(petsc_dir + '/bin')

import PetscBinaryIO
import scipy.io

matrix_name = sys.argv[1]
vec_name = sys.argv[2]

matrix = scipy.io.mmread(matrix_name + '.mtx')
with open(matrix_name + '_petsc', 'w') as petsc_file:
    PetscBinaryIO.PetscBinaryIO().writeMatSciPy(petsc_file, matrix)

vector = scipy.io.mmread(vec_name + '.mtx')
with open(matrix_name + '_petsc', 'a') as petsc_file:
    PetscBinaryIO.PetscBinaryIO().writeVec(petsc_file, vector)
