# -*- coding: utf-8 -*-

import matplotlib.pyplot as plt
from ExtSource import ExternalSource
from Material import Material
from Mesh_1D import Mesh_1D
from FEM_linear import FEM_linear
from GL_quad import GL_Quadrature

#  select which problem to solve:
problem_type='two_zones'

# create quadrature
GL=GL_Quadrature(2, False)

if problem_type=='homogeneous_v1':
    # create mesh: widths/subdivision/matID/srcID
    mesh = Mesh_1D([100.], [100], [0], [0],printout=False)
    # create material: D/Siga
    MAT=Material(mesh,[2.],[0.5])
    # create source
    Q=ExternalSource(mesh,[10.])   
    # create BC
    left  = {"type":"dir", "val":5.}
    right = {"type":"dir", "val":5.}
    bc = [left, right]
elif problem_type=='two_zones':
    # create mesh: widths/subdivision/matID/srcID
    mesh = Mesh_1D([5., 13.], [20, 60], [0, 1], [0, 1],printout=False)
    # create material: D/Siga
    MAT=Material(mesh,[1.3, 43],[1.4, 0.4])
    # create source
    Q=ExternalSource(mesh,[5.2, 3.2])
    # create BC
    left  = {"type":"dir", "val":1.}
    right = {"type":"rob", "val":0.}
    bc = [left, right]
else:
    raise Exception('problem type unknown')

# create FEM
FE = FEM_linear(GL, mesh, MAT, Q, bc)
sol = FE.assemble_system(printout=False)

# close all figure
plt.close('all')

plt.figure(0)
plt.plot(mesh.x, sol, 'b-o', ms=2.5)
plt.xlabel('x')
plt.ylabel('flux')
plt.grid(True)
plt.title('Neutron diffusion solution')
