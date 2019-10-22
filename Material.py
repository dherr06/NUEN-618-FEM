# -*- coding: utf-8 -*-

import numpy as np
import sys

class Material:
    """
    Load material properties
    """

    def __init__(self, MESH, D, Siga, printout=False):

        # checks:
        try:
            if len(D) != len(Siga):
                raise ValueError("D and Siga of different lengths")
        except ValueError as err:
            print(err.args)
            sys.exit(1)
            
        # Assign material arrays 
        self.D = np.asarray(D)
        self.Siga = np.asarray(Siga)
        
        # check of length of D/Siga is compatible with length of the 
        #   unique entries in mesh.iel2mat
        imat_unique = np.unique(MESH.iel2mat)       
        try:
            if len(imat_unique) != len(D):
                raise ValueError("imat_unique and D/Siga of different lengths")
        except ValueError as err:
            print(err.args)
            sys.exit(1)

        