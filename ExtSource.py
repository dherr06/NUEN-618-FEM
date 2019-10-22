# -*- coding: utf-8 -*-

import numpy as np
import sys

class ExternalSource:
    """
    Load source properties
    """

    def __init__(self, MESH, qext, printout=False):

        # Assign qext to array
        self.qext = np.asarray(qext)
        
        # check of length of D/Siga is compatible with length of the 
        #   unique entries in mesh.iel2mat
        isrc_unique = np.unique(MESH.iel2src)       
        try:
            if len(isrc_unique) != len(qext):
                raise ValueError("isrc_unique and qext of different lengths")
        except ValueError as err:
            print(err.args)
            sys.exit(1)

        