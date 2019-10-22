# -*- coding: utf-8 -*-

import numpy as np
import sys

class GL_Quadrature:
    """ 
    Class for computing Gauss-Legendre quadrature in 1D
    """

    def __init__(self, n_qpts=2, printout=False):
        self.n_qpts = n_qpts
        self.printout = printout

        try:
            if self.n_qpts <= 0:
                raise ValueError('quadrature order must be >0','input values is: '+str(self.n_qpts))
        except ValueError as err:
            print(err.args)
            sys.exit(1)
        
        [self.xq, self.wq] = np.polynomial.legendre.leggauss(self.n_qpts)
        if self.printout:
            print('n_qpts '+str(self.n_qpts))
            print('xq')
            print(self.xq)
            print('wq')
            print(self.wq)
