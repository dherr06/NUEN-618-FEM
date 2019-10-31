# -*- coding: utf-8 -*-
"""
Created on Sun Oct  6 23:35:39 2019

@author: ragusa
"""
import sys
import numpy as np
import scipy.sparse
import scipy.sparse.linalg
from Material import Material as MAT

class FEM_nonlinear:
    
    def __init__(self, GL, MESH, BC, porder=1, printout=False):
        self.porder = porder
        self.printout = printout
        self.n = MESH.n_el * porder + 1
        # The entries below are passed by reference, no need to worry about copies/memory
        self.MESH = MESH
        self.GL = GL
        self.BC = BC

        try:
            if porder != MESH.porder:
                raise ValueError("porder in FEM and porder from MESH are different")
        except ValueError as err:
            print(err.args)
            sys.exit(1)

    """
    def b0(x):
        return (1-x)/2.
    
    def b1(x):
        return (1+x)/2.
    
    def db0dx(x):
        return -0.5
    
    def db1dx(x):
        return 0.5
    
    def bf(argument,x):
        switcher = {
        0: b0,
        1: b1
        }
        func = switcher.get(argument, "invalid")
        return func(x)
    def dbfdx(argument,x):
        switcher = {
        0: db0dx,
        1: db1dx
        }
        func = switcher.get(argument, "invalid")
        return func(x)
    """
    
    def assemble_system(self,T):
        printout = self.printout
    
        # hardcoding linear basis functions in [-1,+1]
        bf = lambda i, x: (1-x)/2. if (i==0) else (1+x)/2. 
        dbfdx = lambda i, x : -0.5 if (i==0) else 0.5
        
        # rhs init
        rhs = np.zeros(self.n)
        
        # CSR matrix init
        rows = []  # row positions of matrix entries
        cols = []  # column positions of matrix entries
        entries = []  # matrix entry values
       
        # elementary matrices
        p = self.porder
        b    = np.zeros((self.GL.n_qpts, p+1))
        dbdx = np.zeros((self.GL.n_qpts, p+1))
        for i in range(p+1):
            for q in range(self.GL.n_qpts):
                b[q,i]    = bf(   i, self.GL.xq[q])
                dbdx[q,i] = dbfdx(i, self.GL.xq[q])
        
        # shortcut
        gn = self.MESH.gn

        for iel in range(self.MESH.n_el):
            
            Jac = self.MESH.dx[iel] /2.
            local_A = np.zeros((p+1, p+1))
            for q in range(self.GL.n_qpts):
                for i in range(p+1):
                    for j in range(p+1):
                        local_A[i,j] += 1/Jac * self.GL.wq[q] * MAT.k((T[iel]*b[q,0]) + (T[iel + 1]*b[q,1])) * dbdx[q,i] * dbdx[q,j]
                        rhs[iel+j] += Jac * self.GL.wq[q] * b[q,i] * b[q,j] * 3e4

            # left BC:
            if iel==0:
                bc_ = self.BC[0]
                if bc_.get('type')=='neu':
                    pass
                else:
                    raise Exception("Unknown left BC type")
        
            # right BC:
            if iel==self.MESH.n_el-1:
                bc_ = self.BC[1]
                if bc_.get('type')=='dir':
                    val = bc_.get('val')
                    local_A[1,0] = 0
                    local_A[1,1] = 1
                    rhs[iel+1] = val
                else:
                    raise Exception("Unknown right BC type")
                    
            rows += [ gn[iel,0], gn[iel,0], gn[iel,1], gn[iel,1] ]
            cols += [ gn[iel,0], gn[iel,1], gn[iel,0], gn[iel,1] ]    
            entries += list(local_A.flatten())
            
        A = scipy.sparse.csr_matrix((entries, (rows, cols)), shape=(self.n, self.n))
        res = A.dot(T) - rhs
        if printout:
            print("\nMatrix A:-----------------------")
            print(A)
            print("\RHS b:-----------------------")
            print(rhs)
        return res