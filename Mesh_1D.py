# -*- coding: utf-8 -*-

import numpy as np
import sys

class Mesh_1D:

    def __init__(self, width, zone_subdivision, material_zone, source_zone, porder=1, printout=False):
        self.width = width
        self.n_zones = len(zone_subdivision)
        self.n_el = sum(zone_subdivision)
        self.x = np.zeros([self.n_el + 1])
        self.dx = np.zeros([self.n_el]) 
        self.porder = porder
        self.printout = printout

        # sort user's list of width
        self.width.sort()

        # check whether 0 needs to be included
        if self.width[0] > sys.float_info.epsilon:
            self.width = np.insert(self.width, 0, 0.0)

        # checks:
        try:
            if len(material_zone) != len(source_zone):
                raise ValueError("material_zone and source_zone of different lengths")
        except ValueError as err:
            print(err.args)
            sys.exit(1)
        try:
            if self.n_zones != len(source_zone):
                raise ValueError("zone_subdivision and source_zone of different lengths")
        except ValueError as err:
            print(err.args)
            sys.exit(1)

        # loop over zones to compute cell thicknesses
        x = np.linspace(self.width[0], self.width[1], zone_subdivision[0] + 1)
        for i in range(1, self.n_zones):
            work = np.linspace(self.width[i], self.width[i + 1], zone_subdivision[i] + 1)
            x = np.insert(x, len(x), work[1::])
        self.x = x

        # compute cell thickness
        for i in range(self.n_el):
            self.dx[i] = self.x[i + 1] - self.x[i]
            
        # create connectivity array
        gn = np.zeros((self.n_el, self.porder+1), dtype=int)
        gn[0,:] = np.linspace(0, self.porder, self.porder+1)
        for iel in range(1,self.n_el):
            gn[iel,0] = gn[iel-1,-1]
            gn[iel,1:] = gn[iel-1,1:]+ self.porder
        self.gn = gn
        
        if self.printout:
            print("One-D mesh as follows:")
            print("x")
            print(self.x)
            print("dx")
            print(self.dx)
            print("connectivity array")
            print(gn)

