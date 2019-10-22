# -*- coding: utf-8 -*-

import numpy as np
import sys

class Material:
    """
    Load material properties
    """
    def __init__(self):
        def k(T):
            K = 1.5 + (2510/(215+T))
            return K
        self.k = k
