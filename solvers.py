import numpy as np
import scipy.sparse.linalg as spla
class non_linear:
    
    def Newton(x,func,epsilon=1e-8,lin_solve='direct',tol=1e-5):
        '''
        Performs a non linear solve using a finite difference Jacobian
        '''
        index = 0
        diff = 1
        while diff > tol:
            J = np.zeros((x.shape[0],x.shape[0]))
            for i in range(x.shape[0]):
                e = np.zeros(x.shape[0])
                e[i] = 1
                J[:,i] = (func(x + (e*epsilon)) - func(x))/epsilon
            if lin_solve == 'direct':
                R = np.linalg.solve(J,-func(x))
            elif lin_solve == 'gmres':
                R = spla.gmres(J,-func(x))[0]
            x_new = x + R
            diff = np.abs(np.linalg.norm(x_new) - np.linalg.norm(x))
            x = x_new
            index += 1
            if index == 99:
                print('Max iterations hit, exiting')
                sys.exit()
        return x_new