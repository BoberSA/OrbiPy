#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import numpy as np
import matplotlib.pyplot as plt
from numba import compiler, types
import correction

class detection_tool():
    def __init__(self, model):
        self.model = model
        

    def detect(self):
        L = self.model.lagrange_points
       
        L1 = L[0, 0]
        L2 = L[1, 0]

        leftp = self.model.mu1 + 500000 / self.model.ER
        rightp = L2 + 500000 / self.model.ER
        topp = 1.0
        planes = {'left':leftp, 'right':rightp}
        
        X0km = -277549
        Z0km =  200000
        y0 = np.array([L2 + X0km/self.model.ER, 0, Z0km/self.model.ER, 0, 0, 0])
        #y0 = np.array([1.00817646, 0, 0.0013369, 0, 0, 0])
        
        #print('start with - ', y0)
        #v = correction.correction.findVPlanes(self, self.mu1, y0, 90, planes, 0.2, int_param=int_param)
        #y0[3:5] = v
        #print('after correction - ', y0)
        cor = correction.correction_tool()
        y0[3:5] = cor.findVLimits(self.model, y0, 90, planes, 0.2,retit=False, maxit=100)
        #y0[3:5] = cor.findVLimits(model, y0, beta, lims, dv0, retit=False, maxit=100, **kwargs)
        
        self.model.equation = compiler.compile_isolated(self.model.equation, [types.double, types.double[:], types.double], return_type=types.double[:]).entry_point

        return(self.model.integrator.integrate_ode(self.model, y0, [0, 2*np.pi]))
        
    @staticmethod
    def plot(arr):
        
        #model.model_info()
        

        plt.figure(figsize=(10,10))
        plt.plot(arr[:,0],arr[:,1],'.-')
        plt.axis('equal')
