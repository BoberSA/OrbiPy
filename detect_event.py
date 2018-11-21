#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import numpy as np
import matplotlib.pyplot as plt
from numba import compiler, types
import correction

class detection_tool():

    @staticmethod
    def plot(model):
        
        #model.model_info()
        L = model.lagrange_points
       
        L1 = L[0, 0]
        L2 = L[1, 0]
        rtol = 1e-14 
        nmax = 1e6 
        int_param = {'atol':rtol, 'rtol':rtol, 'nsteps':nmax, 'method':'dop853'}
        
        leftp = model.mu1 + 500000 / model.ER
        rightp = L2 + 500000 / model.ER
        topp = 1.0
        planes = [leftp, rightp, topp]
        
        X0km = -277549
        Z0km =  200000
        y0 = np.array([L2 + X0km/model.ER, 0, Z0km/model.ER, 0, 0, 0])
        
        #print('start with - ', y0)
        #v = correction.correction.findVPlanes(self, self.mu1, y0, 90, planes, 0.2, int_param=int_param)
        #y0[3:5] = v
        #print('after correction - ', y0)
        
        right_part = compiler.compile_isolated(model.equation.ode, [types.double, types.double[:], types.double], return_type=types.double[:]).entry_point

        arr = model.integrator.integrate_ode(right_part, model, y0, [0, 2*np.pi], int_param = int_param)

        plt.figure(figsize=(10,10))
        plt.plot(arr[:,0],arr[:,1],'.-')
        plt.axis('equal')
