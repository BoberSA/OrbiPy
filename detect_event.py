#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import numpy as np
import matplotlib.pyplot as plt
from numba import compiler, types
import correction

class detection_tool():
    def __init__(self, model):
        self.model = model
        

    def detect(self, y0, ev1, ev2):
        events = {'left':[ev1], 'right':[ev2]}

        cor = correction.correction_tool()
        dv = cor.findVLimits(self.model, y0, 90, events, 0.05, retit=False, maxit=100)
        y0[3:5] = dv
        self.model.equation = compiler.compile_isolated(self.model.equation, [types.double, types.double[:], types.double], return_type=types.double[:]).entry_point
        evout= []
        arr = self.model.integrator.integrate_ode(self.model, y0, [0, 2*np.pi], events, evout)
        print("that is evout", evout)
        return(arr, evout)
        
    @staticmethod
    def plot(arr, evout):
        plt.figure(figsize=(10,10))
        plt.plot(arr[:,0],arr[:,1],'.-')
        plt.axis('equal')
        ev_names = ['X:0', 'alpha:120', 'Y:0', 'alpha:60']
        for ie, _, s, _ in evout:
            plt.plot(s[0], s[1], '+k')
            plt.text(s[0], s[1], ' [%d] %s' % (ie, ev_names[ie]))
