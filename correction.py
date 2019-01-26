#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Nov 14 22:56:15 2018

@author: deniszagorodnev
"""
import math
import numpy as np


class correction_tool():
    
    def __init__(self, func=None):
        if func!=None:
            None
            #self.correction = lambda t, s: func(self.model, y0, beta, lims, dv0=0.1, retit=False, maxit=100)
         
        else:
            self.correction = self.findVLimits
        
    
    def prop2Limits(self, model, y0, lims): 
        evout = []
        arr = model.integrator.integrate_ode(model, y0, [0, 3140.0], events = lims['left']+lims['right'], out=evout)
        if(len(evout)!=0):
            #print('evout[-1][0] ', evout[-1][0])
            if evout[-1][0] < len(lims['left']):
                #print('returned ', 0)
                return 0, arr
            else:
                #print('returned ', 1)
                return 1, arr
        else: 
            #print('returned ', 1)
            return 1, arr
        
        
    def findVLimits(self, model, y0, beta, lims, dv0, retit=False, maxit=100):

        y1 = np.asarray(y0).copy()
        vstart = y1[3:5].copy()
        dv = dv0
        dvtol = 1e-16
        rads = math.radians(beta)
        beta_n = np.array([math.cos(rads), math.sin(rads)])
       
        p, _ = self.prop2Limits(model, y1, lims)
        y1[3:5] = vstart + dv * beta_n
        p1, _ = self.prop2Limits(model, y1, lims)
        if p == p1 and p == 1:
            dv = -dv
        v = dv        
        i = 0
        while math.fabs(dv) > dvtol and i < maxit:
            y1[3:5] = vstart + v * beta_n
            p1, _ = self.prop2Limits(model, y1, lims)
     
            if p1 != p:
                v -= dv
                dv *= 0.5

            v += dv
            i += 1
        if retit:
            return v * beta_n, i
        return v * beta_n
    
    def time2Sphere(self, model, y0, beta, lims, dv0=0.1, retit=False, maxit=100):
        import scipy
        def Tsphere(y0, v, model, events):
            evout = []
            y0[4] = v
            model.integrator.integrate_ode(model, y0, [0, 100], events=lims['left']+lims['right'], out=evout)
            if(len(evout)!=0):
                return evout[-1][2][6]
            else:
                return 8.514498422344646

        v_max = scipy.optimize.minimize_scalar(lambda v: -Tsphere(y0, v, model, lims), bracket=(-0.05, 0.05), method='Brent', tol=1e-16)
        return(v_max.x)