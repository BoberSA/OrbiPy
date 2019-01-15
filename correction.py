#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Nov 14 22:56:15 2018

@author: deniszagorodnev
"""
import math
import numpy as np


class correction_tool():
    
    def prop2Limits(self, model, y0, lims): 
        evout = []
        arr = model.integrator.integrate_ode(model, y0, [0, 3140.0], events = lims['left']+lims['right'], out=evout)
        if(len(evout)!=0):
            if evout[-1][0] < len(lims['left']):
                return 0, arr
            else:
                return 1, arr
        else: 
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