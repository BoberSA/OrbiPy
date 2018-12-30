#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Nov 14 22:56:15 2018

@author: deniszagorodnev
"""
import math
import numpy as np
import scipy
import stop_funcs_p as sfp

class correction_tool():
   #def __init__(self):
        #self.algorithm = self.findVPlanes
    @staticmethod
    def prop2Planes(model, mu, S, planes, **kwargs):
 
            def _stopPlanes(t, s, planes=planes):
                    if ((s[0] < planes[0]) or (s[0] > planes[1]) or (math.fabs(s[1]) > planes[2])):
                        return -1
                    return 0
            prop = scipy.integrate.ode(model.equation.ode)
        
            prop.set_initial_value(S, 0)
            prop.set_f_params(*[mu])
            if 'int_param' in kwargs:
                prop.set_integrator('dopri5', **kwargs['int_param'])
            else:
                    prop.set_integrator('dopri5')
            prop.set_solout(lambda t, s:_stopPlanes(t, s, planes))
            s_new = prop.integrate(3140.0)
            if ((s_new[0] > planes[1]) or (math.fabs(s_new[1]) > planes[2])):
                return 1
            return 0

    def findVPlanes(self, model, mu, s0, beta, planes, dv0, **kwargs):
        
        s1 = np.asarray(s0).copy()
        vstart = s1[3:5].copy()
        dv = dv0
        dvtol = kwargs.get('dvtol', 1e-16)

        
        rads = math.radians(beta)
        beta_n = np.array([math.cos(rads), math.sin(rads)])
        #print(planes)
        p = self.prop2Planes(model, mu, s1, planes, **kwargs)
        #print(vstart + dv * beta_n, ' that is vstart + dv * beta_n')
        s1[3:5] = vstart + dv * beta_n
        p1 = self.prop2Planes(model, mu, s1, planes, **kwargs)
        #print(p == p1, ' positive check')
        if p == p1 and p == 1:
            dv = -dv
       
        v = dv        

        while math.fabs(dv) > dvtol:
            s1[3:5] = vstart + v * beta_n
            p1 = self.prop2Planes(model, mu, s1, planes, **kwargs)
            #print(p1)
            if p1 != p:
                #print('i am in if')
                v -= dv
                dv *= 0.5
                print(dv, ' _____ ', dvtol)
           # print('i am not in if')
            v += dv
            
        return v * beta_n

    def prop2Limits(self, model, y0, lims, **kwargs):
        mu1 = model.mu1
        
        evout = []
        arr = model.integrator.integrate_ode(mu1, y0, [0, 3140.0],\
                    events = lims['left']+lims['right'], out=evout)
        #print(y0,evout)
        if evout[-1][0] < len(lims['left']):
            return 0, arr
        else:
            return 1, arr
    def findVLimits(self, model, y0, beta, lims, dv0, retit=False, maxit=100, **kwargs):

        y1 = np.asarray(y0).copy()
        vstart = y1[3:5].copy()
        dv = dv0
        dvtol = kwargs.get('dvtol', 1e-16)
    
        rads = math.radians(beta)
        beta_n = np.array([math.cos(rads), math.sin(rads)])
       
        p, _ = self.prop2Limits(model, y1, lims, **kwargs)
        y1[3:5] = vstart + dv * beta_n
        p1, _ = self.prop2Limits(model, y1, lims, **kwargs)
    
        if p == p1 and p == 1:
            dv = -dv
       
        v = dv        
        i = 0
        while math.fabs(dv) > dvtol and i < maxit:
            y1[3:5] = vstart + v * beta_n
            p1, _ = self.prop2Limits(model, y1, lims, **kwargs)
     
            if p1 != p:
                v -= dv
                dv *= 0.5

            v += dv
            i += 1
#    print('findv iterations:', i)
#    print('%g'%v, end=' ')
        if retit:
            return v * beta_n, i
        return v * beta_n