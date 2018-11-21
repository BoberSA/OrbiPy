#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Nov 14 22:56:15 2018

@author: deniszagorodnev
"""
import math
import numpy as np
import scipy


class correction():
   #def __init__(self):
        #self.algorithm = self.findVPlanes
    

    @staticmethod
    def findVPlanes(model, mu, s0, beta, planes, dv0, **kwargs):
        
        def prop2Planes(model, mu, S, planes, **kwargs):
 
            def _stopPlanes(self, t, s, planes=planes):
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
            prop.set_solout(lambda t, S:_stopPlanes(t, S, planes))
            s_new = prop.integrate(3140.0)
            if ((s_new[0] > planes[1]) or (math.fabs(s_new[1]) > planes[2])):
                return 1
            return 0
        
        
        
        
        
        
        s1 = np.asarray(s0).copy()
        vstart = s1[3:5].copy()
        dv = dv0
        dvtol = kwargs.get('dvtol', 1e-16)

        
        rads = math.radians(beta)
        beta_n = np.array([math.cos(rads), math.sin(rads)])
        print(planes)
        p = prop2Planes(model, mu, s1, planes, **kwargs)
        print(vstart + dv * beta_n, ' that is vstart + dv * beta_n')
        s1[3:5] = vstart + dv * beta_n
        p1 = prop2Planes(model, mu, s1, planes, **kwargs)
        print(p == p1, ' positive check')
        if p == p1 and p == 1:
            dv = -dv
       
        v = dv        

        print(dv, ' - that is dv')
        while math.fabs(dv) > dvtol:
            #print('i am in while')
            s1[3:5] = vstart + v * beta_n
            p1 = prop2Planes(model, mu, s1, planes, **kwargs)
     
            if p1 != p:
                print('i am in if')
                v -= dv
                dv *= 0.5
           # print('i am not in if')
            v += dv
            
        return v * beta_n
