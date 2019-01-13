#!/usr/bin/env python3
# -*- coding: utf-8 -*-
#import stop_functions as sf
import numpy as np
import scipy
from stop_funcs_p import stopNull,correctEvents, stopFunCombined

import sympy as sp
sp.init_printing(use_unicode=True)

class integrator_tool:
   def __init__(self, rtol, nmax, method, stopf):
       
        self.rtol = float(rtol) 
        self.nmax = float(nmax) 
        self.method = method 
        self.atol = float(rtol)
        self.stopf = stopf
        self.int_param = {'atol':self.rtol, 'rtol':self.rtol, 'nsteps':self.nmax, 'method':self.method, 'stopf':self.stopf}
   
   def integrate_ode(self, model, s0, tspan, events, out):
    
    retarr=True
    mu = model.mu1
    prop = scipy.integrate.ode(model.equation)

    if self.int_param != None:
        method = self.method
        prop.set_integrator(method, method=self.method, rtol=self.rtol)
    else:
        prop.set_integrator('dopri5')

    prop.set_initial_value(s0, tspan[0])
    prop.set_f_params(*[mu])

    lst = []
    if 'stopf' in self.int_param:
        prop.set_solout(lambda t, s: self.int_param['stopf'](t, s, lst, events, out))
    else:
        prop.set_solout(lambda t, s: stopNull(t, s, lst))
    prop.integrate(tspan[1])
    del prop
    if len(out) > 0:
        cor_out = correctEvents(model, events, out, None, sn=len(s0), int_param=self.int_param)
        out.clear()
        out.extend(cor_out)
      
    if retarr:
 
        return np.asarray(lst)
       
   def set_param(self):
        print('Please choose rtol/nmax/method or Quit')
        name = input()
        
        if(name == 'rtol'):
            rtol = float(input())
            self.rtol = rtol
            
        elif(name == 'nmax'):
            nmax = float(input())
            self.nmax = nmax
            
        elif(name == 'method'):
            method = input()
            self.method = method
        
        elif(name == 'Quit'):
            return(-1)
        
        return(1)
       
   def integrator_help(self):
        print('integrate_ode is integrator of current system\n')
        print('set_param is for change integrator parameters\n')
        print('integrator_info is for some info staff\n')     