#!/usr/bin/env python3
# -*- coding: utf-8 -*-
#import stop_functions as sf
import numpy as np
import scipy
from stop_funcs_p import stopNull,correctEvents

import sympy as sp
sp.init_printing(use_unicode=True)

class integrator_tool:
   def __init__(self, rtol, nmax, method):
       
        #self.integrator_help()
        
        
            #self.rtol = 1e-14 # integration relative tolerance
            #self.nmax = 1e6 # max number of integration steps
            #self.method = 'dop853'
            
        self.rtol = float(rtol) 
        self.nmax = float(nmax) 
        self.method = method   
        #default case
        
        self.int_param = {'atol':self.rtol, 'rtol':self.rtol, 'nsteps':self.nmax, 'method':self.method}
   
   def integrate_ode(self, model, s0, tspan, retarr=True):
    right_part = model.equation
    mu = model.mu1
    prop = scipy.integrate.ode(right_part)

    if self.int_param != None:
        method = self.method
        prop.set_integrator(method, method=self.method, rtol=self.rtol)
    else:
        prop.set_integrator('dopri5')

    prop.set_initial_value(s0, tspan[0])
    prop.set_f_params(*[mu])

    lst = []
    #print(kwargs)
    if 'stopf' in self.int_param:
        prop.set_solout(lambda t, s: self.int_param['stopf'](t, s, lst, self.int_param))
    else:
        prop.set_solout(lambda t, s: stopNull(t, s, lst))
    prop.integrate(tspan[1])

    del prop
    
    #evout = kwargs.get('out', [])
    #if len(evout) > 0:
     #   events = kwargs.get('events', [])
#        cor_out = correctEvents(events, evout, prop, sn=len(s0),
#                            tol=kwargs['int_param']['atol'])
      #  cor_out = correctEvents(events, evout, None, sn=len(s0),
       #                     tol=kwargs['int_param']['atol'], 
        #                    int_param=kwargs['int_param'],
         #                   mu1=kwargs['mu'])
       # evout.clear()
        #evout.extend(cor_out)
    
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