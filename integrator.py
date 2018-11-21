#!/usr/bin/env python3
# -*- coding: utf-8 -*-
#import stop_functions as sf
import numpy as np
import scipy
from stop_funcs_p import stopNull,correctEvents

import sympy as sp
sp.init_printing(use_unicode=True)

class integrator_tool:
   def __init__(self):
        rtol = 1e-14 # integration relative tolerance
        nmax = 1e6 # max number of integration steps
        self.int_param = {'atol':rtol, 'rtol':rtol, 'nsteps':nmax, 'method':'dop853'}
   @staticmethod 
   def integrate_ode_old(right_part, model, s0, tspan, retarr=True, **kwargs):
        #stop = sf.stop_tool()
        crtbp = right_part
        mu1 = model.mu1
       # crtbp = compiler.compile_isolated(crtbp, [types.double, types.double[:], types.double], return_type=types.double[:]).entry_point
        prop = scipy.integrate.ode(crtbp)

        if 'int_param' in kwargs:
            method = kwargs['int_param'].get('method', 'dopri5')
            prop.set_integrator(method, **kwargs['int_param'])
        else:
            prop.set_integrator('dopri5')
    
        prop.set_initial_value(s0, tspan[0])
        prop.set_f_params(*[mu1])
        i = model.integrator.integrate_ode(right_part, model, s0, tspan, retarr=True, **kwargs)
        lst = []
        kwargs['mu'] = mu1
        #print(kwargs)
        if 'stopf' in kwargs:
            prop.set_solout(lambda t, s: kwargs['stopf'](t, s, lst, **kwargs))
        else:
            prop.set_solout(lambda t, s: stopNull(t, s, lst))
        prop.integrate(tspan[1])
    
        del prop
        
        evout = kwargs.get('out', [])
        if len(evout) > 1:
            events = kwargs.get('events', [])
            
            cor_out = correctEvents(i, events, evout, None, sn=len(s0),
                                tol=kwargs['int_param']['atol'], 
                                int_param=kwargs['int_param'],
                                mu1=kwargs['mu'])
            evout.clear()
            evout.extend(cor_out)
        
        if retarr:
            return np.asarray(lst)
    
   @staticmethod 
   def integrate_ode(right_part, model, s0, tspan, retarr=True, **kwargs):
  
    mu = model.mu1
    prop = scipy.integrate.ode(right_part)

    if 'int_param' in kwargs:
        method = kwargs['int_param'].get('method', 'dopri5')
        prop.set_integrator(method, **kwargs['int_param'])
    else:
        prop.set_integrator('dopri5')

    prop.set_initial_value(s0, tspan[0])
    prop.set_f_params(*[mu])

    lst = []
    kwargs['mu'] = mu
    #print(kwargs)
    if 'stopf' in kwargs:
        prop.set_solout(lambda t, s: kwargs['stopf'](t, s, lst, **kwargs))
    else:
        prop.set_solout(lambda t, s: stopNull(t, s, lst))
    prop.integrate(tspan[1])

    del prop
    
    evout = kwargs.get('out', [])
    if len(evout) > 0:
        events = kwargs.get('events', [])
#        cor_out = correctEvents(events, evout, prop, sn=len(s0),
#                            tol=kwargs['int_param']['atol'])
        cor_out = correctEvents(events, evout, None, sn=len(s0),
                            tol=kwargs['int_param']['atol'], 
                            int_param=kwargs['int_param'],
                            mu1=kwargs['mu'])
        evout.clear()
        evout.extend(cor_out)
    
    if retarr:
        return np.asarray(lst)
