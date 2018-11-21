#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Nov 15 09:43:54 2018

@author: deniszagorodnev
"""
import numpy as np
import matplotlib.pyplot as plt
import warnings
import math
from model import model_tool


class stop_tool(model_tool):
   def newton_int(self, event, t0, s0, mu1, int_param, tol=1e-12, maxiter=50, debug=False):
    
    
    ivar = event['ivar']
    dvar = event['dvar']
    stopval = event['stopval']
    evkwargs = event.get('kwargs', {})
    
    t = t0
    s = s0.copy()

    if debug:
        lst = []
    
    for iter in range(maxiter):
        fder = dvar(t, s, **evkwargs)
        if fder == 0:
            msg = "derivative was zero."
            warnings.warn(msg, RuntimeWarning)
            return t, s
        fval = ivar(t, s, **evkwargs) - stopval
        if debug:
            lst.append((t, fval, fder))

        if abs(fval) < tol:
            return t, s

        newton_step = fval / fder
        t1 = t - newton_step
        s1 = self.integrator.integrate_ode(mu1, s0, [t0, t1], int_param=int_param)[-1, :6]
#        prop.set_initial_value(s0, t0)
#        s1 = prop.integrate(t1)

        t = t1
        s = s1

    if debug:
        lst = np.array(lst)
        fig, ax = plt.subplots(1, 3, figsize=(15,5))
        ax[0].plot(range(len(lst)), lst[:,0])
        ax[0].set_title('steps')
        ax[1].plot(lst[:,0], lst[:,1], 'r')
        ax[1].set_title('fval')
        ax[2].plot(lst[:,0], lst[:,2], 'g')
        ax[2].set_title('fder')
        fig.tight_layout()

#    print('\nt0:', t0, '\ns0:', s0)
#    msg = "Failed to converge after %d iterations, value is %s" % (maxiter, t)
#    raise RuntimeError(msg) 
    return None

   def brent(self, event, t0, t1, s0, mu1, int_param, tol=1e-12, maxiter=50, debug=False):
    import scipy.optimize
    """
    Special edition of Brent method for event-finding process.
    Slower than Newton method but converges better. This method 
    used when Newton's method can't converge.
    """
    
    ivar = event['ivar']
#    dvar = event['dvar']
    stopval = event['stopval']
    evkwargs = event.get('kwargs', {})
    
    s_opt = [0]
    
    def fopt(t, s0, t0):
        if t == t0:
            s = s0.copy()
        else:
            s = self.integrator.integrate_ode(mu1, s0, [t0, t], int_param=int_param)[-1, :6]
        s_opt[0] = s
        fval = ivar(t, s, **evkwargs) - stopval
        return math.fabs(fval)
        
    t_opt = scipy.optimize.brent(fopt, args=(s0, t0), brack=(t0, t1), tol=tol)

    return t_opt, s_opt[0]

   def correctEvents(self, events, evout, prop, tol=1e-12, sn=6, maxiter=50, **kwargs):

        out = []
        for ev in evout[1:]:
            if ev[4] == False:
                out.append([ev[0], ev[1], ev[3][:sn+1], ev[5]])
                continue

            ret = self.newton_int(events[ev[0]], 
                             ev[2][sn], 
                             ev[2][:sn], 
                             tol=tol,
                             maxiter=maxiter,
                             **kwargs)
            if ret is None:
                print('Newton failed -> Calling Brent :)')
                t, s = self.brent(events[ev[0]], 
                             ev[2][sn], 
                             ev[3][sn], 
                             ev[2][:sn], 
                             tol=tol,
                             maxiter=maxiter,
                             **kwargs)
            else:
                t, s = ret

            out.append([ev[0], ev[1], list(s)+[t], ev[5]])
        
        return out
        
   def stopNull(t, s, lst, **kwargs):
        lst.append(np.hstack((s,t)))
        return 0
    
  