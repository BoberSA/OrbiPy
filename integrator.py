#!/usr/bin/env python3
# -*- coding: utf-8 -*-
#import stop_functions as sf
import numpy as np
from scipy import integrate

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
    prop = integrate.ode(model.equation)

    if self.int_param != None:
        method = self.method
        prop.set_integrator(method, method=self.method, rtol=self.rtol)
    else:
        prop.set_integrator('dopri5')

    prop.set_initial_value(s0, tspan[0])
    prop.set_f_params(*[mu])

    lst = []
    if self.int_param['stopf']!=None:
        prop.set_solout(lambda t, s: self.int_param['stopf'](t, s, lst, events, out))
    else:
        prop.set_solout(lambda t, s: self.stopFunCombined(t, s, lst, events, out))
    prop.integrate(tspan[1])
    del prop
    if len(out) > 0:
        cor_out = self.correctEvents(model, events, out, None, sn=len(s0), int_param=self.int_param)
        out.clear()
        out.extend(cor_out)
      
    if retarr:
 
        return np.asarray(lst)
        
        
   def correctEvents(self, model, events, evout, prop, sn, int_param):
    out = []
    tol=1e-12
    maxiter = 50
    for ev in evout[1:]:
        if ev[4] == False:
            out.append([ev[0], ev[1], ev[3][:sn+1], ev[5]])
            continue

        ret = self.newton_int(model, events[ev[0]], 
                         ev[2][sn], 
                         ev[2][:sn], 
                         tol=tol,
                         maxiter=maxiter)
        if ret is None:
            print('Newton failed -> Calling Brent :)')
            t, s = self.brent(model, events[ev[0]], 
                         ev[2][sn], 
                         ev[3][sn], 
                         ev[2][:sn], 
                         tol=tol,
                         maxiter=maxiter)
        else:
            t, s = ret

        out.append([ev[0], ev[1], list(s)+[t], ev[5]])
        
    return out 

   def newton_int(self, model, event, t0, s0, mu1, int_param, tol=1e-12, maxiter=50, debug=False):
    import matplotlib.pyplot as plt
    import warnings
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
        s1 = model.integrator.integrate_ode(model, s0, [t0, t1], int_param=int_param)[-1, :6]
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
    return None
    
   def brent(model, event, t0, t1, s0, mu1, int_param, tol=1e-12, maxiter=50, debug=False):
    import scipy.optimize
    import math
    ivar = event['ivar']
    stopval = event['stopval']
    evkwargs = event.get('kwargs', {})
    
    s_opt = [0]
    
    def fopt(t, s0, t0):
        if t == t0:
            s = s0.copy()
        else:
            s = model.integrator.integrate_ode(model, s0, [t0, t], int_param=int_param)[-1, :6]
        s_opt[0] = s
        fval = ivar(t, s, **evkwargs) - stopval
        return math.fabs(fval)
        
    t_opt = scipy.optimize.brent(fopt, args=(s0, t0), brack=(t0, t1), tol=tol)

    return t_opt, s_opt[0]    
       
   def stopFunCombined(self, t, s, lst, events, out=[], **kwargs):
    if not events:
        return 0
    
    terminal = False
    cur_ivs = []
    sn = s.shape[0] + 1
    #print("that is events", events)
    for event in events:
        if(type(event)==dict):
            ivar = event['ivar']
            evkwargs = event.get('kwargs', {})        
            cur_iv = ivar(t, s, **evkwargs)
            cur_ivs.append(cur_iv)

    if not lst: # fast way to check if lst is empty
        cur_cnt = []
        for event in events:
            if(type(event)==dict):
                cur_cnt.append(event.get('count', -1))
        if not out:
            out.append(cur_cnt)
        else:
            out[0] = cur_cnt
        lst.append([*s,t,*cur_ivs])
        return 0
    lst.append([*s,t,*cur_ivs])
    cur_cnt = out[0]
    for i, event in enumerate(events):
     if(type(event)==dict):
        stopval = event.get('stopval', 0)
        direction = event.get('direction', 0)
        corr = event.get('corr', True)
        isterminal = event.get('isterminal', True)
        init_cnt = event.get('count', -1)
        
        cur_iv = cur_ivs[i]
        prev_iv = lst[-2][sn+i]

        f1 = (prev_iv < stopval) and (cur_iv > stopval) and ((direction == 1) or (direction == 0))
        f2 = (prev_iv > stopval) and (cur_iv < stopval) and ((direction == -1) or (direction == 0))
        if (f1 or f2) and ((cur_cnt[i] == -1) or (cur_cnt[i] > 0)):
            if cur_cnt[i] > 0:
                cur_cnt[i] -= 1
           
            out.append([i, # event index
                        (-1 if cur_cnt[i]==-1 else init_cnt-cur_cnt[i]), # event trigger counter
                        lst[-2].copy(), # state before event
                        lst[-1].copy(), # state after event
                        corr, # if correction is needed
                        isterminal]) # if event is terminal
            if isterminal and ((cur_cnt[i] == -1) or (cur_cnt[i] == 0)):
                terminal = True
            
    if terminal:
        return -1
    
    return 0