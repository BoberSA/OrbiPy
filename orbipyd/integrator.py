#!/usr/bin/env python3
# -*- coding: utf-8 -*-
import numpy as np
from scipy import integrate

import sympy as sp
sp.init_printing(use_unicode=True)

  
def stopFunCombined(t, s, lst, events, out=[]):
    terminal = False
    cur_ivs = []
    sn = s.shape[0] + 1
    for event in events:
        if(type(event)!=str):
            ivar = event.ivar      
            cur_iv = ivar(t, s)
            cur_ivs.append(cur_iv)

    if not lst: # fast way to check if lst is empty
        cur_cnt = []
        for event in events:
            if(type(event)!=str):
                cur_cnt.append(event.count)
        if not out:
            out.append(cur_cnt)
        else:
            out[0] = cur_cnt
        lst.append([*s,t,*cur_ivs])
        return 0
    lst.append([*s,t,*cur_ivs])
    cur_cnt = out[0]
    
    for i, event in enumerate(events):
     if(type(events[i])!=str):
        
        stopval=event.stopval
        direction=event.direction
        corr=event.corr
        isterminal=event.isterminal
        init_cnt=event.count
        
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

class integrator_tool:
   def __init__(self, rtol, nmax, method=None, stopf=None):
       
        self.rtol = float(rtol) 
        self.nmax = float(nmax) 
        self.method = method 
        self.atol = float(rtol)
        self.stopf = stopf
   
   def integrate_ode(self, model, s0, tspan, events=[], out=[]):
    retarr=True
    mu = model.mu1
    prop = integrate.ode(model.equation)

    if self.method != None:
        method = self.method
        prop.set_integrator(method, method=self.method, rtol=self.rtol)
    else:
        prop.set_integrator('dopri5')

    prop.set_initial_value(s0, tspan[0])
    prop.set_f_params(*[mu])

    lst = []
    if self.stopf!=None:
        prop.set_solout(lambda t, s: self.int_param['stopf'](t, s, lst, events, out))
    else:
        
        prop.set_solout(lambda t, s: stopFunCombined(t, s, lst, events, out))
    prop.integrate(tspan[1])
    del prop
    if len(out) > 0:
        cor_out = self.correctEvents(model, events, out, None, sn=len(s0))
        out.clear()
        out.extend(cor_out)
      
    if retarr:
 
        return np.asarray(lst)
        
        
   def correctEvents(self, model, events, evout, prop, sn):
    out = []
    tol=1e-14
    maxiter = 50
    for ev in evout[1:]:
        if ev[4] == False:
            out.append([ev[0], ev[1], ev[3][:sn+1], ev[5]])
            continue
        t, s = self.brent(model, events[ev[0]], ev[2][sn], ev[3][sn], ev[2][:sn], tol=tol, maxiter=maxiter)
        out.append([ev[0], ev[1], list(s)+[t], ev[5]]) 
    return out 

   def brent(self, model, event, t0, t1, s0, tol=1e-12, maxiter=50, debug=False):
    import scipy.optimize
    import math
    ivar = event.ivar
    stopval = event.stopval 
    s_opt = [0]
    
    def fopt(t, s0, t0):
        if t == t0:
            s = s0.copy()
        else:
            s = model.integrator.integrate_ode(model, s0, [t0, t])[-1, :6]
        s_opt[0] = s
        fval = ivar(t, s) - stopval
        return math.fabs(fval)
        
    t_opt = scipy.optimize.brent(fopt, args=(s0, t0), brack=(t0, t1), tol=tol)

    return t_opt, s_opt[0]  

  
   