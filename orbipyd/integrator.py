#!/usr/bin/env python3
# -*- coding: utf-8 -*-
import numpy as np
from scipy import integrate

import sympy as sp
sp.init_printing(use_unicode=True)

  
def stopFunCombined(t, s, lst, events, out=[]):
    """ Universal event detection function that handles multiple events. 
        Intended for scipy.integrate.ode solout application. Provides 
        termination of integration process when first terminate event occur. 
        This happens when independent variable associated with this event 
        goes through defined stopval value in specified direction.
        Uses almost the same ideas as in matlab event functions.
        Can be used for gathering all intergation steps.
        Shoudn't be called directly but through scipy.integrate.ode.
    
    Parameters
    ----------
    t : scalar
        Dimensionless time (same as angle of system rotation)
        
    s : array_like with 6 components
        State vector of massless spacecraft (x,y,z,vx,vy,vz)
    lst : list
        Every call of this function put [*s,t,*cur_ivs] into lst, where 
        t - time (at current integration step),
        s - spacecraft state vector at time t,
        cur_ivs - list of independent variable values at time t.
        
    events : list of dicts
        Each dict consists of necessary information for event:
        {
        
        ivar : function(t, s, **kwargs)
            Should return independent variable from spacecraft state vector.
        
        stopval : double
            stopFun return -1 if independent variable crosses stopval value in
            right direction
        
        direction : integer
            1 : stops integration when independent variable crosses stopval value
                 from NEGATIVE to POSITIVE values
            -1 : stops integration when independent variable crosses stopval value
                 from POSITIVE to NEGATIVE values
            0 : in both cases (like 'direction' argument in matlab's event functions)
            isterminal : integer, bool         
            Terminal event terminates integration process when event occurs.
            
        corr : bool
            Determines whether it is necessary to adjust last state vector or not
        
        count : int
            Number of event occasions.
            If count == 0
            Then event doesnt occur (turned off event).
            
            If count == -1
            Then (possibly) unlimited number of events can occur.
            
            If isterminal == True
                If count == 1
                Then only one terminal event (possibly) occur.
                If count > 1
                Then non-terminal event (possibly) triggers count-1
                times and one terminal event (possibly) occur.
            If isterminal == False
                Event (possibly) occurs count times.
        kwargs : dict
            Other parameters for ivar function
        }
        
    out : list
        If non-terminal event(s) occur in [ti-1, ti] interval, 'out' will
        be filled with [ei, ci, np.array([*s,te,*cur_ivs])], where:
        ei - event index in events list,
        ci - event triggered ci times,
        te - time of event,
        s - state vector at te,
        cur_ivs - values of independent variables at te.
        
        Returns
    -------
    
    -1 : scalar
        When there are terminal event in event list and independent variable \
        assiciated with this event goes through defined stopval value in \
        specified direction. Will be treated by scipy.integrate.ode as it \
        should stop integration process.
        
    0 : scalar
        Otherwise. Will be treated by scipy.integrate.ode as it\
        should continue integration process.
    """
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
   """Класс, реализующий интегратор."""
   def __init__(self, rtol, nmax, method=None, stopf=None):
        """Переменные: Параметры для интегратора: rtol - точность, nmax - максимальное число шагов интегрирования, method - метод интегрирования Предназначение: Задать параметры интегратора int_param для дальнейшей работы с scypy.integrate.   """
       
        self.rtol = float(rtol) 
        self.nmax = float(nmax) 
        self.method = method 
        self.atol = float(rtol)
        self.stopf = stopf
   
   def integrate_ode(self, model, s0, tspan, events=[], out=[]):
    """Переменные: model - собранная модель, s0 - начальный вектор состояния, tspan - время Предназначение: интегрирует вектор состояния на время, возвращает массив векторов.   """
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
    """Calculate corrected event states using newton method with
        same tolerance as integrator used.
        
    Parameters
    ----------
    
    events : list of events
        See stopFunCombined for description
        
    evout : list (generated by stopFunCombined)
        evout elements are lists where (by index):
            0 - event index
            1 - event trigger count
            2 - state before event
            3 - state after event
            4 - correction flag
            5 - terminal flag
            
    sn : integer
        State vector length
        
        
    model : class model
        
      
    """
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
    """Функция, реализующая численный метод Брента """
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

  
   