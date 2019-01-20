#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Jan 17 17:46:12 2019

@author: deniszagorodnev
"""


class station_keeping():
    def __init__(self, model, correction, y0):
        #type of those elements - class model and class correction
        self.model = model
        self.corr = correction
        self.y0 = y0
        
    def orbit_calculate(self, time, ev1, ev2):
        import numpy as np
        events = {'left':[ev1], 'right':[ev2]}
        interval = 0.01
        traectory = []
        col_dv=[]
        for i in range (0, time//interval-1):
            arr = self.model.integrator.integrate_ode(self.model, self.y0, [i*interval, i*interval+1], events['left']+events['right'])
            traectory.append(arr)
            dv = self.correction.findVLimits(self.model, self.y0, 90, events['left']+events['right'], 0.05, retit=False, maxit=100)
            self.y0[3:5] = dv
            col_dv.append(dv)
            
            
        arr = self.model.integrator.integrate_ode(self.model, self.y0, [time//interval*interval, time], events['left']+events['right'])
        traectory.append(arr)   
            
        return(np.array(traectory), np.array(dv))
        
    
    