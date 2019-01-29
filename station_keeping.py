#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Jan 17 17:46:12 2019

@author: deniszagorodnev
"""

import numpy as np


class station_keeping():
    def __init__(self, model, correction, y0):
        #type of those elements - class model and class correction
        self.model = model
        self.corr = correction
        self.y0 = y0
        
    def orbit_calculate(self, time, ev1, ev2):
        events = {'left':[ev1], 'right':[ev2]}
        event_list = events['left']+events['right']
        
        intervals = int(time/(2*np.pi))
        intervals = 7
        print(intervals)
        traectory = []
        col_dv=[]

        initial_state = self.y0
        for i in range (0, intervals):
            evout=[]

            #print ("initial_state = ", initial_state)
            #dv = self.corr.findVLimits(self.model, initial_state, 90, events, 0.05, retit=False, maxit=100)
            dv = self.corr.time2Sphere(self.model, initial_state, 90, events, 0.05, retit=False, maxit=100)
            initial_state[4] = dv
            #print ("initial_state + dv = ", initial_state)
            col_dv.append(dv)


            time_range = [time * i / intervals, time * (i + 1) / intervals]
            #print ("time_range = ", time_range)
            arr = self.model.integrator.integrate_ode(self.model, initial_state, time_range, event_list, out=evout)
            traectory.extend(arr[:-1])
            
            initial_state = arr[-1][:6] 
            
            
        #arr = self.model.integrator.integrate_ode(self.model, self.y0, [int(time//interval)*interval, time], events['left']+events['right'])
        #traectory.extend(arr)   
        
        
        return(np.array(traectory), np.array(col_dv))
        
    
    