#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import matplotlib.pyplot as plt
import correction

class detection_tool():
    def __init__(self, model):
        self.model = model
        

    def detect(self, y0, time, ev1, ev2):
        events = {'left':ev1, 'right':ev2}
        cor = correction.correction_tool()
        #dv = cor.time2Sphere(self.model, y0, 90, events, 0.05, retit=False, maxit=100)
        dv = cor.findVLimits(self.model, y0, 90, events, 0.05, retit=False, maxit=100)
        y0[3:5] = dv
        evout= []
        arr = self.model.integrator.integrate_ode(self.model, y0, [0, time], events['left']+events['right'], evout)
        return(arr, evout)
        
    def plot(self, arr, evout=None):
        plt.figure(figsize=(10,10))
        plt.plot(arr[:,0],arr[:,1],'.-')
        plt.axis('equal')
        ev_names = ['X:0', 'alpha:120', 'Y:0', 'alpha:60']
        if evout!=None:
         for ie, _, s, _ in evout:
            plt.plot(s[0], s[1], '+k')
            plt.text(s[0], s[1], ' [%d] %s' % (ie, ev_names[ie]))
