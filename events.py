#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sun Jan 13 16:42:33 2019

@author: deniszagorodnev
"""

class event():
    def __init__(self, stopval=0, direction=0, isterminal=True, corr=True, count=-1):
        self.stopval = stopval
        self.direction = direction
        self.isterminal = isterminal
        self.corr = corr
        self.count=count
        
    def ivar(self, t, s):
        None
        
class event_X(event):
    def ivar(self, t, s):
        return(s[0])
        
class event_Y(event):
    def ivar(self, t, s):
        return(s[1])
    
class event_Z(event):
    def ivar(self, t, s):
        return(s[2])
    
class event_VX(event):
    def ivar(self, t, s):
        return(s[3])
        
class event_VY(event):
    def ivar(self, t, s):
        return(s[4])
        
class event_VZ(event):
    def ivar(self, t, s):
        return(s[5])