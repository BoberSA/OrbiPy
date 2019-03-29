#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sun Jan 13 16:42:33 2019

@author: deniszagorodnev
"""

class event():
    """A class for set possible events"""
    def __init__(self, stopval=0, direction=0, isterminal=True, corr=True, count=-1):
        ''' 
        isterminal - is the event terminal or not
        '''
        self.stopval = stopval
        self.direction = direction
        self.isterminal = isterminal
        self.corr = corr
        self.count=count
        
    def ivar(self, t, s):
        ''' Independent variable function that returns coordinate 
        of spacecraft state vector.
        Should be used in stopFun and accurateEvent functions
        
        General view of event function: E(t, s) = ivar(t, s) - stopval
        
        For example, for crossing with plane with x coordinate we need to calculate 
        E(t,s) = event_X.ivar(t,s) - x
        
        Parameters
        ----------
        t : scalar
            Dimensionless time (same as angle of system rotation)        
        s : array_like with 6 components
            State vector of massless spacecraft (x,y,z,vx,vy,vz)
        
        Returns
        -------
        Depends on specific event.
        '''
        None
        
class event_X(event):
    ''' 
        Returns
        -------
        X coordinate of state vector s.
        '''
    def ivar(self, t, s):
        return(s[0])
        
class event_Y(event):
    ''' 
        
        Returns
        -------
        Y coordinate of state vector s.
        '''
    def ivar(self, t, s):
        return(s[1])
    
class event_Z(event):
    ''' 
        
        Returns
        -------
        Z coordinate of state vector s.
        '''
    def ivar(self, t, s):
        return(s[2])
    
class event_VX(event):
    ''' 
        
        Returns
        -------
        VX coordinate of state vector s.
        '''
    def ivar(self, t, s):
        return(s[3])
        
class event_VY(event):
    ''' 
        
        Returns
        -------
        VY coordinate of state vector s.
        '''
    def ivar(self, t, s):
        return(s[4])
        
class event_VZ(event):
    ''' 
        
        Returns
        -------
        VZ coordinate of state vector s.
        '''
    def ivar(self, t, s):
        return(s[5])