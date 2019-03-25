#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sun Jan 13 16:42:33 2019

@author: deniszagorodnev
"""

class event():
    """Класс, реализующий возможные события."""
    def __init__(self, stopval=0, direction=0, isterminal=True, corr=True, count=-1):
        ''' 
        isterminal - является ли событие терминальным
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
    
        Parameters
        ----------
        t : scalar
            Dimensionless time (same as angle of system rotation)        
        s : array_like with 6 components
            State vector of massless spacecraft (x,y,z,vx,vy,vz)
        
        Returns
        -------
        coordinate of state vector s.
        '''
        None
        
class event_X(event):
    ''' Independent variable function that returns X coordinate 
        of spacecraft state vector.
        Should be used in stopFun and accurateEvent functions
    
        Parameters
        ----------
        t : scalar
            Dimensionless time (same as angle of system rotation)        
        s : array_like with 6 components
            State vector of massless spacecraft (x,y,z,vx,vy,vz)
        
        Returns
        -------
        X coordinate of state vector s.
        '''
    def ivar(self, t, s):
        return(s[0])
        
class event_Y(event):
    ''' Independent variable function that returns Y coordinate 
        of spacecraft state vector.
        Should be used in stopFun and accurateEvent functions
    
        Parameters
        ----------
        t : scalar
            Dimensionless time (same as angle of system rotation)        
        s : array_like with 6 components
            State vector of massless spacecraft (x,y,z,vx,vy,vz)
        
        Returns
        -------
        Y coordinate of state vector s.
        '''
    def ivar(self, t, s):
        return(s[1])
    
class event_Z(event):
    ''' Independent variable function that returns Z coordinate 
        of spacecraft state vector.
        Should be used in stopFun and accurateEvent functions
    
        Parameters
        ----------
        t : scalar
            Dimensionless time (same as angle of system rotation)        
        s : array_like with 6 components
            State vector of massless spacecraft (x,y,z,vx,vy,vz)
        
        Returns
        -------
        Z coordinate of state vector s.
        '''
    def ivar(self, t, s):
        return(s[2])
    
class event_VX(event):
    ''' Independent variable function that returns VX coordinate 
        of spacecraft state vector.
        Should be used in stopFun and accurateEvent functions
    
        Parameters
        ----------
        t : scalar
            Dimensionless time (same as angle of system rotation)        
        s : array_like with 6 components
            State vector of massless spacecraft (x,y,z,vx,vy,vz)
        
        Returns
        -------
        VX coordinate of state vector s.
        '''
    def ivar(self, t, s):
        return(s[3])
        
class event_VY(event):
    ''' Independent variable function that returns VY coordinate 
        of spacecraft state vector.
        Should be used in stopFun and accurateEvent functions
    
        Parameters
        ----------
        t : scalar
            Dimensionless time (same as angle of system rotation)        
        s : array_like with 6 components
            State vector of massless spacecraft (x,y,z,vx,vy,vz)
        
        Returns
        -------
        VY coordinate of state vector s.
        '''
    def ivar(self, t, s):
        return(s[4])
        
class event_VZ(event):
    ''' Independent variable function that returns VZ coordinate 
        of spacecraft state vector.
        Should be used in stopFun and accurateEvent functions
    
        Parameters
        ----------
        t : scalar
            Dimensionless time (same as angle of system rotation)        
        s : array_like with 6 components
            State vector of massless spacecraft (x,y,z,vx,vy,vz)
        
        Returns
        -------
        VZ coordinate of state vector s.
        '''
    def ivar(self, t, s):
        return(s[5])