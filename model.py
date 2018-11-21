#!/usr/bin/env python3
# -*- coding: utf-8 -*-

#Содержимое: константы, объект integrator, правая часть уравнения, координаты точек либрации

from integrator import integrator_tool
from three_body import crtbp
import lagrange_pts_p as lpts

class model_tool:
    def __init__(self):
        
        self.model_help()
        
        print('default or custom model?')
        info = input()
        if(info == 'custom'):
            flag = 1
            while(flag==1):
                flag = self.ch_params()  
        
            
            
                
            
        self.P1_mass =  1.9891e30 
        self.P2_mass =  5.97e24 
        self.P1 = 'Sun'
        self.P2 = 'Earth'
        self.ER =  1.496e8 
        self.mu1 = self.P1_mass / (self.P1_mass + self.P2_mass)
        
        self.integrator = integrator_tool()
        self.equation = crtbp()
        self.lagrange_points = lpts.lagrange_pts(self.equation.ode, self.mu1)
    
    
    def model_info(self):
        print('P1 body is ', self.P1)
        print('P2 body is ', self.P2)
        print('Radius is ', self.ER)
    
    def ch_params(self):
        
        print('Please choose P1/P2/Integrator/CRTBP or Quit')
        name = input()
        
        if(name == 'P1'):
            P1_mass = 0
            P1 = 'P1'
            #name = input()
            self.P1_mass = P1_mass
            self.P1 = P1
            
        elif(name == 'P2'):
            #name = input()
            P2_mass = 0
            P2 = 'P2'
            #name = input()
            self.P2_mass = P2_mass
            self.P2 = P2
            
        elif(name == 'Integrator'):
            Integrator = 0
            #return(Integrator)
            
        elif(name == 'CRTBP'):
            CRTBP = 0
            #return(Integrator)
        
        elif(name == 'Quit'):
            return(-1)
        
        else:
            None
        
        return(1)
        
    
    def model_help(self):
        print('model_info is for information about current system')
        print('ch_params is for change system parameters')
        
        
        

