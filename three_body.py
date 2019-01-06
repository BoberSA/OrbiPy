#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import numpy as np

class crtbp:
    def __init__(self, ode, STM):
        if ode == None:
           self.ode = self.ode1
           if STM!= None:
               self.ode = self.ode1_STM
        else:
            self.ode = ode

    def ode1(self, t, s, mu):
        x, y, z, x1, y1, z1 = s
        mu2 = 1 - mu

        yz2 = y * y + z * z;
        r13 = ((x + mu2) * (x + mu2) + yz2) ** 1.5;
        r23 = ((x - mu ) * (x - mu ) + yz2) ** 1.5;

        yzcmn = (mu / r13 + mu2 / r23);

        dx1dt = 2 * y1 + x - (mu * (x + mu2) / r13 + mu2 * (x - mu) / r23);
        dy1dt = -2 * x1 + y - yzcmn * y;
        dz1dt = - yzcmn * z;
    
        ds = np.array([x1, y1, z1, dx1dt, dy1dt, dz1dt])
        return(ds)
    
    def ode1_STM(self, t, s, mu):
        x, y, z, x1, y1, z1 = s[:6]

        STM = np.ascontiguousarray(s[6:]).reshape(6, 6)
    
        #mu - примерно равно единице
        mu2 = 1 - mu
    
        yz2 = y * y + z * z;
        r1 = ((x + mu2) * (x + mu2) + yz2) ** 1.5
        r2 = ((x - mu ) * (x - mu ) + yz2) ** 1.5

        yzcmn = (mu / r1 + mu2 / r2);
    
        dx1dt =  2 * y1 + x - (mu * (x + mu2) / r1 + mu2 * (x - mu) / r2)
        dy1dt = -2 * x1 + y - yzcmn * y
        dz1dt =             - yzcmn * z
    
        #ds = np.array([x1, y1, z1, dx1dt, dy1dt, dz1dt])
        ds = np.empty_like(s)

        Uxx = 3.0*mu*(mu2 + x)**2*(y**2 + z**2 + (mu2 + x)**2)**(-2.5) - 1.0*mu*(y**2 + z**2 + (mu2 + x)**2)**(-1.5) + 3.0*mu2*(mu - x)**2*(y**2 + z**2 + (mu - x)**2)**(-2.5) - 1.0*mu2*(y**2 + z**2 + (mu - x)**2)**(-1.5) + 1.0
        Uxy = 3.0*y*(mu*(mu2 + x)*(y**2 + z**2 + (mu2 + x)**2)**(-2.5) - mu2*(mu - x)*(y**2 + z**2 + (mu - x)**2)**(-2.5))

        Uxz = 3.0*z*(mu*(mu2 + x)*(y**2 + z**2 + (mu2 + x)**2)**(-2.5) - mu2*(mu - x)*(y**2 + z**2 + (mu - x)**2)**(-2.5))
        Uyy = 3.0*mu*y**2*(y**2 + z**2 + (mu2 + x)**2)**(-2.5) - 1.0*mu*(y**2 + z**2 + (mu2 + x)**2)**(-1.5) + 3.0*mu2*y**2*(y**2 + z**2 + (mu - x)**2)**(-2.5) - 1.0*mu2*(y**2 + z**2 + (mu - x)**2)**(-1.5) + 1.0
        Uyz = 3.0*y*z*(mu*(y**2 + z**2 + (mu2 + x)**2)**(-2.5) + mu2*(y**2 + z**2 + (mu - x)**2)**(-2.5))
        Uzz = 3.0*mu*z**2*(y**2 + z**2 + (mu2 + x)**2)**(-2.5) - 1.0*mu*(y**2 + z**2 + (mu2 + x)**2)**(-1.5) + 3.0*mu2*z**2*(y**2 + z**2 + (mu - x)**2)**(-2.5) - 1.0*mu2*(y**2 + z**2 + (mu - x)**2)**(-1.5)
    
        A=np.array(((0,0,0,1,0,0),(0,0,0,0,1,0),(0,0,0,0,0,1),(Uxx,Uxy,Uxz,0,2,0),(Uxy,Uyy,Uyz,-2,0,0),(Uxz,Uyz,Uzz,0,0,0)))

        STM = np.dot(A, STM)
        ds[0], ds[1], ds[2], ds[3], ds[4], ds[5] = x1, y1, z1, dx1dt, dy1dt, dz1dt
        ds[6:] = STM.ravel()
        return (ds)
   