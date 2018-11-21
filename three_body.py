#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import numpy as np

class crtbp:
    @staticmethod
    def ode(t, s, mu):
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
    
   