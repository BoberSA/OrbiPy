#!/usr/bin/env python3
# -*- coding: utf-8 -*-
import model
from detect_event import detection_tool
from integrator import integrator_tool as it
import events as evs
import numpy as np

#integrator initial
I = it(1e-14, 1e6, 'dopri5', None)

#model initial
M = model.model_tool('Sun-Earth', I, None, None)

#take lagrange pts from model
L = M.lagrange_points       
L1 = L[0, 0]
L2 = L[1, 0]

#setting initial y0 vector
X0km = -200000
Z0km =  200000
y0 = np.array([L2 + X0km/M.ER, 0, Z0km/M.ER, 0, 0, 0])
leftp = M.mu1 + 500000 / M.ER
rightp = L2 + 500000 / M.ER

#events initial
ev1 = evs.event_X(leftp, 0, True, False).event
ev2 = evs.event_X(rightp, 0, True, False).event

#starting event detection procedure
example = detection_tool(M)
arr, evout = example.detect(y0, ev1, ev2)

#events and orbit visuslisation
example.plot(arr, evout)



 



