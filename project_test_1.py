#!/usr/bin/env python3
# -*- coding: utf-8 -*-
import model
from detect_event import detection_tool
from integrator import integrator_tool 
import events as evs
import numpy as np
import station_keeping as sk
import correction
import matplotlib.pyplot as plt

I = integrator_tool(1e-14, 1e6, 'dopri5', None)
M = model.model_tool('Sun-Earth', I, None, None)
X0km = -200000
Z0km =  200000
y0 = np.array([M.L2 + X0km/M.ER, 0, Z0km/M.ER, 0, 0, 0])
leftp = M.mu1 + 500000 / M.ER
rightp = M.L2 + 500000 / M.ER
ev1 = evs.event_X(leftp)
ev2 = evs.event_X(rightp)
#example = detection_tool(M)
#arr, evout = example.detect(y0, 2*np.pi, [ev1], [ev2])
#example.plot(arr, evout)

cor = correction.correction_tool()
K = sk.station_keeping(M, cor, y0)
arr, dv = K.orbit_calculate(8 * np.pi, ev1, ev2)
#print(arr)
plt.figure(figsize=(10,10))
plt.plot(arr[:,0],arr[:,1],'.-')
plt.axis('equal')





