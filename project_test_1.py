#!/usr/bin/env python3
# -*- coding: utf-8 -*-
import model
from detect_event import detection_tool
from integrator import integrator_tool as it
from stop_funcs_p import stopFunCombined

I = it(1e-14, 1e6, 'dopri5', stopFunCombined)
M = model.model_tool('Sun-Earth', I, None, None)

example = detection_tool(M)
arr = example.detect()
example.plot(arr)



 



