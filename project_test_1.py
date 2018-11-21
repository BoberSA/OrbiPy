#!/usr/bin/env python3
# -*- coding: utf-8 -*-
import model
from detect_event import detection_tool


M = model.model_tool()
M.model_info()
example = detection_tool()
example.plot(M)


#detection_tool().plot()


 



