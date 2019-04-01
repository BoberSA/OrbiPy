#!/usr/bin/env python3
# -*- coding: utf-8 -*-

#events = {'left':[ev1], 'right':[ev2]}
#dv = Corrector.time2Sphere(Model, initial_vector, 90, events, 0.05, retit=False, maxit=100)
#initial_vector[4] = dv
#result_array = Integrator.integrate_ode(Model, initial_vector, [0, 2 * np.pi])
"""
"""
import unittest
import orbipyd as orb
import numpy as np

Integrator = orb.integrator.integrator_tool(1e-14, 1e6, 'dopri5', None)
Model = orb.model.model_tool('Sun-Earth', Integrator, None, None)
X0km = -200000
Z0km =  200000
initial_vector = np.array([Model.L2 + X0km/Model.ER, 0, Z0km/Model.ER, 0, 0, 0])
leftp = Model.mu1 + 500000 / Model.ER
rightp = Model.L2 + 500000 / Model.ER
ev1 = orb.events.event_X(leftp)
ev2 = orb.events.event_X(rightp)
Corrector = orb.correction.correction_tool()
Keeper = orb.station_keeping.station_keeping(Model, Corrector, initial_vector)

class integrator_test(unittest.TestCase):

    

        #result_array, dv, evout = Keeper.orbit_calculate(2 * np.pi, ev1, ev2)

    def test_integrate_ode(self):
        x = np.array([Model.L2 + X0km/Model.ER, 0, Z0km/Model.ER, 0, 0, 0])
        x[4] = 0.0076201674736758925
        y = np.array([1.00871674e+00,  9.75960362e-04,  1.29750112e-03,  3.40512141e-04,\
                      7.33582905e-03, -2.06372494e-04])
        result = Integrator.integrate_ode(Model, x, [0, 2 * np.pi])
        self.assertEqual(np.array(result[-1][:6]).any(), y.any())
        
    def test_brent(self):
        x = np.array([Model.L2 + X0km/Model.ER, 0, Z0km/Model.ER, 0, 0, 0])
        x[4] = 0.0076201674736758925
        y = np.array([ 1.00869484e+00, -3.60291918e-13,  1.33689840e-03, -1.32608140e-13,\
                      7.62016747e-03,  3.39982679e-13])
        t, result = Integrator.brent(Model, ev1, 0, 0.001, x)
        print(result, y)
        self.assertEqual(np.array(result).any(), y.any())
        
if __name__ == '__main__':
    unittest.main()   