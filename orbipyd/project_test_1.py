#!/usr/bin/env python3
# -*- coding: utf-8 -*-

'''
    Исследование возможностей модуля для расчета орбиты космического аппарата
    в ограниченной круговой задаче трех тел.
'''


"""
Необходимые импорты: из установленного модуля импортируем все его классы.
"""
import orbipyd as orb
import numpy as np
import matplotlib.pyplot as plt

"""
Задаем интегратор и модель системы. В нашем примере это модель Солнце-Земля. Все необходимые константы 
передаются в качестве аргументов или получаются встроенными функциями по умолчанию.
"""
Integrator = orb.integrator.integrator_tool(1e-14, 1e6, 'dopri5', None)
Model = orb.model.model_tool('Sun-Earth', Integrator, None, None)

"""
Маленькая гало орбита, координаты по осям X и Z
"""
X0km = -200000
Z0km =  200000

"""
Начальный вектор состояния космического аппарата. В его определении 
использованы координаты точек либрации для заданной системы, расчитанные автоматически 
при инициализации модели.
"""
initial_vector = np.array([Model.L2 + X0km/Model.ER, 0, Z0km/Model.ER, 0, 0, 0])

"""
Расположение соответсвующих ограничивающих плоскостей.
"""
leftp = Model.mu1 + 500000 / Model.ER
rightp = Model.L2 + 500000 / Model.ER

"""
События.
"""
ev1 = orb.events.event_X(leftp)
ev2 = orb.events.event_X(rightp)

"""
Задаем инструмент для коррекции орбиты.
Далее задаем инструмент для расчета орбиты на много оборотов вперед.
"""
Corrector = orb.correction.correction_tool()
#Keeper = station_keeping.station_keeping(Model, Corrector, initial_vector)
#result_array, dv, evout = Keeper.orbit_calculate(8 * np.pi, ev1, ev2)
events = {'left':[ev1], 'right':[ev2]}
dv = Corrector.corrector(Model, initial_vector, 90, events, 0.05, retit=False, maxit=100)
print(initial_vector, dv)
initial_vector[3] = dv[0]
initial_vector[4] = dv[1]
result_array = Integrator.integrate_ode(Model, initial_vector, [0, 2 * np.pi])
#print(result_array[0], result_array[-1])
"""
"""
#ev_names = ['X:0', 'alpha:120', 'Y:0', 'alpha:60']
plt.figure(figsize=(10,10))
plt.plot(result_array[:,0],result_array[:,1],'.-')


#print(initial_vector, result_array[-1])
#for ie, _, s, _ in evout:
#    plt.plot(s[0], s[1], '+k')
#    plt.text(s[0], s[1], ' [%d] %s' % (ie, ev_names[ie]))    
plt.axis('equal')

plt.show()



