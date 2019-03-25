#!/usr/bin/env python3
# -*- coding: utf-8 -*-

#Содержимое: константы, объект integrator, правая часть уравнения, координаты точек либрации

import numpy as np
from scipy import optimize
import math
from numba import compiler, types
import os

class model_tool:
    """Класс, реализующий модель исследуемой системы."""
    def __init__(self, system, integrator, crtbp=None, STM=None, filename=None): 
        """Переменные: system - идентификатор системы трех тел в формате P1P2 (Planet-1, Planet-2) integrator - объект-функция, производящая интегрирование орбиты crtbp - правая часть системы уравнений, описывающих модель (при отсутствии таковой задается предусмотренная - ode1)
        STM - если не None, то правая часть становится ode1_STM, то есть с матрицей. Но пока ode только предусмотренная.

        Предназначение: Получить набор входных данных, необходимых для сборки модели. Вынуть из файла orbi_data.csv константы, описывающие систему. Подсчитать координаты точек либрации, используя константы. Включить входные данные в информационное поле модели."""
        if filename == None:
            data = self.prepare_file(os.path.join(os.path.dirname(__file__), 'orbi_data.csv'))
        else:
            data = self.prepare_file(os.path.join(os.path.dirname(__file__), filename))
        const = self.get_orbi(data, system)
        self.mu1 = const[0]
        self.ER = const[1] 
        self.integrator = integrator 
        if crtbp == None:
           self.equation = self.ode1
           if STM!= None:
               self.equation = self.ode1_STM 
               self.equation = compiler.compile_isolated(self.model.equation, [types.double, types.double[:], types.double], return_type=types.double[:]).entry_point
        else:
            self.equation = crtbp
        self.L1=optimize.root(self.opt, 0.5, args=(self.mu1,)).x[0]
        self.L2=optimize.root(self.opt, 2.0, args=(self.mu1,)).x[0]
        self.L3=optimize.root(self.opt, -1.0, args=(self.mu1,)).x[0]
        self.L4=np.array([self.mu1-0.5, math.sqrt(3)*0.5])
        self.L5=np.array([self.mu1-0.5, -math.sqrt(3)*0.5])

            
        
        
    def opt(self, x, mu):
        y = np.zeros(6)
        y[0] = x[0]
        return self.equation(0., y, mu)[3]
 
    
    @staticmethod
    def ode1(t, s, mu):
        """Переменные: t - время, s - вектор состояния, mu - параметр отношения масс Предназначение: Система обыкновенных дифференциальных уравнений. Описывает модель задачи трех тел. Является значением по умолчанию."""
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
        
    @staticmethod
    def ode1_STM(t, s, mu):
        """Переменные: время, вектор состояния, параметр отношения масс Предназначение: Система обыкновенных дифференциальных уравнений. Описывает модель задачи трех тел. Является значением по умолчанию. Дополняет ode1 использованием матрицы STM."""
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
 
    @staticmethod
    def prepare_file(filename):
     """Функция, осуществляющая выгрузку данных из таблицы с различными моделями для CRTBP."""
        
     import csv
     a = {}
     with open(filename, 'r') as csvfile:
        next(csvfile)
        reader = csv.reader(csvfile, delimiter=';')
        for row in reader:
            a[row[0]] = [float(row[1].replace(',', '.')), float(row[2].replace(',', '.'))]
     return(a)   
    @staticmethod
    def get_orbi(data, system):
        """Переменные: data - словарь, собранный из данных filename, system - имя объекта Предназначение: Осуществляет поиск в словаре объекта с именем system."""
        return(data.get(system))
        
        
        

