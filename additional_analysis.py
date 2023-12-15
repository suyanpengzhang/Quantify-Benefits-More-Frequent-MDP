#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Dec 14 13:21:24 2023

@author: suyanpengzhang
"""

from numpy import linalg as LA
import pandas as pd
import math as _math
import time as _time
import numpy as _np
import scipy.sparse as _sp

#import mdptoolbox.util as _util
import liver_MDP
import warnings

omega_less = 1
omega_more = 1

o = 0.8

pre_transplant = pd.read_csv('pre_transplant_3states.csv')

post_transplant = pd.read_csv('post_transplant_3states.csv')

gamma1 = 1-pre_transplant['state1']
gamma2 = 1-pre_transplant['state2']
gamma3 = 1-pre_transplant['state3']
def compute_a(gamma,o):
    return o*(1-gamma)
def compute_b(gamma,o):
    return 1+(1-o)*(1-gamma)
def compute_c(gamma,o):
    return gamma-1
a = compute_a(gamma3[1],0.7)
b = compute_b(gamma3[1],0.7)
c = compute_c(gamma3[1],0.7)

print((-b+_np.sqrt(b**2-4*a*c))/(2*a))