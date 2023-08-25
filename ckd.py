#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue May  2 15:57:35 2023

@author: suyanpengzhang
"""

from numpy import linalg as LA
from scipy.linalg import eig, sqrtm
import pandas as pd
import math as _math
import time as _time
import numpy as _np
import scipy.sparse as _sp
import pickle
#import mdptoolbox.util as _util
import liver_MDP
import warnings


# Load the array from the pickle file
with open('stage_reward.pickle', 'rb') as f:
    stage_reward = pickle.load(f)
with open('terminal_reward.pickle', 'rb') as f:
    terminal_reward = pickle.load(f)
with open('transitions.pickle', 'rb') as f:
    transitions = pickle.load(f)
terminal_reward = stage_reward[-1,:,1].copy()
lambda_ = 0.97
num_absorbing = 4
T = 20
stage_reward = stage_reward[0:T,:,:]
transitions = transitions[0:T,:,:,:]


#make adjustment to make 3a and 3b absorbing states
for t in range(T):
    stage_reward[t,2:4,:]=0
    #stage_reward[t,1,:]=stage_reward[t,0,:]
    transitions[t,0,2,2] = 1
    transitions[t,0,2,3] = 0
    transitions[t,0,2,5] = 0
    transitions[t,0,3,2] = 1
    transitions[t,0,3,5] = 0
    transitions[t,0,3,3] = 0
    #transitions[t,0] = transitions[t,0]@transitions[t,0]
    #stage_reward[t,:,0] = stage_reward[t,:,0] * 2


lambda_two = 0.97**(1/2)



transitions_two = _np.zeros((2*T,2,6,6))
for t in range(T):
    #eigenvalues, eigenvectors= LA.eig(transitions[t,0])
    #B = _np.dot(_np.dot(eigenvectors, _np.diag(_np.sqrt(eigenvalues))), _np.linalg.inv(eigenvectors))
    # Compute the square root of A using sqrtm
    new = _np.real_if_close(sqrtm(transitions[t,0]))
# =============================================================================
#     new[0,2]=0
#     new[0,3]=0
#     new[1,3]=0
#     new[0] = new[0]/sum(new[0])
#     new[1] = new[1]/sum(new[1])
# =============================================================================
    transitions_two[2*t,0] = new
    transitions_two[2*t+1,0] = new
    transitions_two[2*t,1] = transitions[t,1]
    transitions_two[2*t+1,1] = transitions[t,1]
terminal_reward = stage_reward[-1,:,1].copy()
terminal_reward_two = terminal_reward

stage_reward_two = _np.zeros((2*T,6,2))
ck =20
for t in range(T):
    stage_reward_two[2*t,:,1] = stage_reward[t,:,1]
    stage_reward_two[2*t+1,:,1] = stage_reward[t,:,1]
    stage_reward_two[2*t,:,0] = stage_reward[t,:,0]*0.5
    stage_reward_two[2*t+1,:,0] = stage_reward[t,:,0]*0.5
    stage_reward[t,:,0] = stage_reward_two[2*t,:,0] + lambda_two*transitions_two[2*t,0]@stage_reward_two[2*t+1,:,0]
    stage_reward_two[2*t,:,:] -=ck
    stage_reward_two[2*t+1,:,:] -=ck
    
    

    
ckd_yearly=liver_MDP.MDP_class(transitions,terminal_reward,stage_reward,lambda_,num_absorbing)
pol1,V1 = ckd_yearly._back_induction()
ckd_biyearly=liver_MDP.MDP_class(transitions_two,terminal_reward_two,stage_reward_two,lambda_two,num_absorbing)
pol2,V2 = ckd_biyearly._back_induction()    

print('yearly CKD 1 start in ',pol1[0].shape[0]-sum(pol1[0]))
print('yearly CKD 1 start in ',pol1[1].shape[0]-sum(pol1[1]))

print('twice yearly CKD 1 start in ',(pol2[0].shape[0]-sum(pol2[0]))/2)
print('twice yearly CKD 1 start in ',(pol2[1].shape[0]-sum(pol2[1]))/2)   

diff = _np.zeros((2,T))
for t in range(T):
    diff[0,t] = V2[0,2*t]-V1[0,t]
    diff[1,t] = V2[1,2*t]-V1[1,t]
print(diff)
import matplotlib.pyplot as plt
#diff = diff[:,0:7]
plt.plot(diff[0],color='black')
plt.show()
    