# -*- coding: utf-8 -*-
"""
Created on Sat Mar  7 16:55:22 2020

@author: suyan
"""

# -*- coding: utf-8 -*-
"""
Created on Tue Jul  9 11:15:32 2019
This considers a discrete time finite horizon MDP with time dependent transitions
@author: suyan
"""

import math as _math
import time as _time
import numpy as _np
import scipy.sparse as _sp
#import mdptoolbox.util as _util
import datetime

# eigendecomposition support packages
from numpy import array
from numpy.linalg import eig
from numpy import diag
from numpy import dot
from numpy.linalg import inv
#Key assumptions:
import matplotlib.pyplot as plt

#1. rewards are time homogeneous
class MDP_class(object):
    """
    parameters:
        
    tp: transition probability matrices
    dcf: discount factor
    tr: term reward
    rf: reward function
    num_abs: number of absorbing state, asummedto be the last states in tp
    
    attributes:
    tp:
    tr:
    rf:
    S:
    A:
    T:
    """
    def __init__(self,tp,tr,rf,dcf,num_abs):
        
        if dcf is not None:
            self.dcf = float(dcf)
            assert 0.0 < self.dcf <= 1.0, (
                "Discount rate must be in ]0; 1]"
            )
            if self.dcf == 1:
                print("WARNING: check conditions of convergence. With no "
                      "discount, convergence can not be assumed.")
        # Initially the time taken to perform the computations is set to None
        self.time = None
        # set the initial iteration count to zero
        self.iter = 0
        # V should be stored as a vector ie shape of (S,) or (1, S)
        self.V = None
        # policy can also be stored as a vector
        self.policy = None
        #State space
        self.Ss = tp.shape[2]-num_abs
        #Action space
        self.A = tp.shape[1]
        #Time horizon
        self.T = tp.shape[0]
        #Absorbing state
        self.num_abs=num_abs
        #Transition matrix
        self.tp=tp
        #Reward function
        self.rf=rf
        #term reward
        self.tr=tr
     
    #run backward induction
    def _back_induction(self):
        nn =datetime.datetime.now()
        # this algorithm is for adjustment
        #initialize the optimal policy and value
        self.V=_np.zeros((self.Ss+self.num_abs,self.T+1),dtype=float) 
        self.policy=_np.zeros((self.Ss,self.T),dtype=float) 
        #initialize the boundry conditions & V for absorbing states
        self.V[:,-1]=self.tr
        vs = _np.zeros((self.T,self.Ss,self.A))
        if self.num_abs>0:
            for i in range(self.num_abs):
                #if not changing the epochs, treat as oridinal MDP problem and solve use backward induction   
                self.V[self.Ss+i,:]=_np.full(self.T+1, self.tr[self.Ss+i], dtype=float)
                self.V[self.Ss+i,:]=_np.full(self.T+1, self.tr[self.Ss+i], dtype=float)
        for t in reversed(range(self.T)):
            for i in range(self.Ss):
                v_a = _np.zeros(self.A,dtype=float)
                for a in range(self.A):
                    _help = self.tp[t,a,i,:] @ self.V[:,t+1]
                    v_a[a] = self.rf[t,i,a]+self.dcf*_help
                    vs[t,i,a] = v_a[a]
                # find & record the optimal policy & values
                if abs(min(v_a) - max(v_a))<=10**(-10):
                    max_a=_np.argmax(v_a)
                    self.V[i,t]=v_a[max_a]
                    self.policy[i,t]=9
                else:
                    max_a=_np.argmax(v_a)
                    self.V[i,t]=v_a[max_a]
                    self.policy[i,t]=max_a
        
        print("Optimal value for backward induction:")
        print(self.V[0:self.Ss])
        print("Optimal policy for backward induction:")
        print(self.policy)
        print(datetime.datetime.now()-nn)
        return self.policy,self.V


            

        
            
            
            
            
            
            
            
            
            
            
            
            
            