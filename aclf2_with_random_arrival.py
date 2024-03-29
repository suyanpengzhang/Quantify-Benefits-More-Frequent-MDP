
# -*- coding: utf-8 -*-
"""
Created on Sun Apr  5 11:35:02 2020

@author: Suyanpeng Zhang
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
warnings.filterwarnings("ignore")

pre_transplant = pd.read_csv('pre_transplant_3states.csv')

post_transplant = pd.read_csv('post_transplant_3states.csv')
c_k = 2000
print(post_transplant)

##########################################################
def compute_a(gamma,o,lambda_):
    return o*(1-gamma)*lambda_**2
def compute_b(gamma,o,lambda_):
    return 1*lambda_+(1-o)*(1-gamma)*lambda_**2
def compute_c(gamma,o,lambda_):
    return (gamma-1)*lambda_**2



def mdp(omega,alpha,beta,rr,small,test):
    gamma1 = 1-pre_transplant['state1']
    gamma2 = 1-pre_transplant['state2']
    gamma3 = 1-pre_transplant['state3']
    omega1 = omega[0]#0.9
    omega2 = omega[1]#0.95
    omega3 = omega[2]
    rw = _np.zeros([14*2,11,2])
    new_rw = _np.zeros([14,11,2])
    new_tp = _np.zeros([14,2,11,11])
    tp = _np.zeros([14*2,2,11,11])
    _lambda=0.97**(1/365)
    new_lambda = _lambda*_lambda
    utl_post = _np.array([0.576/12,0.576/12,0.576/12,0.601/12,0.601/12,0.601/12,0.601/12,0.601/12,0.601/12,0.626/12,0.626/12,0.626/12])
    num_absorbing=2
    term_rw=_np.array([0,0,0,0,0,0,0,0,0,0,0],dtype=float)
    for tt in range(14):
        t = 2*tt
        if test ==0:
            omega1 = omega[0]
            omega2 = omega[1]
            omega3 = omega[2]
        else:
            a = compute_a(gamma1[t],alpha[t],_lambda)
            b = compute_b(gamma1[t],alpha[t],_lambda)
            c = compute_c(gamma1[t],alpha[t],_lambda)
            omega1 = (-b+_np.sqrt(b**2-4*a*c))/(2*a)
            a = compute_a(gamma2[t],alpha[t],_lambda)
            b = compute_b(gamma2[t],alpha[t],_lambda)
            c = compute_c(gamma2[t],alpha[t],_lambda)
            omega2 = (-b+_np.sqrt(b**2-4*a*c))/(2*a)
            a = compute_a(gamma3[t],alpha[t],_lambda)
            b = compute_b(gamma3[t],alpha[t],_lambda)
            c = compute_c(gamma3[t],alpha[t],_lambda)
            omega3 = (-b+_np.sqrt(b**2-4*a*c))/(2*a)
        #state 1: ACLF2
        survival1 = (pre_transplant['state1'][t]+pre_transplant['state1'][t])/2
        tp[t,0,0,0] = alpha[t]*(survival1)*omega1
        tp[t,0,0,1] = (1-alpha[t])*(survival1)*omega1
        tp[t,0,0,2] = (survival1)*(1-omega1)
        tp[t,0,0,10] = (1-survival1)
        tp[t,0,1,0] = alpha[t]*(survival1)*omega1
        tp[t,0,1,1] = (1-alpha[t])*(survival1)*omega1
        tp[t,0,1,2] = (survival1)*(1-omega1)
        tp[t,0,1,10] = (1-survival1)
        tp[t,0,2,0] = alpha[t]*(survival1)*omega1
        tp[t,0,2,1] = (1-alpha[t])*(survival1)*omega1
        tp[t,0,2,2] = (survival1)*(1-omega1)
        tp[t,0,2,10] = (1-survival1)
        #state 2: ACLF3
        survival2 = (pre_transplant['state2'][t]+pre_transplant['state2'][t])/2
        tp[t,0,3,0] = (1-pre_transplant['improve'][t])*alpha[t]*(survival2)*omega2
        tp[t,0,3,1] = (1-pre_transplant['improve'][t])*(1-alpha[t])*(survival2)*omega2
        tp[t,0,3,2] = (1-pre_transplant['improve'][t])*(1-omega2)*(survival2)
        tp[t,0,3,3] = (pre_transplant['improve'][t])*alpha[t]*(survival2)*omega2
        tp[t,0,3,4] = (pre_transplant['improve'][t])*(1-alpha[t])*(survival2)*omega2
        tp[t,0,3,5] = (pre_transplant['improve'][t])*(survival2)*(1-omega2)
        tp[t,0,3,10] = (1-survival2)
        tp[t,0,4,0] = (1-pre_transplant['improve'][t])*alpha[t]*(survival2)*omega2
        tp[t,0,4,1] = (1-pre_transplant['improve'][t])*(1-alpha[t])*(survival2)*omega2
        tp[t,0,4,2] = (1-pre_transplant['improve'][t])*(1-omega2)*(survival2)
        tp[t,0,4,3] = (pre_transplant['improve'][t])*alpha[t]*(survival2)*omega2
        tp[t,0,4,4] = (pre_transplant['improve'][t])*(1-alpha[t])*(survival2)*omega2
        tp[t,0,4,5] = (pre_transplant['improve'][t])*(survival2)*(1-omega2)
        tp[t,0,4,10] = (1-survival2)
        tp[t,0,5,0] = (1-pre_transplant['improve'][t])*alpha[t]*(survival2)*omega2
        tp[t,0,5,1] = (1-pre_transplant['improve'][t])*(1-alpha[t])*(survival2)*omega2
        tp[t,0,5,2] = (1-pre_transplant['improve'][t])*(1-omega2)*(survival2)
        tp[t,0,5,3] = (pre_transplant['improve'][t])*alpha[t]*(survival2)*omega2
        tp[t,0,5,4] = (pre_transplant['improve'][t])*(1-alpha[t])*(survival2)*omega2
        tp[t,0,5,5] = (pre_transplant['improve'][t])*(survival2)*(1-omega2)
        tp[t,0,5,10] = (1-survival2)
        #state 3
        survival3 = (pre_transplant['state3'][t]+pre_transplant['state3'][t])/2
        tp[t,0,6,3] = (1-pre_transplant['improve'][t])*alpha[t]*(survival3)*omega3
        tp[t,0,6,4] = (1-pre_transplant['improve'][t])*(1-alpha[t])*(survival3)*omega3
        tp[t,0,6,5] = (1-pre_transplant['improve'][t])*(1-omega3)*(survival3)
        tp[t,0,6,6] = (pre_transplant['improve'][t])*alpha[t]*(survival3)*omega3
        tp[t,0,6,7] = (pre_transplant['improve'][t])*(1-alpha[t])*(survival3)*omega3
        tp[t,0,6,8] = (pre_transplant['improve'][t])*(survival3)*(1-omega3)
        tp[t,0,6,10] = (1-survival3)
        tp[t,0,7,3] = (1-pre_transplant['improve'][t])*alpha[t]*(survival3)*omega3
        tp[t,0,7,4] = (1-pre_transplant['improve'][t])*(1-alpha[t])*(survival3)*omega3
        tp[t,0,7,5] = (1-pre_transplant['improve'][t])*(1-omega3)*(survival3)
        tp[t,0,7,6] = (pre_transplant['improve'][t])*alpha[t]*(survival3)*omega3
        tp[t,0,7,7] = (pre_transplant['improve'][t])*(1-alpha[t])*(survival3)*omega3
        tp[t,0,7,8] = (pre_transplant['improve'][t])*(survival3)*(1-omega3)
        tp[t,0,7,10] = (1-survival3)
        tp[t,0,8,3] = (1-pre_transplant['improve'][t])*alpha[t]*(survival3)*omega3
        tp[t,0,8,4] = (1-pre_transplant['improve'][t])*(1-alpha[t])*(survival3)*omega3
        tp[t,0,8,5] = (1-pre_transplant['improve'][t])*(1-omega3)*(survival3)
        tp[t,0,8,6] = (pre_transplant['improve'][t])*alpha[t]*(survival3)*omega3
        tp[t,0,8,7] = (pre_transplant['improve'][t])*(1-alpha[t])*(survival3)*omega3
        tp[t,0,8,8] = (pre_transplant['improve'][t])*(survival3)*(1-omega3)
        tp[t,0,8,10] = (1-survival3)
        tp[t,0,9,9] = 1
        tp[t,0,10,10] = 1
        #print('Time:::::::::::',t)
        #print(tp[t,0]@tp[t,0])
        #action 1
        tp[t,1,0,9] = 1
        tp[t,1,1,9] = 1
        tp[t,1,2,9] = 1
        tp[t,1,3,9] = 1
        tp[t,1,4,9] = 1
        tp[t,1,5,9] = 1
        tp[t,1,6,9] = 1
        tp[t,1,7,9] = 1
        tp[t,1,8,9] = 1
        tp[t,1,9,9] = 1
        tp[t,1,10,9] = 1
        tp[t,0]/_np.sum(tp[t,0],axis=1)
        tp[t+1] = tp[t]
        new_tp[tt,0] = tp[t,0]@tp[t,0]
        new_tp[tt,1] = tp[t,1]
        print('**********************')
        print(new_tp[tt,0,:,2])
        print('**********************')
        #reward
        new_rw[tt,0:3,0] = 0.4/(365)+_lambda*((1-gamma1[t])*0.4/(365))
        new_rw[tt,3:6,0] = 0.4/(365)+_lambda*((1-gamma2[t])*0.4/(365))
        new_rw[tt,6:9,0] = 0.4/(365)+_lambda*((1-gamma3[t])*0.4/(365))
        #rw[t,0:2,0] = 0.4/(365)
        new_rw[:,0,1] = post_transplant['state1_1'].dot(utl_post)+post_transplant['state1_1'][11]*(15.41*0.6+16.59*0.4)
        new_rw[:,1,1] = post_transplant['state1_1'].dot(utl_post)*rr+post_transplant['state1_1'][11]*(15.41*0.6+16.59*0.4)*rr
        new_rw[:,3,1] = post_transplant['state2_1'].dot(utl_post)+post_transplant['state1_1'][11]*(15.41*0.6+16.59*0.4)
        new_rw[:,4,1] = post_transplant['state2_1'].dot(utl_post)*rr+post_transplant['state1_1'][11]*(15.41*0.6+16.59*0.4)*rr
        new_rw[:,6,1] = post_transplant['state3_1'].dot(utl_post)+post_transplant['state1_1'][11]*(15.41*0.6+16.59*0.4)
        new_rw[:,7,1] = post_transplant['state3_1'].dot(utl_post)*rr+post_transplant['state1_1'][11]*(15.41*0.6+16.59*0.4)*rr
        rw[t,0:9,0] = (0.4/(365))
        rw[t,:,1] = new_rw[tt,:,1]
        rw[t+1,:,:]=rw[t,:,:]
    new_rw= new_rw*50000
    rw = rw*50000
    rw[:,0:9,:] = rw[:,0:9,:]-c_k
    print('old')
    if small==1:
        print("Small")
        a=liver_MDP.MDP_class(tp,term_rw,rw,_lambda,num_absorbing)
        pol1,V = a._back_induction()
        if alpha[0]==0:
            response = 0
        else:
            response = list(pol1[1]).index(1)
        print(V)
    if small==0:
        print('Large')
        a=liver_MDP.MDP_class(new_tp,term_rw,new_rw,new_lambda,num_absorbing)
        pol1,V = a._back_induction()
        if alpha[0]==0:
            response = 0
        else:
            response = list(pol1[1]).index(1)
        #print(V)
    return V,response,tp,new_tp,_lambda,new_lambda,rw,new_rw,pol1



import matplotlib.pyplot as plt
# =============================================================================
# font = {'family' : 'normal',
#         'weight' : 'normal',
#         'size'   : 15}
# 
# plt.rc('font', **font)
# =============================================================================

import pylab as pl
from matplotlib import collections  as mc
#0.6,0.7 to 0.598,0.7
#0.6,0.8 to 0.598,0.8
#0.6,0.9 to 0.598,0.9
#0.7,0.7 to 0.698,0.7
#0.7,0.8 to 0.698,0.8
#0.7,0.9 to 0.698,0.9
alpha = _np.repeat(0.7,28)
alpha2 = _np.repeat(0.7,28)
beta = alpha
beta2 = alpha2
rr=0.5
v1o,r1o,tp,new_tp,_lambda,new_lambda,rw,new_rw,pol=mdp([0.9,0.95,1],alpha,beta,rr,0,0)
v2o,r2o,tp,new_tp,_lambda,new_lambda,rw,new_rw,pol_l=mdp([0.9,0.95,1],alpha2,beta2,rr,1,0)
diff = _np.zeros((3,14))
diff_pol = [[(0,0)for j in range(14)] for i in range(3)]
for t in range(14):
    diff[0,t] = v2o[1,2*t]-v1o[1,t]
    diff[1,t] = v2o[4,2*t]-v1o[4,t]
    diff[2,t] = v2o[7,2*t]-v1o[7,t]
    diff_pol[0][t] = (pol_l[1,2*t],pol[1,t])
    diff_pol[1][t] = (pol_l[4,2*t],pol[4,t])
    diff_pol[2][t] = (pol_l[7,2*t],pol[7,t])
import matplotlib.pyplot as plt
#diff = diff[:,0:7]
plt.plot(diff[0],color='black')
plt.plot(diff[1],color='blue')
plt.plot(diff[2],color='red')
plt.show()

plt.plot(diff[:,0],label = '2')
plt.plot(diff[:,1],label = '4')
plt.plot(diff[:,2],label = '6')
plt.legend()
plt.show()
print('threshold-small')
print(list(pol[1]).index(1)*2)
print(list(pol[4]).index(1)*2)
print(list(pol[7]).index(1)*2)
print('threshold-large')
print(list(pol_l[1]).index(1))
print(list(pol_l[4]).index(1))
print(list(pol_l[7]).index(1))
print('max dt')
print(_np.round(max(diff[0]),0))
print(_np.round(max(diff[1]),0))
print(_np.round(max(diff[2]),0))
print('time max dt')
print(list(diff[0]).index(max(diff[0]))*2+2)
print(list(diff[1]).index(max(diff[1]))*2+2)
print(list(diff[2]).index(max(diff[2]))*2+2)
# case 1 o =0.5 rr=0.7
# case 2 o =0.69 rr=0.9
#pop aclf =3: 1364
#pop aclf>3 889
#pop aclf2 1721+622=2343
print(max(diff[0])*2343+max(diff[1])*1364+max(diff[2])*889)
print(27209*2343+24256*1364+27044*889)

print('########################')
print('check thm 4')

# =============================================================================
# max_ = -99999
# # =============================================================================
# # for t in range(27):
# #     rw[t,2,1] = rw[t,2,0]+_lambda*tp[t,0,2,:]@rw[t+1,:,1]
# #     rw[t,5,1] = rw[t,5,0]+_lambda*tp[t,0,5,:]@rw[t+1,:,1]
# #     rw[t,8,1] = rw[t,8,0]+_lambda*tp[t,0,8,:]@rw[t+1,:,1]
# # =============================================================================
# for t in range(27):
#     for ss in range(9):
#         if rw[t,ss,1]<rw[t,ss,0]+_lambda*tp[t,0,ss,:]@rw[t+1,:,1]:
#             rw[t,ss,1]=rw[t,ss,0]+_lambda*tp[t,0,ss,:]@rw[t+1,:,1]
# s = 7
# for t in range(28):        
#     print(t)
#     print(rw[t,s,0]+_lambda*tp[t,0,s,:]@rw[t,:,1]-v1o[s,int(0.5*t)])
#     print('************')
#     #print(-new_lambda*c_k +new_rw[t,s,0]+new_lambda*new_tp[t,0,s,:]@new_rw[t,:,1]-v1o[s,int(0.5*t)])
#     if rw[t,s,0]+_lambda*tp[t,0,s,:]@rw[t,:,1]-v1o[s,int(0.5*t)]> max_:
#         max_ = rw[t,s,0]+_lambda*tp[t,0,s,:]@rw[t,:,1]-v1o[s,int(0.5*t)]
# print(max_)
# =============================================================================

'''
def upper(state):
    x = new_lambda*new_tp[0,0,state]
    y = _np.array([new_rw[0,0,1]-new_rw[0,0,0],0,0,0])+new_lambda*(new_tp[0,0]@v2o[:,0]-new_tp[-1,0]@v1o[:,6])
    x1 = new_rw[0,1,0]-new_rw[0,1,1]+new_lambda*new_tp[0,0,1]@v2o[:,0]
    return max(x1,x@y)
print('$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$')
print(upper(1))
print('$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$')
'''
# =============================================================================
# gamma1 = 0.1
# gamma2 = 0.2
# gamma3 = 0.3
# tau = 0.1
# o = 0.5
# p = _np.zeros((8,8))
# p[0,0] = o*(1-gamma1)
# p[0,1] = (1-o)*(1-gamma1)
# p[0,7] =  gamma1
# p[1,:] = p[0,:]
# p[2,0] = tau*o*(1-gamma2)
# p[2,1] = tau*(1-o)*(1-gamma2)
# p[2,2] = (1-tau)*o*(1-gamma2)
# p[2,3] = (1-tau)*(1-o)*(1-gamma2)
# p[2,7] = gamma2
# p[3,:] = p[2,:]
# p[4,2] = tau*o*(1-gamma3)
# p[4,3] = tau*(1-o)*(1-gamma3)
# p[4,4] = (1-tau)*o*(1-gamma3)
# p[4,5] = (1-tau)*(1-o)*(1-gamma3)
# p[4,7] = gamma3
# p[5,:] = p[4,:]
# p[6,6] = 1
# p[7,7] = 1
# 
# def compute_P(gamma1,gamma2,gamma3,o,tau):
#     P = _np.zeros((8,8))
#     P[0,0] = o*(1-gamma1)**2
#     P[0,1] = (1-o)*(1-gamma1)**2
#     P[0,7] = 1-(1-gamma1)**2
#     P[1,:] = P[0,:]
#     P[2,0] = tau*o*(1-gamma2)*(1-gamma1)+tau*(1-tau)*o*(1-gamma2)**2
#     P[2,1] = tau*(1-o)*(1-gamma2)*(1-gamma1)+tau*(1-tau)*(1-o)*(1-gamma2)**2
#     P[2,2] = (1-tau)**2*o*(1-gamma2)**2
#     P[2,3] = (1-tau)**2*(1-o)*(1-gamma2)**2
#     P[2,7] = 1-(1-tau)*(1-gamma2)**2-tau*(1-gamma1)*(1-gamma2)
#     P[3,:] = P[2,:]
#     P[4,0] = tau**2*o*(1-gamma2)*(1-gamma3)
#     P[4,1] = tau**2*(1-o)*(1-gamma2)*(1-gamma3)
#     P[4,2] = tau*(1-tau)*o*(1-gamma2)*(1-gamma3)+tau*(1-tau)*o*(1-gamma3)**2
#     P[4,3] = tau*(1-tau)*(1-o)*(1-gamma2)*(1-gamma3)+tau*(1-tau)*(1-o)*(1-gamma3)**2
#     P[4,4] = (1-tau)**2*o*(1-gamma3)**2
#     P[4,5] = (1-tau)**2*(1-o)*(1-gamma3)**2
#     P[4,7] = 1-tau*(1-gamma2)*(1-gamma3) - (1-tau)*(1-gamma3)**2
#     P[5,:] = P[4,:]
#     P[6,6] = 1
#     P[7,7] = 1
#     return P
# P = compute_P(gamma1,gamma2,gamma3,o,tau)
# print(p@p)
# print(P)
# w, v = LA.eig(P)
# for i in range(len(w)):
#     if w[i]<0:
#         w[i]=0
# w = w**(0.5)
# new = v@_np.diag(w)@_np.linalg.inv(v)
# for i in range(len(new)):
#     for j in range(len(new[i])):
#         new[i,j] = max(0,_np.round(new[i,j],6))
# 
# print(new)
# print(p)
# =============================================================================
def sum(i,j):
    sol = 0
    sol2 =0
    j=26-j+1
    for x in range(1,i+1):
        sol+=x
    for x in range(j,i+1):
        sol2+=x
    return sol2/sol