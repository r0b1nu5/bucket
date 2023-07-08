#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Jul  6 11:00:05 2023

@author: glory
"""

import Nat_opn_generator as OPN_GEN

import numpy as np

def type_natural_opinion(Num_ppl,Num_party,NAT_TYPE):
    if NAT_TYPE=='Bigaussian':
        x0=OPN_GEN.Nat_opn_2(mean,Delta,sd,Num_ppl,R=50)
    else:
        x0=OPN_GEN.Random_generation(Num_ppl, Num_party)
    return x0

#Function to compute the outcome of the election
#INPUT:
# x0: natural opinion (Num_ppl,Num_party)
# epsilon: confidence bound parameter (float)
# W:external influence (Num_ppl,Num_party)

def outcome(x0, epsilon, W):
    Num_ppl = len(x0)
    
    A=np.zeros((Num_ppl,Num_ppl))
    for i in range(Num_ppl):
        A[i,:]=(np.linalg.norm(np.tile(x0[i],(Num_ppl,1))-x0,ord=construction_norm,axis=1)<=epsilon)*1
    A=A-np.identity(Num_ppl) # Adjacency Matrix
   
    Mat_inv=np.linalg.inv(np.matmul(np.diag(np.array([1/i if i!=0 else 0 for i in np.sum(A,axis=1)])),np.diag(np.sum(A,axis=1))-A)+np.identity(len(x0)))
    
    y=np.matmul(Mat_inv,x0+W) 
    
    return y
 
if __name__ == "__main__": 
    
    construction_norm=1
   
    Num_ppl= int(input("Enter number of agents: "))
    epsilon = float(input("Enter confidence bound value between 0 to 2 : "))
    
    while epsilon>2 or epsilon<0:
        epsilon = float(input("Invalid entry! Enter confidence bound value between 0 to 2 : "))
    
    Num_party=int(input("Enter number of parties: "))
    if Num_party==2:
        NAT_TYPE=input("Enter the natural opinion type: 'Bigaussian' or 'Uniform' :")
        if NAT_TYPE=='Bigaussian':
            mean=float(input('Enter mean value of the whole distribution within the interval I= [0,0.1] : '))
            while mean<0 or mean>0.1:
                mean=float(input('Invalid entry! Enter mean value of the whole distribution within the interval I= [0,0.1] : '))
    
            Delta=float(input('Enter Delta value within the interval I= [0,1] : '))
            while Delta<0 or Delta>1:
                Delta=float(input('Invalid entry! Enter Delta value within the interval I= [0,1] : '))
            sd=float(input('Enter standard deviation of the distribution : '))
    else:
        NAT_TYPE='Uniform'
    
    W=np.zeros((Num_ppl,Num_party)) # External influence, change the components of W to add external influence
    
    x0=type_natural_opinion(Num_ppl,Num_party,NAT_TYPE) # Generate natural opinion
    y = outcome(x0,epsilon,W) # Compute final opinion
    
    Num_votes=np.sum(np.argsort(x0).argsort()==Num_party-1,axis=0)
    print('Number of votes acc. to natural opinion:', Num_votes)
    
    Num_votes_final=np.sum(np.argsort(y).argsort()==Num_party-1,axis=0)
    print('Number of votes acc. to final opinion:', Num_votes_final)
    
        
    