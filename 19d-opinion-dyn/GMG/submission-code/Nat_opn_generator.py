#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Jul  6 10:41:47 2023

@author: glory
"""

import numpy as np
import math
import random

# Generating gaussian and bigaussian random natural opinion for 2 parties
# INPUTS:
#   mu: mean of the whole distribution
#   Delta: polarization coefficient, i.e., distance between the# INPUT:
# num_ppl:number of agents
# p:number of parties
#OUTPT:
#X:natural opinion (num_ppl,p) two peaks of the distribution
#   sd: standard deviation of the gaussian(s)
#   num_ppl: number of opinions to be generated
#   R: percentage of agents in each peak
# OUTPUTS:
#   x0: array of natural opinions (num_ppl,2)

def Nat_opn_2(mu,Delta,sd,num_ppl,R=50):
    r1 = int(np.round(R*0.01*num_ppl))
    r2 = num_ppl-r1
    y1 = np.random.normal(mu-Delta/2,sd,r1)
    while np.sum(y1<-1) or np.sum(y1>1):
        y1[y1<-1]=np.random.normal(mu-Delta/2,sd,np.sum(y1<-1))
        y1[y1>1]=np.random.normal(mu-Delta/2,sd,np.sum(y1>1))
    y2 = np.random.normal(mu+Delta/2,sd,r2)
    while np.sum(y2<-1) or np.sum(y2>1):
        y2[y2<-1]=np.random.normal(mu+Delta/2,sd,np.sum(y2<-1))
        y2[y2>1]=np.random.normal(mu+Delta/2,sd,np.sum(y2>1))
    x = np.concatenate((y1,y2))
    x0 = np.zeros((num_ppl,2))
    x0[:,0]=0.5-0.5*x
    x0[:,1]=1-x0[:,0]
    return x0

# Draws a random point uniformly from an arbitrary simplex, defined by the columns of S.
# Follows the idea presented in doi.org/10.13140/RG.2.1.3807.6968.
def unit_simplex_arb(S):
    n,p=S.shape
    
    zs=np.random.uniform(0,1,p+1)
    zs[0],zs[p]=1,0
    
    lambda_s=[1.]
    for j in range(1,p):
        lambda_s.append(zs[j]**(1/(p-j)))
    lambda_s.append(0.)
    
    p_lambda_s=np.cumprod(lambda_s)
    
    V = (1-lambda_s[1])*p_lambda_s[0]*S[:,0]
    
    for i in range(1,p):
        V=np.vstack((V, (1-lambda_s[i+1])*p_lambda_s[i]*S[:,i]))
    
    x=np.sum(V,axis=0)
    
    return x

def rand_opinion_pos(p,V,k):
    perm = gen_rand_perm_vect(np.concatenate((np.zeros(k-1),np.ones(p-k))))  # Take one of the permutations of k-1 zeros and p-k ones.
	
    idx=np.zeros(p, dtype=int)
    left=(np.setdiff1d((1-perm)*np.arange(2,p+1,1),np.array((0)))).astype('int') # Indices of the parties at the left of party k in the order of the opinion.
    right=(np.setdiff1d((perm)*np.arange(2,p+1,1),np.array((0)))).astype('int') # Indices of the parties at the right of party k in the order of the opinions.
    
    idx[0]=k
    idx[left-1]=np.arange(k-1,0,-1)
    idx[right-1]=np.arange(k+1,p+1)
    
    # Get the summits of the orthoscheme as summits of the simplex where the opinion is to be drawn.
    S = np.empty((p, p))
    string = ""
    for i in range(p):
        string += str(int(idx[i]))
        S[:, i] = V[''.join(sorted(string))]
    
    # Generate the opinion.
    x = unit_simplex_arb(S)
    
    return x
  
# Draws one opinion uniformly in the admissible space.
# The list of summits is already given in 'V'. One can get 'V' by using the function 'admissible_summits'.
def rand_opinion(p,V):
    k=random.randint(1,p)  # Selects one of the parties
    
    x=rand_opinion_pos(p,V,k)
    return x


# Generates a random permutation of x.
def gen_rand_perm_vect(x):
    n=len(x)
    p=x.copy()
    
    for i in range(0,n-1):
        j=i+random.randint(0,n-i-1)
        if j!=i:
            p[i],p[j]=p[j],p[i]  
    return p    


# Slides the point x (typically in the simplex) towards the centroid of the unitary p-simplex. 
# For α = 1, the summit does not move, and for α = 0, the summit reaches the centroid.
def slide_summit_vec(x,a):
    if a<0:
        print("Warning a < 0")
    elif a >1:
        print("Warning a >1")
        
    p=len(x)
    c=1/p*np.ones(p)   
    
    return (1-a)*c +a*x


# Generates the list of summits of interest in the unitary simplex with 'p' summits.
def vertices(p):
    vertex={}
    if p==1:
        vertex["1"]=[1.]
    else:
        v0=vertices(p-1)
        ks=v0.keys()
        
        s=np.identity(p)[p-1]
        vertex[f"{p}"]=s
        
        for k in ks:
            vertex[k]=np.append(v0[k],0)
            l=len(k)
            vertex[f"{k}{p}"]=(l*vertex[k]+s)/(l+1)
            
    return vertex
    
# Defines the summits of the admissible opinion space.
#INPUT:
# p: number of parties
#OUTPUT:
# V:summits of the admissible opinion space
def admissible_summits(p):
    V=vertices(p)
    
    c=""  # Name of the centroid
    for k in range(1,p+1):
        c=c+str(k)
        
    # Slide each of the summits to the appropriate position to equalize the volume of each party.
	
    for k in range(1,p+1):
        # Party 'k' has exactly 'binomial(p-1,k-1)' admissible orthoschemes. Therefore, extremists have only one. The ratio between the admissible volumes is the given by 'r', which is the reduction factor for each volume
        r=1/math.comb(p-1,k-1)
        V[f"{k}"]=slide_summit_vec(V[f"{k}"],r)
    return V    
      
#Random generation with agents of 'k' th party  
# INPUT:
# num_ppl:number of agents
# p:number of parties
# k:party to which the agents will vote for
#OUTPT:
#X:natural opinion (num_ppl,p)
def Random_generation_pos(num_ppl, p, k):
    
    AS = admissible_summits(p)
    
    X=np.zeros((num_ppl,p))
    for num in range(num_ppl):
        X[num]=rand_opinion_pos(p, AS, k)
        
    return X    

# Random generation 
# INPUT:
# num_ppl:number of agents
# p:number of parties
#OUTPT:
#X:natural opinion (num_ppl,p)
def Random_generation(num_ppl,p):
    
    AS = admissible_summits(p)
    
    X=np.zeros((num_ppl,p))
    for num in range(num_ppl):
        X[num]=rand_opinion(p, AS)
        
    return X  
    