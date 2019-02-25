#!/usr/bin/python
import numpy as np

def ml_zipf(xx):
	x = np.array(xx)
	
	mi = min(x)
	ma = max(x)
	
	l = .01
	u = 10.
	
	s = 10.
	
	f = []
	for i in range(mi,ma+1):
		f.append(sum(x==i)/len(x))
	
	suf = sum(f)
	su = sum(np.multiply(f,np.log(range(mi,ma+1))))
	
	error = 1000
	
	while error > 1e-5:
		ss = np.linspace(l,u,100)
		Hs = []
		for s in ss:
			Hs.append(sum(1/(np.power(range(mi,ma+1),s))))
		
		La = -np.multiply(ss,su) - suf*np.log(Hs)
			
		idx = np.argmin(La)
		l = max(l,ss[max(idx-1,0)])
		u = min(u,ss[min(idx+1,99)])
		error = u - l
	
	s = (l+u)/2
	
	return s
	

def ml_clauset(xx):
	x = np.array(xx)
	
	mi = min(x)
	ma = max(x)
	n = len(x)
	
	l = .01
	u = 10.
	s = 10.
	
	sum_data = sum(np.log(x))
	
	error = 1000
	
	while error > 1e-5:
		ss = np.linspace(l,u,100)
		zeta = []
		for si in ss:
			zeta.append(sum(np.power(range(mi,ma+1),-si)))
		
		L = -n*np.log(zeta) - np.multiply(ss,sum_data)
		
		idx = np.argmax(L)
		l = max(l,ss[max(idx-1,0)])
		u = min(u,ss[min(idx+1,99)])
		error = u - l
	
	s = (l+u)/2
	
	return s


