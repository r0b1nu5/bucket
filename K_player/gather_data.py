#!/usr/bin/python3
import csv
import numpy as np

from journals import *

np = ["prd","chaos-18"]
#for j in journals_short:
for j in np:
	dat = list(csv.reader(open("./data/"+j+".txt","r"), delimiter="\t"))
	num = []
	for i in range(1,len(dat)-2):
		num.append(int(dat[i][1]))
	
	num = np.array(num)
	
	vals = np.array(list(set(num)))
	
	distr = []
	for i in vals:
		distr.append(sum(num==i)/len(num))
	
	distr = np.array(distr)
	
	with open("./data/"+j+"_vals.csv",mode="w") as f:
		wr = csv.writer(f,delimiter=",")
		wr.writerow(vals)
	
	with open("./data/"+j+"_distr.csv",mode="w") as f:
		wr = csv.writer(f,delimiter=",")
		wr.writerow(distr)
	

