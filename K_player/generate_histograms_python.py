import matplotlib.pyplot as plt
import numpy as np
import csv
import math
from scipy.special import zeta

from ml_zipf import *
from journals import *

nbins = 300
logbins = True
linbins = False
zipf_plot = False
tail_plot = True

js = ["prl",]

a1 = 0; a2 = 0; b1 = 0; b2 = 0

for j in js:
# fetch data
	print("Collecting data")
	dat = list(csv.reader(open("./data/"+j+".txt","r"), delimiter="\t"))
	num = []
	for i in range(1,len(dat)-2):
		num.append(int(dat[i][1]))
	
	num = np.array(num)
	mi = min(num)
	ma = max(num)
	vals = list(set(num))
	
	if logbins:
		print("Entering logbins")
		plt.figure(j)
	
# plot the data	
		for i in vals:
			plt.plot([i,i],[1e-8,sum(num==i)/len(num)],"-b",linewidth=2,color=journals_colors[j][0])
	
# determine the power law exponent and the normalizing coefficient
		print("Computing exponent")
		s = ml_clauset(num)
		C = 1/zeta(s,mi)

# plot the power law fit and the max number of articles
		Hs = sum(1/np.power(range(mi,ma+1),-s))
		z = C*np.power(range(mi,ma+1,10),-s)
		
		plt.plot(range(mi,ma+1,10),z,"--k",label="Zipf's law fit, s = "+str(round(s,3)),linewidth=2)
		max_k = math.ceil(np.power(len(num)*C,1/s))
		plt.plot([max_k,max_k],[.5/len(num),2*max(z)],":k",linewidth=2)
		
		plt.title(journals_code[j]+", max. # articles predicted: "+str(math.ceil(max_k)))
		plt.ylabel("Number of authors")
		plt.xlabel("Number of articles")
		plt.axis([.9,max([2*ma,2*max_k]),.5/len(num),2*max(z)])
		plt.loglog()
		plt.legend()
	
	if linbins:
		# TODO ???
		1 + 1
	
	if zipf_plot:
		vals = []
		for i in range(mi,ma+1):
			vals.append(sum(num==i))
		
		plt.figure(j+" 3")
		plt.plot([mi,ma],[0,0],"--k")
		plt.plot(range(mi,ma+1),np.abs(vals/sum(vals) - z),"oy",label="Zipf")
		plt.semlogy()
		plt.title(journals_codes[j]+", goodness-of-fit")
		plt.xlabel("Number of authors")
		plt.ylabel("Error of Zipf's fit")
		plt.legend()
	
	if tail_plot:
		print("Entering tail_plot")
#		s = ml_clauset(num)
#		C = 1/zeta(s,mi)
		print("Computing cdf's")
		cdf = [C*mi**(-s),]
		ecdf = [sum(num==mi)/len(num),]
		for i in range(mi+1,ma+1):
			cdf.append(cdf[len(cdf)-1] + C*i**(-s))
			ecdf.append(sum(num<=i)/len(num))
		
		ecdf = np.array(ecdf)
		cdf = np.array(cdf)
		print("Plotting cdf's")
		plt.figure(j+" 666")
		plt.plot(range(mi,ma+1),1+1e-8-ecdf,color=journals_colors[j][0],linewidth=2)
		plt.plot(range(mi,ma+1),1+1e-8-cdf,"--k",linewidth=2)
		plt.loglog()
		plt.axis([mi,2*ma,min([(1-max(ecdf[0:len(ecdf)-2])),(1-max(cdf[0:len(cdf)-2]))])/2,1])
		plt.title(journals_code[j]+": Tails comparison")
	

plt.show()

print("All done.")	
