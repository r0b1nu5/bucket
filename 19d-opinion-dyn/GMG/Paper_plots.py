#!/usr/bin/env python
# coding: utf-8

# In[9]:


import numpy as np
import matplotlib.pyplot as plt

#import geopandas
import pandas as pd
import json
import plotly.express as px
import plotly.io as pio
pio.renderers.default = 'browser'

import plotly.graph_objects as go
from plotly.subplots import make_subplots
import requests
from scipy import stats
import kaleido


# In[ ]:





# In[ ]:





# In[3]:


path='PAPER_DATA/VALIDATION/US_CONST/'

Influence_array1=np.array((0.25,0.5,1,2,3,4,5,6,7,8,9,10,12,15,20,25,30,35,40,45,50,75,100))
neps=25
#neps=151
years=5
NOS=100
Num_type=4
year_start=2012
year_end=2020
Num_districts=435

YEARS=np.linspace(year_start,year_end,years,dtype=int)
Districts_influenced=78
#NUM_DISTRICTS_CHANGED=np.zeros((len(Influence_Array),Num_type, years, NOS, neps,Num_districts))
SUCCESS=np.zeros((years,Num_type,NOS,neps,len(Influence_array1),Num_districts))
for num_type in range(Num_type):
    i=0
    for year in YEARS:
        data=np.load(path+'usrepub_'+str(num_type+1)+'_'+str(year)+'.npz')
        e_range=data['e_range']
        SUCCESS[i,num_type]=data['SUCCESS']
        Influence_percentage_array=data['Influence_percentage_array']
        i=i+1 
Success=np.mean(np.sum(SUCCESS,axis=2),axis=0)  

df=pd.read_csv('PAPER_DATA/HOR_DATA/Election_results.csv', converters={'CD_NUM (2010)': str, 'STATE_FP (2010)': str, 'ID (2010)' : str, 'CD_NUM': str , 'STATE_FP': str , 'ID' : str})
years=np.linspace(2012,2020,5, dtype=int)
change=np.zeros((Num_districts))
for i in range(len(years)-1):
    change+=abs((df['Republican '+str(years[i])]>50)*1-(df['Republican '+str(years[i+1])]>50)*1)


Correlation_Matrix=np.zeros((Num_type, neps, len(Influence_array1)))
for i in range(Num_type):
    for j in range(neps):
        for k in range(len(Influence_array1)):
            Correlation_Matrix[i,j,k]=np.corrcoef(Success[i,j,k,:], change)[0,1]
    


# In[4]:
###################################################################
############################ figsup ###############################
###################################################################

c=['r','b','g','m']
LABEL=['D1','D2','D3','D4']
#plt.figure(figsize=(10,4))
Array1=np.zeros((Num_type,len(Influence_array1)))
Array2=np.zeros((Num_type,len(Influence_array1)))
Array3=np.zeros((Num_type,len(Influence_array1)))

for inf in range(len(Influence_array1)):
    for num_type in range(Num_type):
        Array1[num_type,inf]=np.mean(Correlation_Matrix,axis=1)[num_type,inf].round(decimals=3)
        Array2[num_type,inf]=np.std(Correlation_Matrix,axis=1)[num_type,inf]
        Array3[num_type,inf]=np.max(Correlation_Matrix,axis=1)[num_type,inf].round(decimals=3)
          

plt.figure(figsize=(10,5))
Num_type=3
for num_type in range(Num_type):
    plt.plot(Influence_array1[0:16],Array1[num_type][0:16],'o',color=c[num_type],label=LABEL[num_type])
    plt.errorbar(Influence_array1[0:16],Array1[num_type][0:16],Array2[num_type][0:16],color=c[num_type],ls='none')
plt.xlabel('Influence percentage')

plt.ylabel('r\u0305')   
plt.legend()

Max_inf=Influence_array1[np.argmax(Array1,axis=1)]
for num_type in range(Num_type):
    plt.axvline(Max_inf[num_type],linestyle='--',color=c[num_type])
    
#plt.savefig('/home/glory/Dropbox/Opinion Dynamics/paper/Figures/figs1.pdf',bbox_inches='tight')


# In[5]:


data=df['Districts']
Data = pd.DataFrame(data, columns=['Districts'])
Data.insert(len(Data.columns), 'ID', df['ID'])
Data.insert(len(Data.columns), 'Actual change', change)
Num_type=4
for type in range(Num_type):
    for inf in range(len(Influence_array1)):
        eps_pos=np.argmax(Correlation_Matrix,axis=1)[type,inf]
        epsilon=str(e_range[np.argmax(Correlation_Matrix,axis=1)[type,inf]].round(decimals=3))
        corr=str(np.max(Correlation_Matrix,axis=1)[type,inf].round(decimals=3))
        S=Success[type,eps_pos,inf]
        name ='Inf: '+str(Influence_array1[inf])+' Type: '+str(type+1)+' eps:'+epsilon+' Corr: '+ corr
        Data.insert(len(Data.columns),name,S)        


# In[6]:
###################################################################
####################### fig1ab ######################################
###################################################################

geojson = json.load(open("PAPER_DATA/Congressional-Districts/src/eclair/US_Congressional_Districts.geodata", 'r'))

Max_inf= np.array((5,5,5))
ABC=['(b)','(b)','(d)']
subtitles=['(a) Historic Data']
Num_type=3
for num_type in range(Num_type):
    inf=np.where(Influence_array1==Max_inf[num_type])[0][0]
    subtitles.append(ABC[num_type]+' \u03B5 = '+str(e_range[np.argmax(Correlation_Matrix,axis=1)[num_type,inf]].round(decimals=3))+', r = '+ str(np.max(Correlation_Matrix,axis=1)[num_type,inf].round(decimals=3))+', IB = '+str(Influence_array1[inf]*0.1))

white_space='                                                                                                   '
ABCD=['(a)'+white_space, '(b)'+white_space, '(c)'+white_space, '(d)'+white_space]

for feature in geojson['features']:
    feature['properties']['ID'] = feature['properties']['STATEFP'].__str__()+feature['properties']['CD115FP'].__str__()

fig = make_subplots(
rows=1, cols=2,
specs=[
    [{"type": "choropleth"},{"type": "choropleth"}],
],
#subplot_titles = [subtitles[0],subtitles[1],subtitles[2],subtitles[3]],
#subplot_titles = [subtitles[0],subtitles[2]],    
vertical_spacing = 0.15,  
horizontal_spacing = 0.00000015,
)

fig.add_trace(trace=go.Choropleth(
    geojson=geojson,
    featureidkey='properties.ID',
    locations=Data['ID'],
    z=Data['Actual change']*25,
    zmin=0,
    zmax=100,
    colorscale='ylorrd',
), row=1, col=1)

ini=np.array((0.48,0.73,1))

num_type=0
inf=np.where(Influence_array1==Max_inf[num_type])[0][0]
eps_pos=np.argmax(Correlation_Matrix,axis=1)[num_type,inf]
epsilon=str(e_range[np.argmax(Correlation_Matrix,axis=1)[num_type,inf]].round(decimals=3))
corr=str(np.max(Correlation_Matrix,axis=1)[num_type,inf].round(decimals=3))
item ='Inf: '+str(Influence_array1[inf])+' Type: '+str(num_type+1)+' eps:'+epsilon+' Corr: '+ corr
fig.add_trace(trace=go.Choropleth(
    geojson=geojson,
    featureidkey='properties.ID',
    locations=Data['ID'],
    z=Data[item],
    zmin=0,
    zmax=100,
    colorscale='ylorrd', 
), row=1, col=2)

fig.update_geos(scope='usa')

fig.update_layout(
    title='', title_x=0.5, width=1500, height=600)
fig.show()   

#pio.write_image(fig,'fig1.pdf')

# In[7]:
###################################################################
######################## fig1absup ##############################
###################################################################

geojson = json.load(open("PAPER_DATA/Congressional-Districts/src/eclair/US_Congressional_Districts.geodata", 'r'))

#Max_inf=Influence_array1[np.argmax(Array1,axis=1)]

Max_inf= np.array((5,5,5))
ABC=['(b)','(b)','(d)']
subtitles=['(a) Historic Data']
Num_type=3
for num_type in range(Num_type):
    inf=np.where(Influence_array1==Max_inf[num_type])[0][0]
    subtitles.append(ABC[num_type]+' \u03B5 = '+str(e_range[np.argmax(Correlation_Matrix,axis=1)[num_type,inf]].round(decimals=3))+', r = '+ str(np.max(Correlation_Matrix,axis=1)[num_type,inf].round(decimals=3))+', IB = '+str(Influence_array1[inf]*0.1))
#print(subtitles)   

white_space='                                                                                                   '
ABCD=['(a)'+white_space, '(b)'+white_space, '(c)'+white_space, '(d)'+white_space]

for feature in geojson['features']:
    feature['properties']['ID'] = feature['properties']['STATEFP'].__str__()+feature['properties']['CD115FP'].__str__()

fig = make_subplots(
rows=1, cols=2,
specs=[
    [{"type": "choropleth"},{"type": "choropleth"}],
],
#subplot_titles = [subtitles[0],subtitles[1],subtitles[2],subtitles[3]],
#subplot_titles = [subtitles[1],subtitles[2]],    
vertical_spacing = 0.15,  
horizontal_spacing = 0.00015,
)

ini=np.array((0.48,0.73,1))

num_type=1
inf=np.where(Influence_array1==Max_inf[num_type])[0][0]
eps_pos=np.argmax(Correlation_Matrix,axis=1)[num_type,inf]
epsilon=str(e_range[np.argmax(Correlation_Matrix,axis=1)[num_type,inf]].round(decimals=3))
corr=str(np.max(Correlation_Matrix,axis=1)[num_type,inf].round(decimals=3))
item ='Inf: '+str(Influence_array1[inf])+' Type: '+str(num_type+1)+' eps:'+epsilon+' Corr: '+ corr
fig.add_trace(trace=go.Choropleth(
    geojson=geojson,
    featureidkey='properties.ID',
    locations=Data['ID'],
    z=Data[item],
    zmin=0,
    zmax=100,
    colorscale='ylorrd',
    colorbar = dict(x=0.95,y=0.5,len=1.2,title='Per. of change')    
), row=1, col=1)


num_type=2
inf=np.where(Influence_array1==Max_inf[num_type])[0][0]
eps_pos=np.argmax(Correlation_Matrix,axis=1)[num_type,inf]
epsilon=str(e_range[np.argmax(Correlation_Matrix,axis=1)[num_type,inf]].round(decimals=3))
corr=str(np.max(Correlation_Matrix,axis=1)[num_type,inf].round(decimals=3))
item ='Inf: '+str(Influence_array1[inf])+' Type: '+str(num_type+1)+' eps:'+epsilon+' Corr: '+ corr
fig.add_trace(trace=go.Choropleth(
    geojson=geojson,
    featureidkey='properties.ID',
    locations=Data['ID'],
    z=Data[item],
    zmin=0,
    zmax=100,
    colorscale='ylorrd',
    colorbar = dict(x=0.95,y=0.5,len=1.2,title='Per. of change')    
), row=1, col=2)

fig.update_geos(scope='usa')

fig.update_layout(
    title='', title_x=0.5, width=1500, height=600)
fig.show()   

#pio.write_image(fig,'fig1sup.pdf')

# In[8]:
###################################################################
################ fig1csup ###########################################
###################################################################

Districts=np.arange(1,436,1)
Max_inf=Influence_array1[np.argmax(Array1,axis=1)]

ABC=['(b)','(c)','(d)','(e)']
subtitles=['(a) Historic Data']
Num_type=4
S=np.zeros((Num_type,Num_districts))
for type in range(Num_type):
    inf=np.where(Influence_array1==Max_inf[type])[0][0]
    eps_pos=np.argmax(Correlation_Matrix,axis=1)[type,inf]
    subtitles.append(ABC[type]+' \u03B5 = '+str(e_range[np.argmax(Correlation_Matrix,axis=1)[type,inf]].round(decimals=3))+', r = '+ str(np.max(Correlation_Matrix,axis=1)[type,inf].round(decimals=3))+', Per. = '+str(Influence_array1[inf]))
    S[type]=Success[type,eps_pos,inf]

Districts=np.arange(1,436,1)
num_type=0
plt.figure(figsize=(13,5))
plt.plot(Districts,(change/4)*100,'--',color='k',label='Historic Data',linewidth=2)
plt.plot(Districts,S[type],color='r',label='Simulation')
plt.legend(loc='upper right')
plt.xlabel('Districts')
plt.ylabel('Per. of change')


# In[8]:


fig,ax=plt.subplots(1,2,figsize=(18,5))
for num_type in range(1,3): 
    ax[num_type-1].plot(Districts,(change/4)*100,'--',color='k',label='Historic Data',linewidth=2)
    ax[num_type-1].plot(Districts,S[type],color='r',label='Simulation')
    ax[num_type-1].legend(loc='upper right')
    ax[num_type-1].set_xlabel('Districts')
    ax[num_type-1].set_ylabel('Per. of change [%]')

fig.show()


# In[10]:
####################################################
################## fig2 ############################
####################################################




data=np.load('PAPER_DATA/VARIATION_MU_DELTA_P/VAR_FILES/Polarization_delta.npz')

Num_ppl=data['Num_ppl']
e_range=data['e_range']
DELTA_vals=data['DELTA_vals']
MEAN_vals=data['MEAN_vals']
NUM_AGENTS_INFLUENCED=data['NUM_AGENTS_INFLUENCED']

plt.figure(figsize=(6,5))

xlist = e_range
ylist = DELTA_vals
X, Y = np.meshgrid(xlist, ylist)
Z=(np.mean(NUM_AGENTS_INFLUENCED,axis=2)[:,0,:]/Num_ppl)*100
cp = plt.contourf(X, Y, Z,len(DELTA_vals),cmap='YlOrRd')
#ax.set_zlim(.4,1.4)


plt.ylabel('$\Delta$',fontsize='15')
plt.xlabel('$\epsilon$',fontsize='15')
#plt.colorbar(cp,label='$\\xi$') # Add a colorbar to a plot  
cbar1 = plt.colorbar()
#cbar1.ax.set_title('Per. of agents inf.') 
cbar1.ax.set_title('%') 
plt.xlim(right=1.25)
plt.show()

#plt.savefig('/home/glory/Dropbox/Opinion Dynamics/paper/Figures/Variation_delta.pdf',bbox_inches='tight')


# In[11]:
####################################################
################## fig2sup ############################
####################################################


fig,axes=plt.subplots(1,3,figsize=(20,5))

data=np.load('PAPER_DATA/VARIATION_MU_DELTA_P/VAR_FILES/Mean_variation.npz')
Num_ppl=data['Num_ppl']
e_range=data['e_range']
DELTA_vals=data['DELTA_vals']
MEAN_vals=data['MEAN_vals']
NUM_AGENTS_INFLUENCED=data['NUM_AGENTS_INFLUENCED']

xlist=e_range
ylist=MEAN_vals
X1, Y1 = np.meshgrid(xlist, ylist)
Z1=(np.mean(NUM_AGENTS_INFLUENCED,axis=2)[0]/Num_ppl)*100
cp = axes[0].contourf(X1, Y1, Z1,len(MEAN_vals),cmap='YlOrRd')
axes[0].set_ylabel('$\mu$',fontsize='15')
axes[0].set_xlabel('$\epsilon$',fontsize='15')
axes[0].text(0.01+min(axes[0].get_xlim()),0.9*max(axes[0].get_ylim()), '(a)',fontsize=22)
cbar1 = fig.colorbar(cp,ax=axes[0])
cbar1.ax.set_title('%') 
axes[0].set_xlim(right=1.25)

X2, Y2 = np.meshgrid(xlist, ylist)
Z2=(np.mean(NUM_AGENTS_INFLUENCED,axis=2)[1]/Num_ppl)*100
cp = axes[1].contourf(X2, Y2, Z2,len(MEAN_vals),cmap='YlOrRd')
axes[1].set_ylabel('$\mu$',fontsize='15')
axes[1].set_xlabel('$\epsilon$',fontsize='15')
axes[1].text(0.01+min(axes[1].get_xlim()),0.9*max(axes[1].get_ylim()), '(b)',fontsize=22)
cbar1 = fig.colorbar(cp,ax=axes[1])
cbar1.ax.set_title('%') 
axes[1].set_xlim(right=1.25)

data=np.load('PAPER_DATA/VARIATION_MU_DELTA_P/VAR_FILES/Proportion_variation.npz')
Num_ppl=data['Num_ppl']
e_range=data['e_range']
DELTA_vals=data['DELTA_vals']
MEAN_vals=data['MEAN_vals']
NUM_AGENTS_INFLUENCED=data['NUM_AGENTS_INFLUENCED']*0.1
Array_pos=data['Array_pos']

xlist=e_range
ylist=Array_pos/Num_ppl
X3, Y3 = np.meshgrid(xlist, ylist)
Z3=(np.mean(NUM_AGENTS_INFLUENCED,axis=2)[0]/Num_ppl)*100
cp = axes[2].contourf(X3, Y3, Z3,len(Array_pos),cmap='YlOrRd')
axes[2].set_ylabel('$p$',fontsize='15')
axes[2].set_xlabel('$\epsilon$',fontsize='15')
axes[2].text(0.01+min(axes[2].get_xlim()),0.99*max(axes[2].get_ylim()), '(c)',fontsize=22)
cbar1 = fig.colorbar(cp,ax=axes[2])
cbar1.ax.set_title('%') 
axes[2].set_xlim(right=1.25)

#plt.savefig('/home/glory/Dropbox/Opinion Dynamics/paper/Figures/Variation.pdf', bbox_inches='tight')


# In[37]:
####################################################

####################################################

import numpy as np
import matplotlib.pyplot as plt

NP_district=101

path='PAPER_DATA/BI-PARTY/'
data=np.load(path+'Multi_States/Initial_considerations.npz', allow_pickle=True)
NUM_STATES=data['NUM_STATES']
NUM_SEATS_PER_STATE=data['NUM_SEATS_PER_STATE']
Num_simulation=data['Num_simulation']
Num_type=3

data=np.load(path+'US/RealDout_'+str(1)+'_robustness_'+str(0)+'_.npz')
e_range=data['e_range']
neps=len(e_range)
EFF_SR=np.zeros((len(NUM_STATES),Num_type,Num_simulation,neps))
EFF_PR=np.zeros((len(NUM_STATES),Num_type,Num_simulation,neps))
EFF_WTAR=np.zeros((len(NUM_STATES),Num_type,Num_simulation,neps))

NEFF_SR=np.zeros((len(NUM_STATES),Num_type,Num_simulation,neps))
NEFF_PR=np.zeros((len(NUM_STATES),Num_type,Num_simulation,neps))
NEFF_WTAR=np.zeros((len(NUM_STATES),Num_type,Num_simulation,neps))

PEFF_SR=np.zeros((len(NUM_STATES),Num_type,Num_simulation,neps))
PEFF_PR=np.zeros((len(NUM_STATES),Num_type,Num_simulation,neps))
PEFF_WTAR=np.zeros((len(NUM_STATES),Num_type,Num_simulation,neps))

P_SR=np.zeros((len(NUM_STATES)*Num_simulation,Num_type,neps))
P_PR=np.zeros((len(NUM_STATES)*Num_simulation,Num_type,neps))
P_WTAR=np.zeros((len(NUM_STATES)*Num_simulation,Num_type,neps))


AGENTS=[]
Agents_distribution=[]

for N_S in range(len(NUM_STATES)):
    data=np.load(path+'US/RealDout_'+str(1)+'_robustness_'+str(N_S)+'_.npz')
    Num_districts_per_state=data['Num_districts_perstate']
    Agents=Num_districts_per_state*NP_district
    Agents[Agents%2==0]+=1
    Agents_distribution.append(Agents)
    Total_agents=np.sum(Agents)
    AGENTS.append(Total_agents)

for N_S in range(len(NUM_STATES)):
    Num_states=NUM_STATES[N_S]
    Num_districts_perstate=NUM_SEATS_PER_STATE[N_S]
    Num_districts=np.sum(Num_districts_perstate).astype('int')
    Agents=Num_districts_perstate*NP_district
    Agents[Agents%2==0]+=1
    Proportion=(Num_districts_perstate/Agents) 
    
    for num_type in range(Num_type):
        data=np.load(path+'US/RealDout_'+str(num_type+1)+'_robustness_'+str(N_S)+'_.npz')
        EFF_SR[N_S,num_type]=data['E_SR']
        EFF_PR[N_S,num_type]=data['E_PR']
        EFF_WTAR[N_S,num_type]=data['E_WTAR']
        
        NEFF_SR[N_S,num_type]=(data['E_SR']*0.1)/AGENTS[N_S]
        NEFF_PR[N_S,num_type]=(data['E_PR']*0.1)/AGENTS[N_S]
        NEFF_WTAR[N_S,num_type]=(data['E_WTAR']*0.1)/AGENTS[N_S]
        
        PEFF_SR[N_S,num_type]=(data['E_SR']/AGENTS[N_S])*100
        PEFF_PR[N_S,num_type]=(data['E_PR']/AGENTS[N_S])*100
        PEFF_WTAR[N_S,num_type]=(data['E_WTAR']/AGENTS[N_S])*100
        
        P_SR[N_S*Num_simulation:N_S*Num_simulation+Num_simulation,num_type]=(data['E_SR']/AGENTS[N_S])*100
        P_PR[N_S*Num_simulation:N_S*Num_simulation+Num_simulation,num_type]=(data['E_PR']/AGENTS[N_S])*100
        P_WTAR[N_S*Num_simulation:N_S*Num_simulation+Num_simulation,num_type]=(data['E_WTAR']/AGENTS[N_S])*100
        
AD=[]
AG=[]
for N_S in range(len(NUM_STATES)):
    AG.append(np.sum(Agents_distribution[N_S]))
    for item in Agents_distribution[N_S]:
        AD.append(item)
  
fig,axes=plt.subplots(1,Num_type, figsize=(25,5))   
#fig.suptitle('Minimum '+str(min(AD))+', Maximum '+str(max(AD))+', Average '+str(np.mean(np.array(AD)).round(decimals=3))+'      Total N max:'+str(max(AG))+', Total N min:'+str(min(AG))+', Total N mean:'+str(np.mean(np.array(AG)).round(decimals=3)))
ABC=['(aa)','(b)','(c)']
for num_type in range(Num_type):
    axes[num_type].plot(e_range,np.mean(np.mean(PEFF_PR,axis=2),axis=0)[num_type],color='blue',label='PR')
    axes[num_type].plot(e_range,np.mean(np.mean(PEFF_WTAR,axis=2),axis=0)[num_type],color='green',label='WTAR')
    axes[num_type].plot(e_range,np.mean(np.mean(PEFF_SR,axis=2),axis=0)[num_type],color='red',label='SR')
    axes[num_type].legend(loc='upper right')
    axes[num_type].set_xlabel('$\epsilon$',fontsize=15)
    axes[num_type].set_ylabel('${\%}$',fontsize=15)
    axes[num_type].text(0.1*min(axes[num_type].get_xlim()),0.9*max(axes[num_type].get_ylim()), ABC[num_type],fontsize=20)
plt.title('Effort vs. epsilon, percentage of agents')

#plt.savefig('/home/glory/Dropbox/Opinion Dynamics/paper/Figures/Fig_3g.pdf',bbox_inches='tight')         


# In[38]:
'''

fig,axes=plt.subplots(1,Num_type, figsize=(25,5))   
#fig.suptitle('Minimum '+str(min(AD))+', Maximum '+str(max(AD))+', Average '+str(np.mean(np.array(AD)).round(decimals=3))+'      Total N max:'+str(max(AG))+', Total N min:'+str(min(AG))+', Total N mean:'+str(np.mean(np.array(AG)).round(decimals=3)))
ABC=['(a)','(bb)','(c)']
for num_type in range(Num_type):
    axes[num_type].plot(e_range,np.mean(P_PR,axis=0)[num_type],color='blue',label='PR')
    axes[num_type].plot(e_range,np.mean(P_WTAR,axis=0)[num_type],color='green',label='WTAR')
    axes[num_type].plot(e_range,np.mean(P_SR,axis=0)[num_type],color='red',label='SR')
    axes[num_type].legend(loc='upper right')
    axes[num_type].set_xlabel('$\epsilon$',fontsize=15)
    axes[num_type].set_ylabel('${\%}$',fontsize=15)
    axes[num_type].text(0.1*min(axes[num_type].get_xlim()),0.9*max(axes[num_type].get_ylim()), ABC[num_type],fontsize=20)
plt.title('Effort vs. epsilon, percentage of agents')
''' 


# In[14]:


fig,axes=plt.subplots(1,Num_type, figsize=(25,5))   
#fig.suptitle('Minimum '+str(min(AD))+', Maximum '+str(max(AD))+', Average '+str(np.mean(np.array(AD)).round(decimals=3))+'      Total N max:'+str(max(AG))+', Total N min:'+str(min(AG))+', Total N mean:'+str(np.mean(np.array(AG)).round(decimals=3)))
ABC=['(a)','(b)','(cc)']
for num_type in range(Num_type):
    
    ymean = np.mean(P_PR,axis=0)[num_type]
    ymin = np.mean(P_PR,axis=0)[num_type]-np.std(P_PR,axis=0)[num_type]
    ymax = np.mean(P_PR,axis=0)[num_type]+np.std(P_PR,axis=0)[num_type]
    
    yerror = np.stack((ymean-ymin, ymax-ymean))

    axes[num_type].fill_between(e_range, ymin, ymax, color='b',alpha=0.1, label='error band')
    #axes[num_type].errorbar(e_range, ymean, yerror, errorevery=5, color='tab:blue', ecolor='tab:blue',capsize=3, linewidth=1)
    axes[num_type].plot(e_range,np.mean(P_PR,axis=0)[num_type],color='blue',label='PR')
    
    ymean = np.mean(P_SR,axis=0)[num_type]
    ymin = np.mean(P_SR,axis=0)[num_type]-np.std(P_PR,axis=0)[num_type]
    ymax = np.mean(P_SR,axis=0)[num_type]+np.std(P_PR,axis=0)[num_type]
    axes[num_type].fill_between(e_range, ymin, ymax,color='r',alpha=0.1, label='error band')
    #axes[num_type].errorbar(e_range, ymean, yerror, color='tab:red', ecolor='tab:red',capsize=3, linewidth=1)
    axes[num_type].plot(e_range,np.mean(P_SR,axis=0)[num_type],color='red',label='SR')
    
    ymean = np.mean(P_WTAR,axis=0)[num_type]
    ymin = np.mean(P_WTAR,axis=0)[num_type]-np.std(P_PR,axis=0)[num_type]
    ymax = np.mean(P_WTAR,axis=0)[num_type]+np.std(P_PR,axis=0)[num_type]
    axes[num_type].fill_between(e_range, ymin, ymax,color='g',alpha=0.1, label='error band')
    #axes[num_type].errorbar(e_range, ymean, yerror, color='tab:green', ecolor='tab:green', capsize=3, linewidth=1)
    axes[num_type].plot(e_range,np.mean(P_WTAR,axis=0)[num_type],color='green',label='WTAR')
    axes[num_type].set_xlabel('$\epsilon$',fontsize=15)
    axes[num_type].set_ylabel('${\%}$',fontsize=15)
    axes[num_type].text(0.1*min(axes[num_type].get_xlim()),0.9*max(axes[num_type].get_ylim()), ABC[num_type],fontsize=20)
   
plt.title('Effort vs. epsilon, percentage of agents and std.')


# In[15]:


fig,axes=plt.subplots(1,Num_type, figsize=(25,5))   
#fig.suptitle('Minimum '+str(min(AD))+', Maximum '+str(max(AD))+', Average '+str(np.mean(np.array(AD)).round(decimals=3))+'      Total N max:'+str(max(AG))+', Total N min:'+str(min(AG))+', Total N mean:'+str(np.mean(np.array(AG)).round(decimals=3)))
ABC=['(aa)','(bb)','(c)']
for num_type in range(Num_type):
    ymean = np.mean(P_PR,axis=0)[num_type]
    ymin = np.min(P_PR,axis=0)[num_type]
    ymax = np.max(P_PR,axis=0)[num_type]
    yerror = np.stack((ymean-ymin, ymax-ymean))
    
    axes[num_type].fill_between(e_range, ymin, ymax, alpha=0.2)
    axes[num_type].plot(e_range,ymean,color='blue',linewidth=2)
    axes[num_type].errorbar(e_range, ymean, yerror, color='tab:blue', ecolor='tab:blue',capsize=3, linewidth=1)

    ymean = np.mean(P_SR,axis=0)[num_type]
    ymin = np.min(P_SR,axis=0)[num_type]
    ymax = np.max(P_SR,axis=0)[num_type]
    yerror = np.stack((ymean-ymin, ymax-ymean))

    axes[num_type].fill_between(e_range, ymin, ymax, alpha=0.2)
    axes[num_type].plot(e_range,ymean,color='red',linewidth=2)
    axes[num_type].errorbar(e_range, ymean, yerror, color='tab:red', ecolor='tab:red',capsize=3, linewidth=1)
    
    ymean = np.mean(P_WTAR,axis=0)[num_type]
    ymin = np.min(P_WTAR,axis=0)[num_type]
    ymax = np.max(P_WTAR,axis=0)[num_type]
    yerror = np.stack((ymean-ymin, ymax-ymean))

    axes[num_type].fill_between(e_range, ymin, ymax, alpha=0.2)
    axes[num_type].plot(e_range,ymean,color='green',linewidth=2)
    axes[num_type].errorbar(e_range, ymean, yerror, color='tab:green', ecolor='tab:green',capsize=3, linewidth=1)
    
    axes[num_type].set_xlabel('$\epsilon$',fontsize=15)
    axes[num_type].set_ylabel('${\%}$',fontsize=15)
    axes[num_type].text(0.1*min(axes[num_type].get_xlim()),0.9*max(axes[num_type].get_ylim()), ABC[num_type],fontsize=20)
   
plt.title('Effort vs. epsilon, percentage of agents, min-max values')


# In[16]:


fig,axes=plt.subplots(1,Num_type, figsize=(25,5))   
#fig.suptitle('Minimum '+str(min(AD))+', Maximum '+str(max(AD))+', Average '+str(np.mean(np.array(AD)).round(decimals=3))+'      Total N max:'+str(max(AG))+', Total N min:'+str(min(AG))+', Total N mean:'+str(np.mean(np.array(AG)).round(decimals=3)))
ABC=['(aa)','(b)','(cc)']
for num_type in range(Num_type):
    ymean = np.mean(P_PR,axis=0)[num_type]
    ymin = ymean-np.percentile(P_PR,50,axis=0)[num_type]
    ymax = ymean+np.percentile(P_PR,50,axis=0)[num_type]
    yerror = np.stack((ymean-ymin, ymax-ymean))
    axes[num_type].fill_between(e_range, ymin, ymax, alpha=0.2)
    axes[num_type].plot(e_range,ymean,color='blue',linewidth=2)
    axes[num_type].errorbar(e_range, ymean, yerror, color='tab:blue', ecolor='tab:blue',capsize=3, linewidth=1)
    
    ymean = np.mean(P_SR,axis=0)[num_type]
    ymin = ymean-np.percentile(P_SR,50,axis=0)[num_type]
    ymax = ymean+np.percentile(P_SR,50,axis=0)[num_type]
    yerror = np.stack((ymean-ymin, ymax-ymean))
    axes[num_type].fill_between(e_range, ymin, ymax, alpha=0.2)
    axes[num_type].plot(e_range,ymean,color='red',linewidth=2)
    axes[num_type].errorbar(e_range, ymean, yerror, color='tab:blue', ecolor='tab:blue',capsize=3, linewidth=1)
    
    ymean = np.mean(P_WTAR,axis=0)[num_type]
    ymin = ymean-np.percentile(P_WTAR,50,axis=0)[num_type]
    ymax = ymean+np.percentile(P_WTAR,50,axis=0)[num_type]
    yerror = np.stack((ymean-ymin, ymax-ymean))
    axes[num_type].fill_between(e_range, ymin, ymax, alpha=0.2)
    axes[num_type].plot(e_range,ymean,color='green',linewidth=2)
    axes[num_type].errorbar(e_range, ymean, yerror, color='tab:blue', ecolor='tab:blue',capsize=3, linewidth=1)
    
    axes[num_type].set_xlabel('$\epsilon$',fontsize=15)
    axes[num_type].set_ylabel('${\%}$',fontsize=15)
    axes[num_type].text(0.1*min(axes[num_type].get_xlim()),0.9*max(axes[num_type].get_ylim()), ABC[num_type],fontsize=20)
    
plt.title('Effort vs. epsilon, percentage of agents and 50-percentiles')

# In[72]:

#######################################################

P_PRn = np.zeros((1500,3,25))
P_WTARn = np.zeros((1500,3,25))
P_SRn = np.zeros((1500,3,25))
for i in range(1500):
    for j in range(3):
        for k in range(25):
            P_PRn[i,j,k] = 1.
            P_WTARn[i,j,k] = P_WTAR[i,j,k]/P_PR[i,j,k]
            P_SRn[i,j,k] = P_SR[i,j,k]/P_PR[i,j,k]

RANKING=np.zeros((Num_type,Num_simulation*len(NUM_STATES)*len(e_range),3))
for num_type in range(Num_type):
    num=0
    for num_sim in range(Num_simulation*len(NUM_STATES)):
        for num_eps in range(len(e_range)):
            Arr=np.zeros((3))
            Arr[0],Arr[1],Arr[2]=P_PR[num_sim,num_type,num_eps],P_WTAR[num_sim,num_type,num_eps],P_SR[num_sim,num_type,num_eps]
            RANKING[num_type,num]=np.argsort(-Arr).argsort()
            num+=1

fig,axes=plt.subplots(1,Num_type,figsize=(20,5))

for num_type in range(Num_type):
    axes[num_type].hist(RANKING[num_type]+1)

##############################################################################
from mpl_toolkits.mplot3d import Axes3D

fig = plt.figure()
ax = fig.add_subplot(111,projection='3d')
nbins = 100
for z in range(25):
    y1s = RANKING[0,range(z*1500,(z+1)*1500),0]
    y2s = RANKING[0,range(z*1500,(z+1)*1500),1]
    y3s = RANKING[0,range(z*1500,(z+1)*1500),2]
    hist1, bins1 = np.histogram(y1s,bins=nbins)
    hist2, bins2 = np.histogram(y2s,bins=nbins)
    hist3, bins3 = np.histogram(y3s,bins=nbins)
    x1s = (bins1[:-1] + bins1[1:])/2
    x2s = (bins2[:-1] + bins2[1:])/2
    x3s = (bins3[:-1] + bins3[1:])/2
    ax.bar(x1s,hist1,zs=z,color='C0',zdir='y',alpha=.8)
    ax.bar(x2s,hist2,zs=z,color='C1',zdir='y',alpha=.8)
    ax.bar(x3s,hist3,zs=z,color='C2',zdir='y',alpha=.8)
plt.show()

########################################################################
# In[94]:


NP_SR=np.zeros((Num_simulation*len(NUM_STATES),Num_type,len(e_range)))
NP_WTAR=np.zeros((Num_simulation*len(NUM_STATES),Num_type,len(e_range)))
NP_PR=np.zeros((Num_simulation*len(NUM_STATES),Num_type,len(e_range)))

for num_type in range(Num_type):
    for num_eps in range(len(e_range)):
        for num_sim in range(Num_simulation*len(NUM_STATES)):
            Arr=np.zeros((3))
            Arr[0],Arr[1],Arr[2]=P_PR[num_sim,num_type,num_eps],P_WTAR[num_sim,num_type,num_eps],P_SR[num_sim,num_type,num_eps]
            NP_PR[num_sim,num_type,num_eps]=P_PR[num_sim,num_type,num_eps]/np.max(Arr)
            NP_SR[num_sim,num_type,num_eps]=P_SR[num_sim,num_type,num_eps]/np.max(Arr)
            NP_WTAR[num_sim,num_type,num_eps]=P_WTAR[num_sim,num_type,num_eps]/np.max(Arr)
            


# In[96]:


fig,axes=plt.subplots(1,Num_type,figsize=(20,5))
fig.suptitle('Normalized with maximum of (SR, PR, WTAR) effort')
for num_type in range(Num_type):
    '''
    axes[num_type].plot(e_range, np.mean(NP_PR,axis=0)[num_type],color='blue',label='PR')
    axes[num_type].plot(e_range, np.mean(NP_WTAR,axis=0)[num_type],color='green',label='WTAR')
    axes[num_type].plot(e_range, np.mean(NP_SR,axis=0)[num_type],color='red',label='SR')
    axes[num_type].set_xlabel('$\epsilon$')
    axes[num_type].legend()
    axes[num_type].set_ylabel('Normalised %')
    ''' 
    axes[num_type].plot(e_range, np.median(P_PRn,axis=0)[num_type],color='blue',label='PR')
    axes[num_type].plot(e_range, np.median(P_WTARn,axis=0)[num_type],color='green',label='WTAR')
    axes[num_type].plot(e_range, np.percentile(P_WTARn,40,axis=0)[num_type],'--g')
    axes[num_type].plot(e_range, np.percentile(P_WTARn,60,axis=0)[num_type],'--g')
    axes[num_type].plot(e_range, np.median(P_SRn,axis=0)[num_type],color='red',label='SR')
    axes[num_type].plot(e_range, np.percentile(P_SRn,25,axis=0)[num_type],'--r')
    axes[num_type].plot(e_range, np.percentile(P_SRn,75,axis=0)[num_type],'--r')
    axes[num_type].set_xlabel('epsilon')
    axes[num_type].set_ylabel('relative effort')
    axes[num_type].legend()


# In[43]:


fig,axes=plt.subplots(1,Num_type, figsize=(25,5))   
#fig.suptitle('Minimum '+str(min(AD))+', Maximum '+str(max(AD))+', Average '+str(np.mean(np.array(AD)).round(decimals=3))+'      Total N max:'+str(max(AG))+', Total N min:'+str(min(AG))+', Total N mean:'+str(np.mean(np.array(AG)).round(decimals=3)))
ABC=['(a)','(b)','(c)']
for num_type in range(Num_type):
    axes[num_type].plot(e_range,np.median(P_PR,axis=0)[num_type],color='blue',label='PR')
    axes[num_type].plot(e_range,np.median(P_WTAR,axis=0)[num_type],color='green',label='WTAR')
    axes[num_type].plot(e_range,np.median(P_SR,axis=0)[num_type],color='red',label='SR')
    axes[num_type].legend(loc='upper right')
    axes[num_type].set_xlabel('$\epsilon$',fontsize=15)
    axes[num_type].set_ylabel('${\%}$',fontsize=15)
    axes[num_type].text(0.1*min(axes[num_type].get_xlim()),0.9*max(axes[num_type].get_ylim()), ABC[num_type],fontsize=20)
    


# In[17]:


path='PAPER_DATA/MULTI-PARTY/EPS_CRITICAL/'

BIAS_VAL=np.linspace(0,5,11)*0.05
BIAS_VAL=BIAS_VAL[0:1]
Num_simulation=500
neps=151

Total_agents=2001

color=plt.get_cmap("tab10")

Num_party=2
EFFORT_2=np.zeros((len(BIAS_VAL),Num_simulation, neps, Num_party-1))
NUM_VOTES_2=np.zeros((len(BIAS_VAL),Num_simulation,neps,Num_party))
PP2=np.zeros((len(BIAS_VAL),neps,Num_party))
NPP2=np.zeros((len(BIAS_VAL),neps,Num_party))
PP2_2=np.zeros((len(BIAS_VAL),neps,Num_party))
NPP2_2=np.zeros((len(BIAS_VAL),neps,Num_party))
N_VOTES_2=np.zeros((len(BIAS_VAL),Num_simulation,neps,Num_party))
NUM_AGENTS_INFLUENCED_2=np.zeros((len(BIAS_VAL),Num_simulation,neps,Num_party-1))
for b_val in range(len(BIAS_VAL)):
    data=np.load(path+'Ord_all_with_eq_area_'+str(b_val)+'_'+str(Num_party)+'.npz')
    e_range=data['e_range']
    NUM_VOTES_2[b_val]=data['NUM_VOTES']
    Spread_cent_min_2=data['Spread_cent_min']
    EFFORT_2[b_val]=data['Effort']
    NUM_AGENTS_INFLUENCED_2[b_val]=data['NUM_AGENTS_INFLUENCED']
    PP2[b_val]=((np.sum((np.argsort(data['NUM_VOTES']).argsort()==Num_party-1)*1,axis=0))/Num_simulation)*100  
    NPP2[b_val]=(np.sum((np.argsort(data['NUM_VOTES']).argsort()==Num_party-1)*1,axis=0))
    PP2_2[b_val]=((np.sum((np.argsort(data['NUM_VOTES']).argsort()==Num_party-2)*1,axis=0))/Num_simulation)*100   
    NPP2_2[b_val]=(np.sum((np.argsort(data['NUM_VOTES']).argsort()==Num_party-2)*1,axis=0))
    for num in range(Num_simulation):
        for eps in range(len(e_range)):
            Num_votes=NUM_VOTES_2[b_val,num,eps]
            WParg=np.argsort(Num_votes)
            party=Num_party-2
            winner=Num_party-1
            win=int(WParg[winner])
            N_VOTES_2[b_val,num,eps,0]=Num_votes[win]
            for j in range(Num_party-1):  
                wp=int(WParg[party])
                N_VOTES_2[b_val,num,eps,j+1]=Num_votes[wp]
                party=party-1
                
Num_party=3
MAX_NEIGHBOURS_3=np.zeros((len(BIAS_VAL),Num_simulation,neps,Num_party,Num_party))
EDGE_CONN_3=np.zeros((len(BIAS_VAL),Num_simulation,neps,Num_party,Num_party))
EFFORT_3=np.zeros((len(BIAS_VAL),Num_simulation, neps, Num_party-1))
NUM_VOTES_3=np.zeros((len(BIAS_VAL),Num_simulation,neps,Num_party))
PP3=np.zeros((len(BIAS_VAL),neps,Num_party))
NPP3=np.zeros((len(BIAS_VAL),neps,Num_party))
PP3_2=np.zeros((len(BIAS_VAL),neps,Num_party))
NPP3_2=np.zeros((len(BIAS_VAL),neps,Num_party))
N_VOTES_3=np.zeros((len(BIAS_VAL),Num_simulation,neps,Num_party))
NUM_AGENTS_INFLUENCED_3=np.zeros((len(BIAS_VAL),Num_simulation,neps,Num_party-1))
eps_critical_3=np.zeros((len(BIAS_VAL)))

E_EXT_3=np.zeros((len(BIAS_VAL),neps))
E_CENT_3=np.zeros((len(BIAS_VAL),neps))

N_AG_EXT_3=np.zeros((len(BIAS_VAL),neps))
N_AG_CENT_3=np.zeros((len(BIAS_VAL),neps))

NE_EXT_3=np.zeros((len(BIAS_VAL),neps))
NE_CENT_3=np.zeros((len(BIAS_VAL),neps))


for b_val in range(len(BIAS_VAL)):
    data=np.load(path+'Ord_all_with_eq_area_'+str(b_val)+'_'+str(Num_party)+'.npz')
    e_range=data['e_range']
    MAX_NEIGHBOURS_3[b_val]=data['MAX_NEIGHBOURS']
    EDGE_CONN_3[b_val]=data['EDGE_CONN']
    Spread_cent_min_3=data['Spread_cent_min']
    NUM_VOTES_3[b_val]=data['NUM_VOTES']
    EFFORT_3[b_val]=data['Effort']
    NUM_AGENTS_INFLUENCED_3[b_val]=data['NUM_AGENTS_INFLUENCED']
    PP3[b_val]=((np.sum((np.argsort(data['NUM_VOTES']).argsort()==Num_party-1)*1,axis=0))/Num_simulation)*100 
    NPP3[b_val]=(np.sum((np.argsort(data['NUM_VOTES']).argsort()==Num_party-1)*1,axis=0))
    PP3_2[b_val]=((np.sum((np.argsort(data['NUM_VOTES']).argsort()==Num_party-2)*1,axis=0))/Num_simulation)*100 
    PP3_2[b_val]=(np.sum((np.argsort(data['NUM_VOTES']).argsort()==Num_party-2)*1,axis=0))
    
    for num in range(Num_simulation):
        for eps in range(len(e_range)):
            Num_votes=NUM_VOTES_3[b_val,num,eps]
            WParg=np.argsort(Num_votes)
            party=Num_party-2
            winner=Num_party-1
            win=int(WParg[winner])
            N_VOTES_3[b_val,num,eps,0]=Num_votes[win]
            
            wp=int(WParg[party])
            if wp ==0 or wp==Num_party-1:
                E_EXT_3[b_val,eps]+=EFFORT_3[b_val,num,eps,0]
                N_AG_EXT_3[b_val,eps]+=NUM_AGENTS_INFLUENCED_3[b_val,num,eps,0]
                NE_EXT_3[b_val,eps]+=1
            elif wp==1:
                E_CENT_3[b_val,eps]+=EFFORT_3[b_val,num,eps,0]
                N_AG_EXT_3[b_val,eps]+=NUM_AGENTS_INFLUENCED_3[b_val,num,eps,0]
                NE_CENT_3[b_val,eps]+=1
            
            for j in range(Num_party-1):  
                wp=int(WParg[party])
                N_VOTES_3[b_val,num,eps,j+1]=Num_votes[wp]
                party=party-1
    E_EXT_3[b_val]=np.nan_to_num(E_EXT_3[b_val]/NE_EXT_3[b_val])
    E_CENT_3[b_val]=np.nan_to_num(E_CENT_3[b_val]/NE_CENT_3[b_val])
    
    N_AG_EXT_3[b_val]=(np.nan_to_num(N_AG_EXT_3[b_val]/NE_EXT_3[b_val])/Total_agents)*100
    N_AG_CENT_3[b_val]=(np.nan_to_num(N_AG_CENT_3[b_val]/NE_CENT_3[b_val])/Total_agents)*100
    
    
for eps in range(len(e_range)):
    if NPP3_2[b_val,eps,0]+NPP3_2[b_val,eps,Num_party-1]==0:
        eps_critical_3[b_val]=e_range[eps-1]    
        break
        
Num_party=4
MAX_NEIGHBOURS_4=np.zeros((len(BIAS_VAL),Num_simulation,neps,Num_party,Num_party))
EDGE_CONN_4=np.zeros((len(BIAS_VAL),Num_simulation,neps,Num_party,Num_party))
EFFORT_4=np.zeros((len(BIAS_VAL),Num_simulation, neps, Num_party-1))
NUM_VOTES_4=np.zeros((len(BIAS_VAL),Num_simulation,neps,Num_party))
PP4=np.zeros((len(BIAS_VAL),neps,Num_party))
NPP4=np.zeros((len(BIAS_VAL),neps,Num_party))
PP4_2=np.zeros((len(BIAS_VAL),neps,Num_party))
NPP4_2=np.zeros((len(BIAS_VAL),neps,Num_party))

N_VOTES_4=np.zeros((len(BIAS_VAL),Num_simulation,neps,Num_party))
NUM_AGENTS_INFLUENCED_4=np.zeros((len(BIAS_VAL),Num_simulation,neps,Num_party-1))
eps_critical_4=np.zeros((len(BIAS_VAL)))

E_EXT_4=np.zeros((len(BIAS_VAL),neps))
E_CENT_4=np.zeros((len(BIAS_VAL),neps))

N_AG_EXT_4=np.zeros((len(BIAS_VAL),neps))
N_AG_CENT_4=np.zeros((len(BIAS_VAL),neps))

NE_EXT_4=np.zeros((len(BIAS_VAL),neps))
NE_CENT_4=np.zeros((len(BIAS_VAL),neps))

for b_val in range(len(BIAS_VAL)):
    data=np.load(path+'Ord_all_with_eq_area_'+str(b_val)+'_'+str(Num_party)+'.npz')
    e_range=data['e_range']
    MAX_NEIGHBOURS_4[b_val]=data['MAX_NEIGHBOURS']
    EDGE_CONN_4[b_val]=data['EDGE_CONN']   
    NUM_VOTES_4[b_val]=data['NUM_VOTES']
    Spread_cent_min_4=data['Spread_cent_min']
    EFFORT_4[b_val]=data['Effort']
    NUM_AGENTS_INFLUENCED_4[b_val]=data['NUM_AGENTS_INFLUENCED']
    PP4[b_val]=((np.sum((np.argsort(data['NUM_VOTES']).argsort()==Num_party-1)*1,axis=0))/Num_simulation)*100  
    NPP4[b_val]=(np.sum((np.argsort(data['NUM_VOTES']).argsort()==Num_party-1)*1,axis=0))
    PP4_2[b_val]=((np.sum((np.argsort(data['NUM_VOTES']).argsort()==Num_party-2)*1,axis=0))/Num_simulation)*100 
    NPP4_2[b_val]=(np.sum((np.argsort(data['NUM_VOTES']).argsort()==Num_party-2)*1,axis=0))
    for num in range(Num_simulation):
        for eps in range(len(e_range)):
            Num_votes=NUM_VOTES_4[b_val,num,eps]
            WParg=np.argsort(Num_votes)
            party=Num_party-2
            winner=Num_party-1
            win=int(WParg[winner])
            N_VOTES_4[b_val,num,eps,0]=Num_votes[win]
            
            wp=int(WParg[party])
            if wp ==0 or wp==Num_party-1:
                E_EXT_4[b_val,eps]+=EFFORT_4[b_val,num,eps,0]
                N_AG_EXT_4[b_val,eps]+=NUM_AGENTS_INFLUENCED_4[b_val,num,eps,0]
                NE_EXT_4[b_val,eps]+=1
            elif wp==1 or wp==Num_party-2:
                E_CENT_4[b_val,eps]+=EFFORT_4[b_val,num,eps,0]
                N_AG_CENT_4[b_val,eps]+=NUM_AGENTS_INFLUENCED_4[b_val,num,eps,0]
                NE_CENT_4[b_val,eps]+=1
            
            for j in range(Num_party-1):  
                wp=int(WParg[party])
                N_VOTES_4[b_val,num,eps,j+1]=Num_votes[wp]
                party=party-1
    
    E_EXT_4[b_val]=np.nan_to_num(E_EXT_4[b_val]/NE_EXT_4[b_val])
    E_CENT_4[b_val]=np.nan_to_num(E_CENT_4[b_val]/NE_CENT_4[b_val])
    
    N_AG_EXT_4[b_val]=(np.nan_to_num(N_AG_EXT_4[b_val]/NE_EXT_4[b_val])/Total_agents)*100
    N_AG_CENT_4[b_val]=(np.nan_to_num(N_AG_CENT_4[b_val]/NE_CENT_4[b_val])/Total_agents)*100
              
for eps in range(len(e_range)-1,0,-1):
    if NPP4_2[b_val,eps,0]+NPP4_2[b_val,eps,Num_party-1]!=0:
        eps_critical_4[b_val]=e_range[eps]    
        break

Num_party=5
MAX_NEIGHBOURS_5=np.zeros((len(BIAS_VAL),Num_simulation,neps,Num_party,Num_party))
EDGE_CONN_5=np.zeros((len(BIAS_VAL),Num_simulation,neps,Num_party,Num_party))
EFFORT_5=np.zeros((len(BIAS_VAL),Num_simulation, neps, Num_party-1))
NUM_VOTES_5=np.zeros((len(BIAS_VAL),Num_simulation,neps,Num_party))
PP5=np.zeros((len(BIAS_VAL),neps,Num_party))
NPP5=np.zeros((len(BIAS_VAL),neps,Num_party))
PP5_2=np.zeros((len(BIAS_VAL),neps,Num_party))
NPP5_2=np.zeros((len(BIAS_VAL),neps,Num_party))
N_VOTES_5=np.zeros((len(BIAS_VAL),Num_simulation,neps,Num_party))
NUM_AGENTS_INFLUENCED_5=np.zeros((len(BIAS_VAL),Num_simulation,neps,Num_party-1))
eps_critical_5=np.zeros((len(BIAS_VAL)))

E_EXT_5=np.zeros((len(BIAS_VAL),neps))
E_MID_5=np.zeros((len(BIAS_VAL),neps))
E_CENT_5=np.zeros((len(BIAS_VAL),neps))

N_AG_EXT_5=np.zeros((len(BIAS_VAL),neps))
N_AG_MID_5=np.zeros((len(BIAS_VAL),neps))
N_AG_CENT_5=np.zeros((len(BIAS_VAL),neps))

NE_EXT_5=np.zeros((len(BIAS_VAL),neps))
NE_MID_5=np.zeros((len(BIAS_VAL),neps))
NE_CENT_5=np.zeros((len(BIAS_VAL),neps))

for b_val in range(len(BIAS_VAL)):
    data=np.load(path+'Ord_all_with_eq_area_'+str(b_val)+'_'+str(Num_party)+'.npz')
    e_range=data['e_range']
    MAX_NEIGHBOURS_5[b_val]=data['MAX_NEIGHBOURS']
    EDGE_CONN_5[b_val]=data['EDGE_CONN']   
    NUM_VOTES_5[b_val]=data['NUM_VOTES']
    Spread_cent_min_5=data['Spread_cent_min']
    EFFORT_5[b_val]=data['Effort']
    NUM_AGENTS_INFLUENCED_5[b_val]=data['NUM_AGENTS_INFLUENCED']
    PP5[b_val]=((np.sum((np.argsort(data['NUM_VOTES']).argsort()==Num_party-1)*1,axis=0))/Num_simulation)*100 
    NPP5[b_val]=(np.sum((np.argsort(data['NUM_VOTES']).argsort()==Num_party-1)*1,axis=0))
    PP5_2[b_val]=((np.sum((np.argsort(data['NUM_VOTES']).argsort()==Num_party-2)*1,axis=0))/Num_simulation)*100  
    NPP5_2[b_val]=(np.sum((np.argsort(data['NUM_VOTES']).argsort()==Num_party-2)*1,axis=0))
    
    for num in range(Num_simulation):
        for eps in range(len(e_range)):
            Num_votes=NUM_VOTES_5[b_val,num,eps]
            WParg=np.argsort(Num_votes)
            party=Num_party-2
            winner=Num_party-1
            win=int(WParg[winner])
            N_VOTES_5[b_val,num,eps,0]=Num_votes[win]
            
            wp=int(WParg[party])
            if wp ==0 or wp==Num_party-1:
                E_EXT_5[b_val,eps]+=EFFORT_5[b_val,num,eps,0]
                N_AG_EXT_5[b_val,eps]+=NUM_AGENTS_INFLUENCED_5[b_val,num,eps,0]
                NE_EXT_5[b_val,eps]+=1
            elif wp==1 or wp==Num_party-2:
                E_MID_5[b_val,eps]+=EFFORT_5[b_val,num,eps,0]
                N_AG_MID_5[b_val,eps]+=NUM_AGENTS_INFLUENCED_5[b_val,num,eps,0]
                NE_MID_5[b_val,eps]+=1
            elif wp==2:
                E_CENT_5[b_val,eps]+=EFFORT_5[b_val,num,eps,0]
                N_AG_CENT_5[b_val,eps]+=NUM_AGENTS_INFLUENCED_5[b_val,num,eps,0]
                NE_CENT_5[b_val,eps]+=1
            
            for j in range(Num_party-1):  
                wp=int(WParg[party])
                N_VOTES_5[b_val,num,eps,j+1]=Num_votes[wp]
                party=party-1
                
    E_EXT_5[b_val]=np.nan_to_num(E_EXT_5[b_val]/NE_EXT_5[b_val])
    E_CENT_5[b_val]=np.nan_to_num(E_CENT_5[b_val]/NE_CENT_5[b_val])
    E_MID_5[b_val]=np.nan_to_num(E_MID_5[b_val]/NE_MID_5[b_val])
    
    N_AG_EXT_5[b_val]=(np.nan_to_num(N_AG_EXT_5[b_val]/NE_EXT_5[b_val])/Total_agents)*100
    N_AG_CENT_5[b_val]=(np.nan_to_num(N_AG_CENT_5[b_val]/NE_CENT_5[b_val])/Total_agents)*100
    N_AG_MID_5[b_val]=(np.nan_to_num(N_AG_MID_5[b_val]/NE_MID_5[b_val])/Total_agents)*100
    
    
for eps in range(len(e_range)-1,0,-1):
    if NPP5_2[b_val,eps,0]+NPP5_2[b_val,eps,Num_party-1]!=0:
        eps_critical_5[b_val]=e_range[eps]    
        break        
        
        
Num_party=6
MAX_NEIGHBOURS_6=np.zeros((len(BIAS_VAL),Num_simulation,neps,Num_party,Num_party))
EDGE_CONN_6=np.zeros((len(BIAS_VAL),Num_simulation,neps,Num_party,Num_party))
EFFORT_6=np.zeros((len(BIAS_VAL),Num_simulation, neps, Num_party-1))
NUM_VOTES_6=np.zeros((len(BIAS_VAL),Num_simulation,neps,Num_party))
PP6=np.zeros((len(BIAS_VAL),neps,Num_party))
NPP6=np.zeros((len(BIAS_VAL),neps,Num_party))
PP6_2=np.zeros((len(BIAS_VAL),neps,Num_party))
NPP6_2=np.zeros((len(BIAS_VAL),neps,Num_party))

E_EXT_6=np.zeros((len(BIAS_VAL),neps))
E_MID_6=np.zeros((len(BIAS_VAL),neps))
E_CENT_6=np.zeros((len(BIAS_VAL),neps))

N_AG_EXT_6=np.zeros((len(BIAS_VAL),neps))
N_AG_MID_6=np.zeros((len(BIAS_VAL),neps))
N_AG_CENT_6=np.zeros((len(BIAS_VAL),neps))

NE_EXT_6=np.zeros((len(BIAS_VAL),neps))
NE_MID_6=np.zeros((len(BIAS_VAL),neps))
NE_CENT_6=np.zeros((len(BIAS_VAL),neps))

N_VOTES_6=np.zeros((len(BIAS_VAL),Num_simulation,neps,Num_party))
NUM_AGENTS_INFLUENCED_6=np.zeros((len(BIAS_VAL),Num_simulation,neps,Num_party-1))
eps_critical_6=np.zeros((len(BIAS_VAL)))

for b_val in range(len(BIAS_VAL)):
    data=np.load(path+'Ord_all_with_eq_area_'+str(b_val)+'_'+str(Num_party)+'.npz')
    e_range=data['e_range']
    MAX_NEIGHBOURS_6[b_val]=data['MAX_NEIGHBOURS']
    EDGE_CONN_6[b_val]=data['EDGE_CONN']   
    NUM_VOTES_6[b_val]=data['NUM_VOTES']
    EFFORT_6[b_val]=data['Effort']
    NUM_AGENTS_INFLUENCED_6[b_val]=data['NUM_AGENTS_INFLUENCED']
    Spread_cent_min_6=data['Spread_cent_min']
    PP6[b_val]=((np.sum((np.argsort(data['NUM_VOTES']).argsort()==Num_party-1)*1,axis=0))/Num_simulation)*100 
    NPP6[b_val]=(np.sum((np.argsort(data['NUM_VOTES']).argsort()==Num_party-1)*1,axis=0))
    PP6_2[b_val]=((np.sum((np.argsort(data['NUM_VOTES']).argsort()==Num_party-2)*1,axis=0))/Num_simulation)*100   
    NPP6_2[b_val]=(np.sum((np.argsort(data['NUM_VOTES']).argsort()==Num_party-2)*1,axis=0))
    for num in range(Num_simulation):
        for eps in range(len(e_range)):
            Num_votes=NUM_VOTES_6[b_val,num,eps]
            WParg=np.argsort(Num_votes)
            party=Num_party-2
            winner=Num_party-1
            win=int(WParg[winner])
            N_VOTES_6[b_val,num,eps,0]=Num_votes[win]
            
            wp=int(WParg[party])
            if wp ==0 or wp==Num_party-1:
                E_EXT_6[b_val,eps]+=EFFORT_6[b_val,num,eps,0]
                N_AG_EXT_6[b_val,eps]+=NUM_AGENTS_INFLUENCED_6[b_val,num,eps,0]
                NE_EXT_6[b_val,eps]+=1
            elif wp==1 or wp==Num_party-2:
                E_MID_6[b_val,eps]+=EFFORT_6[b_val,num,eps,0]
                N_AG_MID_6[b_val,eps]+=NUM_AGENTS_INFLUENCED_6[b_val,num,eps,0]
                NE_MID_6[b_val,eps]+=1
            elif wp==2 or Num_party-3:
                E_CENT_6[b_val,eps]+=EFFORT_6[b_val,num,eps,0]
                N_AG_CENT_6[b_val,eps]+=NUM_AGENTS_INFLUENCED_6[b_val,num,eps,0]
                NE_CENT_6[b_val,eps]+=1
            
            for j in range(Num_party-1):  
                wp=int(WParg[party])
                N_VOTES_6[b_val,num,eps,j+1]=Num_votes[wp]
                party=party-1
    E_EXT_6[b_val]=np.nan_to_num(E_EXT_6[b_val]/NE_EXT_6[b_val])
    E_CENT_6[b_val]=np.nan_to_num(E_CENT_6[b_val]/NE_CENT_6[b_val])
    E_MID_6[b_val]=np.nan_to_num(E_MID_6[b_val]/NE_MID_6[b_val])
    
    N_AG_EXT_6[b_val]=(np.nan_to_num(N_AG_EXT_6[b_val]/NE_EXT_6[b_val])/Total_agents)*100
    N_AG_CENT_6[b_val]=(np.nan_to_num(N_AG_CENT_6[b_val]/NE_CENT_6[b_val])/Total_agents)*100
    N_AG_MID_6[b_val]=(np.nan_to_num(N_AG_MID_6[b_val]/NE_MID_6[b_val])/Total_agents)*100
           
for eps in range(len(e_range)-1,0,-1):
    if NPP6_2[b_val,eps,0]+NPP6_2[b_val,eps,Num_party-1]!=0:
        eps_critical_6[b_val]=e_range[eps]    
        break
      
Num_party=7
MAX_NEIGHBOURS_7=np.zeros((len(BIAS_VAL),Num_simulation,neps,Num_party,Num_party))
EDGE_CONN_7=np.zeros((len(BIAS_VAL),Num_simulation,neps,Num_party,Num_party))
EFFORT_7=np.zeros((len(BIAS_VAL),Num_simulation, neps, Num_party-1))
NUM_VOTES_7=np.zeros((len(BIAS_VAL),Num_simulation,neps,Num_party))
PP7=np.zeros((len(BIAS_VAL),neps,Num_party))
NPP7=np.zeros((len(BIAS_VAL),neps,Num_party))
PP7_2=np.zeros((len(BIAS_VAL),neps,Num_party))
NPP7_2=np.zeros((len(BIAS_VAL),neps,Num_party))
N_VOTES_7=np.zeros((len(BIAS_VAL),Num_simulation,neps,Num_party))
NUM_AGENTS_INFLUENCED_7=np.zeros((len(BIAS_VAL),Num_simulation,neps,Num_party-1))
eps_critical_7=np.zeros((len(BIAS_VAL)))

E_EXT_7=np.zeros((len(BIAS_VAL),neps))
E_MID_EXT_7=np.zeros((len(BIAS_VAL),neps))
E_MID_CENT_7=np.zeros((len(BIAS_VAL),neps))
E_CENT_7=np.zeros((len(BIAS_VAL),neps))

N_AG_EXT_7=np.zeros((len(BIAS_VAL),neps))
N_AG_MID_EXT_7=np.zeros((len(BIAS_VAL),neps))
N_AG_MID_CENT_7=np.zeros((len(BIAS_VAL),neps))
N_AG_CENT_7=np.zeros((len(BIAS_VAL),neps))

NE_EXT_7=np.zeros((len(BIAS_VAL),neps))
NE_MID_EXT_7=np.zeros((len(BIAS_VAL),neps))
NE_MID_CENT_7=np.zeros((len(BIAS_VAL),neps))
NE_CENT_7=np.zeros((len(BIAS_VAL),neps))

for b_val in range(len(BIAS_VAL)):
    data=np.load(path+'Ord_all_with_eq_area_'+str(b_val)+'_'+str(Num_party)+'.npz')
    e_range=data['e_range']
    MAX_NEIGHBOURS_7[b_val]=data['MAX_NEIGHBOURS']
    EDGE_CONN_7[b_val]=data['EDGE_CONN']   
    NUM_VOTES_7[b_val]=data['NUM_VOTES']
    EFFORT_7[b_val]=data['Effort']
    Spread_cent_min_7=data['Spread_cent_min']
    NUM_AGENTS_INFLUENCED_7[b_val]=data['NUM_AGENTS_INFLUENCED']
    PP7[b_val]=((np.sum((np.argsort(data['NUM_VOTES']).argsort()==Num_party-1)*1,axis=0))/Num_simulation)*100   
    NPP7[b_val]=(np.sum((np.argsort(data['NUM_VOTES']).argsort()==Num_party-1)*1,axis=0)) 
    PP7_2[b_val]=((np.sum((np.argsort(data['NUM_VOTES']).argsort()==Num_party-2)*1,axis=0))/Num_simulation)*100   
    NPP7_2[b_val]=((np.sum((np.argsort(data['NUM_VOTES']).argsort()==Num_party-2)*1,axis=0))/Num_simulation)*100   
    for num in range(Num_simulation):
        for eps in range(len(e_range)):
            Num_votes=NUM_VOTES_7[b_val,num,eps]
            WParg=np.argsort(Num_votes)
            party=Num_party-2
            winner=Num_party-1
            win=int(WParg[winner])
            N_VOTES_7[b_val,num,eps,0]=Num_votes[win]
            
            wp=int(WParg[party])
            if wp ==0 or wp==Num_party-1:
                E_EXT_7[b_val,eps]+=EFFORT_7[b_val,num,eps,0]
                N_AG_EXT_7[b_val,eps]+=NUM_AGENTS_INFLUENCED_7[b_val,num,eps,0]
                NE_EXT_7[b_val,eps]+=1
            elif wp==1 or wp==Num_party-2:
                E_MID_EXT_7[b_val,eps]+=EFFORT_7[b_val,num,eps,0]
                N_AG_MID_EXT_7[b_val,eps]+=NUM_AGENTS_INFLUENCED_7[b_val,num,eps,0]
                NE_MID_EXT_7[b_val,eps]+=1
            elif wp==2 or Num_party-3:
                E_MID_CENT_7[b_val,eps]+=EFFORT_7[b_val,num,eps,0]
                N_AG_MID_CENT_7[b_val,eps]+=NUM_AGENTS_INFLUENCED_7[b_val,num,eps,0]
                NE_MID_CENT_7[b_val,eps]+=1
            elif wp==3:
                E_CENT_7[b_val,eps]+=EFFORT_7[b_val,num,eps,0]
                N_AG_CENT_7[b_val,eps]+=NUM_AGENTS_INFLUENCED_7[b_val,num,eps,0]
                NE_CENT_7[b_val,eps]+=1
            
            for j in range(Num_party-1):  
                wp=int(WParg[party])
                N_VOTES_7[b_val,num,eps,j+1]=Num_votes[wp]
                party=party-1
    E_EXT_7[b_val]=np.nan_to_num(E_EXT_7[b_val]/NE_EXT_7[b_val])
    E_CENT_7[b_val]=np.nan_to_num(E_CENT_7[b_val]/NE_CENT_7[b_val])
    E_MID_EXT_7[b_val]=np.nan_to_num(E_MID_EXT_7[b_val]/NE_MID_EXT_7[b_val])
    E_MID_CENT_7[b_val]=np.nan_to_num(E_MID_CENT_7[b_val]/NE_MID_CENT_7[b_val])

    N_AG_EXT_7[b_val]=(np.nan_to_num(N_AG_EXT_7[b_val]/NE_EXT_7[b_val])/Total_agents)*100
    N_AG_CENT_7[b_val]=(np.nan_to_num(N_AG_CENT_7[b_val]/NE_CENT_7[b_val])/Total_agents)*100
    N_AG_MID_EXT_7[b_val]=(np.nan_to_num(N_AG_MID_EXT_7[b_val]/NE_MID_EXT_7[b_val])/Total_agents)*100
    N_AG_MID_CENT_7[b_val]=(np.nan_to_num(N_AG_MID_CENT_7[b_val]/NE_MID_CENT_7[b_val])/Total_agents)*100

for eps in range(len(e_range)-1,0,-1):
    if NPP7_2[b_val,eps,0]+NPP7_2[b_val,eps,Num_party-1]!=0:
        eps_critical_7[b_val]=e_range[eps]    
        break
            


# In[18]:


fig,axes=plt.subplots(2,3,figsize=(20,10))
width_stack_plot=0.006
ABC=['(a)','(b)','(c)','(d)','(e)','(ff)']

i,j=0,0
axes[i,j].plot(e_range,N_AG_EXT_3[b_val]*(NE_EXT_3[b_val]/Num_simulation),linestyle='--',label='${\%}_{P_1} + {\%}_{P_3}$')
axes[i,j].plot(e_range,N_AG_CENT_3[b_val]*(NE_CENT_3[b_val]/Num_simulation),linestyle='--',label='${\%}_{P_2}$')
#axes[i,j].plot(e_range,100*(np.mean(NUM_AGENTS_INFLUENCED_3[b_val],axis=0)[:,0]/Total_agents),color='k',linewidth=2,label='$N^{\%}_{P_1} +N^{\%}_{P_2} + N^{\%}_{P_3}$')
axes[i,j].plot(e_range,100*(np.mean(NUM_AGENTS_INFLUENCED_3[b_val],axis=0)/Total_agents)[:,0],color='k',linewidth=2,label='$ \sum_{i=1}^p {\%}_{P_i}$')
axes[i,j].axvline(e_range[np.argmax(np.mean(NUM_AGENTS_INFLUENCED_3[b_val],axis=0)[:,0])],color='k')
#axes[i,j].axvline(eps_critical_3[b_val],color='m',label='$\epsilon_c$')
axes[i,j].axvline(e_range[np.argmax(N_AG_EXT_3[b_val]*(NE_EXT_3[b_val]/Num_simulation))],color='r',linestyle='--',linewidth=2)
#axes[i,j].axvline(e_range[np.argmax(E_EXT_3[b_val]*(NE_EXT_3[b_val]/Num_simulation))],color='g',linestyle='--', label='$\epsilon_{ext}(extr.)$')
axes[i,j].text(e_range[np.argmax(100*(np.mean(NUM_AGENTS_INFLUENCED_3[b_val],axis=0)[:,0]/Total_agents))], np.max(100*(np.mean(NUM_AGENTS_INFLUENCED_3[b_val],axis=0)[:,0]/Total_agents)),'$\epsilon_m$',color='k',fontsize=15)
axes[i,j].text(e_range[np.argmax(N_AG_EXT_3[b_val]*(NE_EXT_3[b_val]/Num_simulation))],0.5*np.max(N_AG_EXT_3[b_val]*(NE_EXT_3[b_val]/Num_simulation)),'$\epsilon_{ext}$', color='r', fontsize=15)
axes[i,j].text(0.1+min(axes[i,j].get_xlim()),0.9*max(axes[i,j].get_ylim()), ABC[0],fontsize=20)
axes[i,j].set_ylabel('${\%}_{P_i}$',fontsize='15')
axes[i,j].set_xlabel('$\epsilon$',fontsize='15')
axes[i,j].legend(loc='upper right')
axes[i,j].set_xlim([0.01,1.25])

i,j=0,1
axes[i,j].plot(e_range,N_AG_EXT_4[b_val]*(NE_EXT_4[b_val]/Num_simulation),linestyle='--',label='${\%}_{P_1} +{\%}_{P_4}$')
axes[i,j].plot(e_range,N_AG_CENT_4[b_val]*(NE_CENT_4[b_val]/Num_simulation),linestyle='--',label='${\%}_{P_2} +{\%}_{P_3}$')
#axes[i,j].plot(e_range,100*(np.mean(NUM_AGENTS_INFLUENCED_4[b_val],axis=0)/Total_agents)[:,0],color='k',linewidth=2,label='$N^{\%}_{P_1} +N^{\%}_{P_2}+N^{\%}_{P_3} +N^{\%}_{P_4}$')
axes[i,j].plot(e_range,100*(np.mean(NUM_AGENTS_INFLUENCED_4[b_val],axis=0)/Total_agents)[:,0],color='k',linewidth=2,label='$ \sum_{i=1}^p {\%}_{P_i}$')
axes[i,j].axvline(e_range[np.argmax(100*(np.mean(NUM_AGENTS_INFLUENCED_4[b_val],axis=0)/Total_agents)[:,0])],color='k')
axes[i,j].axvline(eps_critical_4[b_val],color='b')
axes[i,j].axvline(e_range[np.argmax(N_AG_EXT_4[b_val]*(NE_EXT_4[b_val]/Num_simulation))],color='r',linewidth=2,linestyle='--')
#axes[i,j].axvline(e_range[np.argmax(E_EXT_5[b_val]*(NE_EXT_5[b_val]/Num_simulation))],color='g',linestyle='--', label='$\epsilon_{ext}(extr.)$')
axes[i,j].text(e_range[np.argmax(100*(np.mean(NUM_AGENTS_INFLUENCED_4[b_val],axis=0)/Total_agents)[:,0])], np.max(100*(np.mean(NUM_AGENTS_INFLUENCED_4[b_val],axis=0)/Total_agents)[:,0]),'$\epsilon_m$',color='k',fontsize=15)
axes[i,j].text(e_range[np.argmax(N_AG_EXT_4[b_val]*(NE_EXT_4[b_val]/Num_simulation))],0.5*np.max(N_AG_EXT_4[b_val]*(NE_EXT_4[b_val]/Num_simulation)),'$\epsilon_{ext}$', color='r', fontsize=15)
axes[i,j].text(eps_critical_4[b_val], 0.75*np.max(100*(np.mean(NUM_AGENTS_INFLUENCED_4[b_val],axis=0)/Total_agents)[:,0]), '$\epsilon_c$', color='b',fontsize=15)
axes[i,j].text(0.1+min(axes[i,j].get_xlim()),0.9*max(axes[i,j].get_ylim()), ABC[1],fontsize=20)
axes[i,j].set_ylabel('${\%}_{P_i}$',fontsize='15')
axes[i,j].set_xlabel('$\epsilon$',fontsize='15')
axes[i,j].legend()
axes[i,j].set_xlim([0.01,1.25])

i,j=0,2
axes[i,j].plot(e_range,N_AG_EXT_5[b_val]*(NE_EXT_5[b_val]/Num_simulation),linestyle='--',label='${\%}_{P_1} +{\%}_{P_5}$')
axes[i,j].plot(e_range,N_AG_MID_5[b_val]*(NE_MID_5[b_val]/Num_simulation),linestyle='--',label='${\%}_{P_2} +{\%}_{P_4}$')
axes[i,j].plot(e_range,N_AG_CENT_5[b_val]*(NE_CENT_5[b_val]/Num_simulation),linestyle='--',label='${\%}_{P_3}$')
#axes[i,j].plot(e_range,100*(np.mean(NUM_AGENTS_INFLUENCED_5[b_val],axis=0)/Total_agents)[:,0],color='k',linewidth=2,label='$N^{\%}_{P_1} +N^{\%}_{P_2}+N^{\%}_{P_3} +N^{\%}_{P_4} + N^{\%}_{P_5}$')
axes[i,j].plot(e_range,100*(np.mean(NUM_AGENTS_INFLUENCED_5[b_val],axis=0)/Total_agents)[:,0],color='k',linewidth=2,label='$ \sum_{i=1}^p {\%}_{P_i}$')
axes[i,j].axvline(e_range[np.argmax(100*(np.mean(NUM_AGENTS_INFLUENCED_5[b_val],axis=0)/Total_agents)[:,0])],color='k')
axes[i,j].axvline(eps_critical_5[b_val],color='b')
axes[i,j].axvline(e_range[np.argmax(N_AG_EXT_5[b_val]*(NE_EXT_5[b_val]/Num_simulation))],color='r',linewidth=2,linestyle='--')
axes[i,j].text(e_range[np.argmax(100*(np.mean(NUM_AGENTS_INFLUENCED_5[b_val],axis=0)/Total_agents)[:,0])], np.max(100*(np.mean(NUM_AGENTS_INFLUENCED_5[b_val],axis=0)/Total_agents)[:,0]),'$\epsilon_m$',color='k',fontsize=15)
axes[i,j].text(e_range[np.argmax(N_AG_EXT_5[b_val]*(NE_EXT_5[b_val]/Num_simulation))],0.5*np.max(N_AG_EXT_5[b_val]*(NE_EXT_5[b_val]/Num_simulation)),'$\epsilon_{ext}$', color='r', fontsize=15)
axes[i,j].text(eps_critical_5[b_val], 0.75*np.max(100*(np.mean(NUM_AGENTS_INFLUENCED_5[b_val],axis=0)/Total_agents)[:,0]), '$\epsilon_c$', color='b',fontsize=15)
axes[i,j].text(0.1+min(axes[i,j].get_xlim()),0.9*max(axes[i,j].get_ylim()), ABC[2],fontsize=20)
#axes[i,j].axvline(e_range[np.argmax(E_EXT_5[b_val]*(NE_EXT_5[b_val]/Num_simulation))],color='g',linestyle='--', label='$\epsilon_{ext}(extr.)$')
axes[i,j].set_xlabel('$\epsilon$',fontsize='15')
axes[i,j].set_ylabel('${\%}_{P_i}$',fontsize='15')
axes[i,j].legend()
axes[i,j].set_xlim([0.01,1.25])

i,j=1,0
axes[i,j].plot(e_range,N_AG_EXT_6[b_val]*(NE_EXT_6[b_val]/Num_simulation),linestyle='--',label='${\%}_{P_1} +{\%}_{P_6}$')
axes[i,j].plot(e_range,N_AG_MID_6[b_val]*(NE_MID_6[b_val]/Num_simulation),linestyle='--',label='${\%}_{P_2} +{\%}_{P_5}$')
axes[i,j].plot(e_range,N_AG_CENT_6[b_val]*(NE_CENT_6[b_val]/Num_simulation),linestyle='--',label='${\%}_{P_3} +{\%}_{P_4}$')
#axes[i,j].plot(e_range,100*(np.mean(NUM_AGENTS_INFLUENCED_6[b_val],axis=0)/Total_agents)[:,0],color='k',linewidth=2,label='$N^{\%}_{P_1} +N^{\%}_{P_2} + N^{\%}_{P_3} +N^{\%}_{P_4} + N^{\%}_{P_5} +\\xi_{P_6}$')
axes[i,j].plot(e_range,100*(np.mean(NUM_AGENTS_INFLUENCED_6[b_val],axis=0)/Total_agents)[:,0],color='k',linewidth=2,label='$ \sum_{i=1}^p {\%}_{P_i}$')
axes[i,j].axvline(e_range[np.argmax(100*(np.mean(NUM_AGENTS_INFLUENCED_6[b_val],axis=0)/Total_agents)[:,0])],color='k')
axes[i,j].axvline(eps_critical_6[b_val],color='b')
axes[i,j].axvline(e_range[np.argmax(N_AG_EXT_6[b_val]*(NE_EXT_6[b_val]/Num_simulation))],color='r',linewidth=2,linestyle='--')
axes[i,j].text(e_range[np.argmax(100*(np.mean(NUM_AGENTS_INFLUENCED_6[b_val],axis=0)/Total_agents)[:,0])], np.max(100*(np.mean(NUM_AGENTS_INFLUENCED_6[b_val],axis=0)/Total_agents)[:,0]),'$\epsilon_m$',color='k',fontsize=15)
axes[i,j].text(e_range[np.argmax(N_AG_EXT_6[b_val]*(NE_EXT_6[b_val]/Num_simulation))],0.5*np.max(N_AG_EXT_6[b_val]*(NE_EXT_6[b_val]/Num_simulation)),'$\epsilon_{ext}$', color='r', fontsize=15)
axes[i,j].text(eps_critical_6[b_val], 0.75*np.max(100*(np.mean(NUM_AGENTS_INFLUENCED_6[b_val],axis=0)/Total_agents)[:,0]), '$\epsilon_c$', color='b',fontsize=15)
axes[i,j].text(0.1+min(axes[i,j].get_xlim()),0.9*max(axes[i,j].get_ylim()), ABC[3],fontsize=20)
#axes[i,j].axvline(e_range[np.argmax(E_EXT_6[b_val]*(NE_EXT_6[b_val]/Num_simulation))],color='g',linestyle='--', label='$\epsilon_{ext}(extr.)$')
axes[i,j].set_xlabel('$\epsilon$',fontsize='15')    
axes[i,j].set_ylabel('${\%}_{P_i}$',fontsize='15')
axes[i,j].legend()
axes[i,j].set_xlim([0.01,1.25])

i,j=1,1
axes[i,j].plot(e_range,N_AG_EXT_7[b_val]*(NE_EXT_7[b_val]/Num_simulation),linestyle='--',label='${\%}_{P_1} +{\%}_{P_7}$')
axes[i,j].plot(e_range,N_AG_MID_EXT_7[b_val]*(NE_MID_EXT_7[b_val]/Num_simulation),linestyle='--',label='${\%}_{P_2} +{\%}_{P_6}$')
axes[i,j].plot(e_range,N_AG_MID_CENT_7[b_val]*(NE_MID_CENT_7[b_val]/Num_simulation),linestyle='--',label='${\%}_{P_3} +{\%}_{P_4}$')
axes[i,j].plot(e_range,N_AG_CENT_7[b_val]*(NE_CENT_7[b_val]/Num_simulation),linestyle='--',label='${\%}_{P_4}$')
#axes[i,j].plot(e_range,100*(np.mean(NUM_AGENTS_INFLUENCED_7[b_val],axis=0)/Total_agents)[:,0],color='k',linewidth=2,label='$N^{\%}_{P_1} +N^{\%}_{P_2} + N^{\%}_{P_3} +N^{\%}_{P_4} + N^{\%}_{P_5} +N^{\%}_{P_6} + N^{\%}_{P_7}$')
axes[i,j].plot(e_range,100*(np.mean(NUM_AGENTS_INFLUENCED_7[b_val],axis=0)/Total_agents)[:,0],color='k',linewidth=2,label='$ \sum_{i=1}^p {\%}_{P_i}$')
axes[i,j].axvline(e_range[np.argmax(100*(np.mean(NUM_AGENTS_INFLUENCED_7[b_val],axis=0)/Total_agents)[:,0])],color='k')
axes[i,j].axvline(eps_critical_7[b_val],color='b')
axes[i,j].axvline(e_range[np.argmax(N_AG_EXT_7[b_val]*(NE_EXT_7[b_val]/Num_simulation))],color='r',linewidth=2,linestyle='--')
axes[i,j].text(e_range[np.argmax(100*(np.mean(NUM_AGENTS_INFLUENCED_7[b_val],axis=0)/Total_agents)[:,0])], np.max(100*(np.mean(NUM_AGENTS_INFLUENCED_7[b_val],axis=0)/Total_agents)[:,0]),'$\epsilon_m$',color='k',fontsize=15)
axes[i,j].text(e_range[np.argmax(N_AG_EXT_7[b_val]*(NE_EXT_7[b_val]/Num_simulation))],0.5*np.max(E_EXT_7[b_val]*(NE_EXT_7[b_val]/Num_simulation)),'$\epsilon_{ext}$', color='r', fontsize=15)
axes[i,j].text(eps_critical_7[b_val], 0.75*np.max(100*(np.mean(NUM_AGENTS_INFLUENCED_7[b_val],axis=0)/Total_agents)[:,0]), '$\epsilon_c$', color='b',fontsize=15)
axes[i,j].text(0.1+min(axes[i,j].get_xlim()),0.9*max(axes[i,j].get_ylim()), ABC[4],fontsize=20)
#axes[i,j].axvline(e_range[np.argmax(E_EXT_6[b_val]*(NE_EXT_6[b_val]/Num_simulation))],color='g',linestyle='--', label='$\epsilon_{ext}(extr.)$')
axes[i,j].set_xlabel('$\epsilon$',fontsize='15')   
axes[i,j].set_ylabel('${\%}_{P_i}$',fontsize='15')
axes[i,j].legend()
axes[i,j].set_xlim([0.01,1.25])

i,j=1,2
Num_party=6
bottom=0
axes[i,j].bar(e_range,PP6_2[b_val,:,0],width=width_stack_plot,label='P '+str(0),color=color(0))
for num in range(Num_party-1):
    bottom=bottom+PP6_2[b_val,:,num]
    axes[i,j].bar(e_range,PP6_2[b_val,:,num+1],width=width_stack_plot,bottom=bottom,label='P '+str(num+1),color=color(num+1))
axes[i,j].set_xlabel('$\epsilon$',fontsize='15')
axes[i,j].text(0.1+min(axes[i,j].get_xlim()),0.9*max(axes[i,j].get_ylim()), ABC[5],fontsize=20)
axes[i,j].set_ylabel('% of 1st runner-up')
axes[i,j].legend(bbox_to_anchor =(0.94,0.97))
axes[i,j].set_xlim([0,1.25])

#fig.savefig('/home/glory/Dropbox/Opinion Dynamics/paper/Figures/Mult_sing_analysis_per.pdf',bbox_inches='tight')


# In[19]:


fig,axes=plt.subplots(1,5, figsize=(35,6))
ABC=['(a)','(ba)','(c)','(d)','(e)','(f)']
fig.text(0.5, 0.04,'$\epsilon$',fontsize='15', ha='center')
i=0
Num_party=3
bottom=0
axes[i].bar(e_range,PP3[b_val,:,0],width=width_stack_plot,label='P '+str(1),color=color(0))
for num in range(Num_party-1):
    bottom=bottom+PP3[b_val,:,num]
    axes[i].bar(e_range,PP3[b_val,:,num+1],width=width_stack_plot,bottom=bottom,label='P '+str(num+2),color=color(num+1))
#axes[i].set_xlabel('$\epsilon$',fontsize='15')
axes[i].text(0.1+min(axes[i].get_xlim()),0.9*max(axes[i].get_ylim()), ABC[0],fontsize=20)
axes[i].set_ylabel('% of win',fontsize='15')
axes[i].legend(bbox_to_anchor =(0.97, 0.98))
axes[i].set_xlim([0,1.25])
#axes[i].set_title('p = 3')

i=1
Num_party=4
bottom=0
axes[i].bar(e_range,PP4[b_val,:,0],width=width_stack_plot,label='P '+str(1),color=color(0))
for num in range(Num_party-1):
    bottom=bottom+PP4[b_val,:,num]
    axes[i].bar(e_range,PP4[b_val,:,num+1],width=width_stack_plot,bottom=bottom,label='P '+str(num+2),color=color(num+1))
axes[i].text(0.1+min(axes[i].get_xlim()),0.9*max(axes[i].get_ylim()), ABC[1],fontsize=20)    
#axes[i].set_xlabel('$\epsilon$',fontsize='15')
#axes[i].set_ylabel('% of win')
axes[i].legend(bbox_to_anchor =(0.97, 0.98))
axes[i].set_xlim([0,1.25])
#axes[i].set_title('p = 4')

i=2
Num_party=5
bottom=0
axes[i].bar(e_range,PP5[b_val,:,0],width=width_stack_plot,label='P '+str(1),color=color(0))
for num in range(Num_party-1):
    bottom=bottom+PP5[b_val,:,num]
    axes[i].bar(e_range,PP5[b_val,:,num+1],width=width_stack_plot,bottom=bottom,label='P '+str(num+2),color=color(num+1))
axes[i].text(0.1+min(axes[i].get_xlim()),0.9*max(axes[i].get_ylim()), ABC[2],fontsize=20)    
#axes[i].set_xlabel('$\epsilon$',fontsize='15')
#axes[i].set_ylabel('% of win')
axes[i].legend(bbox_to_anchor =(0.97, 0.98))
axes[i].set_xlim([0,1.25])
#axes[i].set_title('p = 5')

i=3
Num_party=6
bottom=0
axes[i].bar(e_range,PP6[b_val,:,0],width=width_stack_plot,label='P '+str(1),color=color(0))
for num in range(Num_party-1):
    bottom=bottom+PP6[b_val,:,num]
    axes[i].bar(e_range,PP6[b_val,:,num+1],width=width_stack_plot,bottom=bottom,label='P '+str(num+2),color=color(num+1))
axes[i].text(0.1+min(axes[i].get_xlim()),0.9*max(axes[i].get_ylim()), ABC[3],fontsize=20)    
#axes[i].set_xlabel('$\epsilon$',fontsize='15')
#axes[i].set_ylabel('% of win')
axes[i].legend(bbox_to_anchor =(0.97, 0.98))
axes[i].set_xlim([0,1.25])
#axes[i].set_title('p = 6')

i=4
Num_party=7
bottom=0
axes[i].bar(e_range,PP7[b_val,:,0],width=width_stack_plot,label='P '+str(1),color=color(0))
for num in range(Num_party-1):
    bottom=bottom+PP7[b_val,:,num]
    axes[i].bar(e_range,PP7[b_val,:,num+1],width=width_stack_plot,bottom=bottom,label='P '+str(num+2),color=color(num+1))
axes[i].text(0.1+min(axes[i].get_xlim()),0.9*max(axes[i].get_ylim()), ABC[4],fontsize=20)    
#axes[i].set_xlabel('$\epsilon$',fontsize='15')
#axes[i].set_ylabel('% of win')
axes[i].legend(bbox_to_anchor =(0.97, 0.98))
#axes[i,j].legend()
axes[i].set_xlim([0,1.25])
#axes[i].set_title('p = 7')

#plt.savefig('/home/glory/Dropbox/Opinion Dynamics/paper/Figures/Per_of_win.pdf',bbox_inches='tight')


# In[ ]:





# In[ ]:





# In[20]:


path='PAPER_DATA/MULTI-PARTY/Multi_States/'

#path='MULTI-DIMENSIONAL/HOMOGENEOUS_NO/NAT_OPN_SING_COND/Multi_States_1/'
data=np.load(path+'Initial_considerations.npz',allow_pickle=True)
val_1=0
NUM_STATES_ORIGINAL=data['NUM_STATES']
val_2=len(NUM_STATES_ORIGINAL)
change=1
NUM_STATES=data['NUM_STATES'][val_1:val_2:change]
NUM_STATES_PER_STATE=data['NUM_SEATS_PER_STATE'][val_1:val_2:change]

Num_party=3
num=0
b_val=0
data=np.load(path+'Pol_Multiparty_diff_elec_sys_'+str(Num_party)+'_'+str(num)+'.npz')
e_range=data['e_range']
BIAS_VAL=np.linspace(0,1,1)
Num_simulation=data['Num_simulation']

Num_party=2
E_SR_2=np.zeros((len(BIAS_VAL),len(NUM_STATES),Num_simulation,len(e_range),Num_party-1))
E_PR_2=np.zeros((len(BIAS_VAL),len(NUM_STATES),Num_simulation,len(e_range),Num_party-1))
E_WTAR_2=np.zeros((len(BIAS_VAL),len(NUM_STATES),Num_simulation,len(e_range),Num_party-1))

PE_SR_2=np.zeros((len(BIAS_VAL),len(NUM_STATES),Num_simulation,len(e_range),Num_party-1))
PE_PR_2=np.zeros((len(BIAS_VAL),len(NUM_STATES),Num_simulation,len(e_range),Num_party-1))
PE_WTAR_2=np.zeros((len(BIAS_VAL),len(NUM_STATES),Num_simulation,len(e_range),Num_party-1))


for b_val in range(len(BIAS_VAL)):
    num=0
    for num_state in range(val_1,val_2,change):
        data=np.load(path+'Pol_Multiparty_diff_elec_sys_'+str(Num_party)+'_'+str(num_state)+'.npz')
        e_range=data['e_range']
        Agents=NUM_STATES_PER_STATE[num]*101
        Agents[Agents%2==0]+=1
        Total_agents=np.sum(Agents)
        E_SR_2[b_val,num]=data['EFF_SR']/Total_agents
        E_PR_2[b_val,num]=data['EFF_PR']/Total_agents
        E_WTAR_2[b_val,num]=data['EFF_WTAR']/Total_agents
        
        PE_SR_2[b_val,num]=(data['NUM_AGENTS_INFLUENCED_SR']/Total_agents)*100
        PE_PR_2[b_val,num]=(data['NUM_AGENTS_INFLUENCED_PR']/Total_agents)*100
        PE_WTAR_2[b_val,num]=(data['NUM_AGENTS_INFLUENCED_WTAR']/Total_agents)*100
        
        num+=1

Num_party=3        
E_SR_3=np.zeros((len(BIAS_VAL),len(NUM_STATES),Num_simulation,len(e_range),Num_party-1))
E_PR_3=np.zeros((len(BIAS_VAL),len(NUM_STATES),Num_simulation,len(e_range),Num_party-1))
E_WTAR_3=np.zeros((len(BIAS_VAL),len(NUM_STATES),Num_simulation,len(e_range),Num_party-1))

PE_SR_3=np.zeros((len(BIAS_VAL),len(NUM_STATES),Num_simulation,len(e_range),Num_party-1))
PE_PR_3=np.zeros((len(BIAS_VAL),len(NUM_STATES),Num_simulation,len(e_range),Num_party-1))
PE_WTAR_3=np.zeros((len(BIAS_VAL),len(NUM_STATES),Num_simulation,len(e_range),Num_party-1))


for b_val in range(len(BIAS_VAL)):
    num=0
    for num_state in range(val_1,val_2,change):
        data=np.load(path+'Pol_Multiparty_diff_elec_sys_'+str(Num_party)+'_'+str(num_state)+'.npz')
        e_range=data['e_range']
        Agents=NUM_STATES_PER_STATE[num]*101
        Agents[Agents%2==0]+=1
        Total_agents=np.sum(Agents)
        E_SR_3[b_val,num]=data['EFF_SR']/Total_agents
        E_PR_3[b_val,num]=data['EFF_PR']/Total_agents
        E_WTAR_3[b_val,num]=data['EFF_WTAR']/Total_agents
       
        PE_SR_3[b_val,num]=(data['NUM_AGENTS_INFLUENCED_SR']/Total_agents)*100
        PE_PR_3[b_val,num]=(data['NUM_AGENTS_INFLUENCED_PR']/Total_agents)*100
        PE_WTAR_3[b_val,num]=(data['NUM_AGENTS_INFLUENCED_WTAR']/Total_agents)*100 
        num+=1

Num_party=4
E_SR_4=np.zeros((len(BIAS_VAL),len(NUM_STATES),Num_simulation,len(e_range),Num_party-1))
E_PR_4=np.zeros((len(BIAS_VAL),len(NUM_STATES),Num_simulation,len(e_range),Num_party-1))
E_WTAR_4=np.zeros((len(BIAS_VAL),len(NUM_STATES),Num_simulation,len(e_range),Num_party-1))

PE_SR_4=np.zeros((len(BIAS_VAL),len(NUM_STATES),Num_simulation,len(e_range),Num_party-1))
PE_PR_4=np.zeros((len(BIAS_VAL),len(NUM_STATES),Num_simulation,len(e_range),Num_party-1))
PE_WTAR_4=np.zeros((len(BIAS_VAL),len(NUM_STATES),Num_simulation,len(e_range),Num_party-1))


for b_val in range(len(BIAS_VAL)):
    num=0
    for num_state in range(val_1,val_2,change):
        data=np.load(path+'Pol_Multiparty_diff_elec_sys_'+str(Num_party)+'_'+str(num_state)+'.npz')
        e_range=data['e_range']
        Agents=NUM_STATES_PER_STATE[num]*101
        Agents[Agents%2==0]+=1
        Total_agents=np.sum(Agents)
        E_SR_4[b_val,num]=data['EFF_SR']/Total_agents
        E_PR_4[b_val,num]=data['EFF_PR']/Total_agents
        E_WTAR_4[b_val,num]=data['EFF_WTAR']/Total_agents
        
        PE_SR_4[b_val,num]=(data['NUM_AGENTS_INFLUENCED_SR']/Total_agents)*100
        PE_PR_4[b_val,num]=(data['NUM_AGENTS_INFLUENCED_PR']/Total_agents)*100
        PE_WTAR_4[b_val,num]=(data['NUM_AGENTS_INFLUENCED_WTAR']/Total_agents)*100
        num+=1

Num_party=5
E_SR_5=np.zeros((len(BIAS_VAL),len(NUM_STATES),Num_simulation,len(e_range),Num_party-1))
E_PR_5=np.zeros((len(BIAS_VAL),len(NUM_STATES),Num_simulation,len(e_range),Num_party-1))
E_WTAR_5=np.zeros((len(BIAS_VAL),len(NUM_STATES),Num_simulation,len(e_range),Num_party-1))

PE_SR_5=np.zeros((len(BIAS_VAL),len(NUM_STATES),Num_simulation,len(e_range),Num_party-1))
PE_PR_5=np.zeros((len(BIAS_VAL),len(NUM_STATES),Num_simulation,len(e_range),Num_party-1))
PE_WTAR_5=np.zeros((len(BIAS_VAL),len(NUM_STATES),Num_simulation,len(e_range),Num_party-1))

for b_val in range(len(BIAS_VAL)):
    num=0
    for num_state in range(val_1,val_2,change):
        data=np.load(path+'Pol_Multiparty_diff_elec_sys_'+str(Num_party)+'_'+str(num_state)+'.npz')
        e_range=data['e_range']
        Agents=NUM_STATES_PER_STATE[num]*101
        Agents[Agents%2==0]+=1
        Total_agents=np.sum(Agents)
        E_SR_5[b_val,num]=data['EFF_SR']/Total_agents
        E_PR_5[b_val,num]=data['EFF_PR']/Total_agents
        E_WTAR_5[b_val,num]=data['EFF_WTAR']/Total_agents
        
        PE_SR_5[b_val,num]=(data['NUM_AGENTS_INFLUENCED_SR']/Total_agents)*100
        PE_PR_5[b_val,num]=(data['NUM_AGENTS_INFLUENCED_PR']/Total_agents)*100
        PE_WTAR_5[b_val,num]=(data['NUM_AGENTS_INFLUENCED_WTAR']/Total_agents)*100
        num+=1

Num_party=6
E_SR_6=np.zeros((len(BIAS_VAL),len(NUM_STATES),Num_simulation,len(e_range),Num_party-1))
E_PR_6=np.zeros((len(BIAS_VAL),len(NUM_STATES),Num_simulation,len(e_range),Num_party-1))
E_WTAR_6=np.zeros((len(BIAS_VAL),len(NUM_STATES),Num_simulation,len(e_range),Num_party-1))

PE_SR_6=np.zeros((len(BIAS_VAL),len(NUM_STATES),Num_simulation,len(e_range),Num_party-1))
PE_PR_6=np.zeros((len(BIAS_VAL),len(NUM_STATES),Num_simulation,len(e_range),Num_party-1))
PE_WTAR_6=np.zeros((len(BIAS_VAL),len(NUM_STATES),Num_simulation,len(e_range),Num_party-1))

for b_val in range(len(BIAS_VAL)):
    num=0
    for num_state in range(val_1,val_2,change):
        data=np.load(path+'Pol_Multiparty_diff_elec_sys_'+str(Num_party)+'_'+str(num_state)+'.npz')
        e_range=data['e_range']
        Agents=NUM_STATES_PER_STATE[num]*101
        Agents[Agents%2==0]+=1
        Total_agents=np.sum(Agents)
        E_SR_6[b_val,num]=data['EFF_SR']/Total_agents
        E_PR_6[b_val,num]=data['EFF_PR']/Total_agents
        E_WTAR_6[b_val,num]=data['EFF_WTAR']/Total_agents
        
        PE_SR_6[b_val,num]=(data['NUM_AGENTS_INFLUENCED_SR']/Total_agents)*100
        PE_PR_6[b_val,num]=(data['NUM_AGENTS_INFLUENCED_PR']/Total_agents)*100
        PE_WTAR_6[b_val,num]=(data['NUM_AGENTS_INFLUENCED_WTAR']/Total_agents)*100
        num+=1


Num_party=7
E_SR_7=np.zeros((len(BIAS_VAL),len(NUM_STATES),Num_simulation,len(e_range),Num_party-1))
E_PR_7=np.zeros((len(BIAS_VAL),len(NUM_STATES),Num_simulation,len(e_range),Num_party-1))
E_WTAR_7=np.zeros((len(BIAS_VAL),len(NUM_STATES),Num_simulation,len(e_range),Num_party-1))

PE_SR_7=np.zeros((len(BIAS_VAL),len(NUM_STATES),Num_simulation,len(e_range),Num_party-1))
PE_PR_7=np.zeros((len(BIAS_VAL),len(NUM_STATES),Num_simulation,len(e_range),Num_party-1))
PE_WTAR_7=np.zeros((len(BIAS_VAL),len(NUM_STATES),Num_simulation,len(e_range),Num_party-1))

for b_val in range(len(BIAS_VAL)):
    num=0
    for num_state in range(val_1,val_2,change):
        data=np.load(path+'Pol_Multiparty_diff_elec_sys_'+str(Num_party)+'_'+str(num_state)+'.npz')
        e_range=data['e_range']
        Agents=NUM_STATES_PER_STATE[num]*101
        Agents[Agents%2==0]+=1
        Total_agents=np.sum(Agents)
        E_SR_7[b_val,num]=data['EFF_SR']/Total_agents
        E_PR_7[b_val,num]=data['EFF_PR']/Total_agents
        E_WTAR_7[b_val,num]=data['EFF_WTAR']/Total_agents
        
        PE_SR_7[b_val,num]=(data['NUM_AGENTS_INFLUENCED_SR']/Total_agents)*100
        PE_PR_7[b_val,num]=(data['NUM_AGENTS_INFLUENCED_PR']/Total_agents)*100
        PE_WTAR_7[b_val,num]=(data['NUM_AGENTS_INFLUENCED_WTAR']/Total_agents)*100
        num+=1


# In[21]:


fig,axes=plt.subplots(2,3,figsize=(20,10))
ABC=['(ab)','(b)','(c)','(d)','(e)','(f)']

i,j=0,0
axes[i,j].plot(e_range,np.mean(np.mean(PE_SR_3,axis=2),axis=1)[b_val,:,0],color='r', label='SR')
axes[i,j].plot(e_range,np.mean(np.mean(PE_PR_3,axis=2),axis=1)[b_val,:,0],color='b', label='PR')
axes[i,j].plot(e_range,np.mean(np.mean(PE_WTAR_3,axis=2),axis=1)[b_val,:,0],color='g', label='WTAR')
axes[i,j].set_xlabel('$\epsilon$',fontsize='15')
#axes[i,j].set_ylabel('Per. of agents inf.')
axes[i,j].set_ylabel('${\%}$',fontsize='15')
axes[i,j].text(0.03+min(axes[i,j].get_xlim()),0.9*max(axes[i,j].get_ylim()), ABC[0],fontsize=20)
axes[i,j].legend(loc = 'upper right')
axes[i,j].set_xlim(right=1.25)

i,j=0,1
axes[i,j].plot(e_range,np.mean(np.mean(PE_SR_4,axis=2),axis=1)[b_val,:,0],color='r', label='SR')
axes[i,j].plot(e_range,np.mean(np.mean(PE_PR_4,axis=2),axis=1)[b_val,:,0],color='b', label='PR')
axes[i,j].plot(e_range,np.mean(np.mean(PE_WTAR_4,axis=2),axis=1)[b_val,:,0],color='g', label='WTAR')
axes[i,j].text(0.03+min(axes[i,j].get_xlim()),0.91*max(axes[i,j].get_ylim()), ABC[1],fontsize=20)
axes[i,j].set_xlabel('$\epsilon$',fontsize='15')
#axes[i,j].set_ylabel('Per. of agents inf.')
axes[i,j].set_ylabel('${\%}$',fontsize='15')
axes[i,j].legend()
axes[i,j].set_xlim(right=1.25)

i,j=0,2
axes[i,j].plot(e_range,np.mean(np.mean(PE_SR_5,axis=2),axis=1)[b_val,:,0],color='r', label='SR')
axes[i,j].plot(e_range,np.mean(np.mean(PE_PR_5,axis=2),axis=1)[b_val,:,0],color='b', label='PR')
axes[i,j].plot(e_range,np.mean(np.mean(PE_WTAR_5,axis=2),axis=1)[b_val,:,0],color='g', label='WTAR')
axes[i,j].set_xlabel('$\epsilon$',fontsize='15')
#axes[i,j].set_ylabel('Per. of agents inf.')
axes[i,j].set_ylabel('${\%}$',fontsize='15')
axes[i,j].text(0.03+min(axes[i,j].get_xlim()),0.91*max(axes[i,j].get_ylim()), ABC[2],fontsize=20)
axes[i,j].legend()
axes[i,j].set_xlim(right=1.25)

i,j=1,0
axes[i,j].plot(e_range,np.mean(np.mean(PE_SR_6,axis=2),axis=1)[b_val,:,0],color='r', label='SR')
axes[i,j].plot(e_range,np.mean(np.mean(PE_PR_6,axis=2),axis=1)[b_val,:,0],color='b', label='PR')
axes[i,j].plot(e_range,np.mean(np.mean(PE_WTAR_6,axis=2),axis=1)[b_val,:,0],color='g', label='WTAR')
axes[i,j].set_xlabel('$\epsilon$',fontsize='15')
#axes[i,j].set_ylabel('Per. of agents inf.')
axes[i,j].set_ylabel('${\%}$',fontsize='15')
axes[i,j].text(0.03+min(axes[i,j].get_xlim()),0.9*max(axes[i,j].get_ylim()), ABC[3],fontsize=20)
axes[i,j].legend()
axes[i,j].set_xlim(right=1.25)

i,j=1,1
axes[i,j].plot(e_range,np.mean(np.mean(PE_SR_7,axis=2),axis=1)[b_val,:,0],color='r', label='SR')
axes[i,j].plot(e_range,np.mean(np.mean(PE_PR_7,axis=2),axis=1)[b_val,:,0],color='b', label='PR')
axes[i,j].plot(e_range,np.mean(np.mean(PE_WTAR_7,axis=2),axis=1)[b_val,:,0],color='g', label='WTAR')
axes[i,j].set_xlabel('$\epsilon$',fontsize='15')
#axes[i,j].set_ylabel('Per. of agents inf.')
axes[i,j].set_ylabel('${\%}$',fontsize='15')
axes[i,j].text(0.03+min(axes[i,j].get_xlim()),0.9*max(axes[i,j].get_ylim()), ABC[4],fontsize=20)
axes[i,j].legend()
axes[i,j].set_xlim(right=1.25)

#plt.savefig('/home/glory/Dropbox/Opinion Dynamics/paper/Figures/Fig_mult.pdf',bbox_inches='tight')    


# In[ ]:





# In[ ]:





# In[ ]:





# In[ ]:





# In[22]:


def gaussian(x, sigma, mu):
    return 1/(sigma * np.sqrt(2 * np.pi)) * np.exp(-0.5*((x-mu)/sigma)**2)


def main():
    x = np.linspace(-1, 1, 10000)
    mu1 = -0.4
    mu2 = 0.5
    sigma = 0.15
    ratio = 0.8
    linestyle = "-"
    color = "C1"
    y1 = ratio * gaussian(mu1 + sigma, sigma, mu1)
    y2 = gaussian(mu2 + sigma, sigma, mu2)
    g1 = ratio * gaussian(x, sigma, mu1)
    g2 = gaussian(x, sigma, mu2)
    m1 = np.max(g1)
    m2 = np.max(g2)
    g = g1 + g2
    fig, ax = plt.subplots()
    ax.plot(x, g)
    ymin, ymax = ax.get_ylim()
    ax.plot([mu1, mu1], [ymin, m1], ls=linestyle, c=color)
    ax.plot(mu1, m1, c=color, marker=6)
    ax.plot([mu1, mu1+sigma], [y1, y1], ls=linestyle, c=color)
    ax.plot(mu1, y1, c=color, marker=4)
    ax.plot(mu1+sigma, y1, c=color, marker=5)
    ax.plot([mu1, mu2], [m1, m1], ls=linestyle, c=color)
    ax.plot(mu1, m1, c=color, marker=4)
    ax.plot(mu2, m1, c=color, marker=5)
    ax.axvline((mu2+mu1)/2, ls=linestyle, c=color)
    ax.plot([mu2, mu2], [ymin, m2], ls=linestyle, c=color)
    ax.plot(mu2, m2, c=color, marker=6)
    ax.plot([mu2, mu2+sigma], [y2, y2], ls=linestyle, c=color)
    ax.plot(mu2, y2, c=color, marker=4)
    ax.plot(mu2+sigma, y2, c=color, marker=5)
    ax.text(mu1 + 0.5 * sigma,  y1, r"$\sigma$",ha="center", c=color, va="bottom")

    ax.text(-0.05, m1, r"$\Delta$", ha="center", c=color, va="bottom")
    ax.text(mu1 + 0.01, 0.3 * m1, "p", c=color)
    ax.text(mu2 + 0.01, 0.3 * m1, "(1-p)", c=color)
    ax.text(mu2 + 0.5 * sigma,  y2, r"$\sigma$", ha="center", c=color, va="bottom")

    ax.text((mu2+mu1)/2-0.01, 0.5*m1, r"$\mu$", c=color, ha="right")
    ax.set_ylim(ymin, ymax)
    fig.tight_layout()
    plt.show()

main()


# In[ ]:





# In[97]:


import numpy as np
import matplotlib.pyplot as plt

data=np.load('PAPER_DATA/COMP_MM_MP_SING/OD1/Comp_proportion.npz')
E_MM=data['PE_MM']
E_MP=data['PE_MP']
plt.plot(E_MM[0,0:5,3],E_MP[0,0:5,3],'o')
plt.plot(E_MM[1,0:5,3],E_MP[1,0:5,3],'o')
plt.plot(E_MM[2,0:5,3],E_MP[2,0:5,3],'o')
plt.plot(E_MM[0,0:5,10],E_MP[0,0:5,10],'*')
plt.plot(E_MM[1,0:5,10],E_MP[1,0:5,10],'*')
plt.plot(E_MM[2,0:5,10],E_MP[2,0:5,10],'*')
data=np.load('PAPER_DATA/COMP_MM_MP_SING/OD1/Comp_polarization_0.0.npz')
E_MM=data['PE_MM']
E_MP=data['PE_MP']
plt.plot(E_MM[0,0:5,3],E_MP[0,0:5,3],'o')
plt.plot(E_MM[1,0:5,3],E_MP[1,0:5,3],'o')
plt.plot(E_MM[2,0:5,3],E_MP[2,0:5,3],'o')
plt.plot(E_MM[0,0:5,10],E_MP[0,0:5,10],'*')
plt.plot(E_MM[1,0:5,10],E_MP[1,0:5,10],'*')
plt.plot(E_MM[2,0:5,10],E_MP[2,0:5,10],'*')
data=np.load('PAPER_DATA/COMP_MM_MP_SING/OD1/Comp_polarization_0.05.npz')
E_MM=data['PE_MM']
E_MP=data['PE_MP']
plt.plot(E_MM[0,0:5,3],E_MP[0,0:5,3],'o')
plt.plot(E_MM[1,0:5,3],E_MP[1,0:5,3],'o')
plt.plot(E_MM[2,0:5,3],E_MP[2,0:5,3],'o')
plt.plot(E_MM[0,0:5,10],E_MP[0,0:5,10],'*')
plt.plot(E_MM[1,0:5,10],E_MP[1,0:5,10],'*')
plt.plot(E_MM[2,0:5,10],E_MP[2,0:5,10],'*')
x=np.linspace(0,10,10)
plt.plot(x,x,'--', color='k')
plt.xlabel('Minimum Majority')
plt.ylabel('Minimum People')

#plt.savefig('/home/glory/Dropbox/Opinion Dynamics/paper/Figures/figs3.pdf', bbox_inches='tight')



fig,axes=plt.subplots(1,3,figsize=(20,5))

i,j = 0,0
axes[j].plot(e_range,N_AG_EXT_6[b_val]*(NE_EXT_6[b_val]/Num_simulation),linestyle='--',label='${\%}_{P_1} +{\%}_{P_6}$')
axes[j].plot(e_range,N_AG_MID_6[b_val]*(NE_MID_6[b_val]/Num_simulation),linestyle='--',label='${\%}_{P_2} +{\%}_{P_5}$')
axes[j].plot(e_range,N_AG_CENT_6[b_val]*(NE_CENT_6[b_val]/Num_simulation),linestyle='--',label='${\%}_{P_3} +{\%}_{P_4}$')
#axes[j].plot(e_range,100*(np.mean(NUM_AGENTS_INFLUENCED_6[b_val],axis=0)/Total_agents)[:,0],color='k',linewidth=2,label='$N^{\%}_{P_1} +N^{\%}_{P_2} + N^{\%}_{P_3} +N^{\%}_{P_4} + N^{\%}_{P_5} +\\xi_{P_6}$')
axes[j].plot(e_range,100*(np.mean(NUM_AGENTS_INFLUENCED_6[b_val],axis=0)/Total_agents)[:,0],color='k',linewidth=2,label='$ \sum_{i=1}^p {\%}_{P_i}$')
axes[j].axvline(e_range[np.argmax(100*(np.mean(NUM_AGENTS_INFLUENCED_6[b_val],axis=0)/Total_agents)[:,0])],color='k')
axes[j].axvline(eps_critical_6[b_val],color='b')
axes[j].axvline(e_range[np.argmax(N_AG_EXT_6[b_val]*(NE_EXT_6[b_val]/Num_simulation))],color='r',linewidth=2,linestyle='--')
axes[j].text(e_range[np.argmax(100*(np.mean(NUM_AGENTS_INFLUENCED_6[b_val],axis=0)/Total_agents)[:,0])], np.max(100*(np.mean(NUM_AGENTS_INFLUENCED_6[b_val],axis=0)/Total_agents)[:,0]),'$\epsilon_m$',color='k',fontsize=15)
axes[j].text(e_range[np.argmax(N_AG_EXT_6[b_val]*(NE_EXT_6[b_val]/Num_simulation))],0.5*np.max(N_AG_EXT_6[b_val]*(NE_EXT_6[b_val]/Num_simulation)),'$\epsilon_{ext}$', color='r', fontsize=15)
axes[j].text(eps_critical_6[b_val], 0.75*np.max(100*(np.mean(NUM_AGENTS_INFLUENCED_6[b_val],axis=0)/Total_agents)[:,0]), '$\epsilon_c$', color='b',fontsize=15)
axes[j].text(0.1+min(axes[j].get_xlim()),0.9*max(axes[j].get_ylim()), ABC[3],fontsize=20)
#axes[j].axvline(e_range[np.argmax(E_EXT_6[b_val]*(NE_EXT_6[b_val]/Num_simulation))],color='g',linestyle='--', label='$\epsilon_{ext}(extr.)$')
axes[j].set_xlabel('$\epsilon$',fontsize='15')    
axes[j].set_ylabel('${\%}_{P_i}$',fontsize='15')
axes[j].legend()
axes[j].set_xlim([0.01,1.25])

i,j = 0,1
i = j
Num_party=6
bottom=0
axes[i].bar(e_range,PP6[b_val,:,0],width=.01,label='P '+str(1),color=color(0))
for num in range(Num_party-1):
    bottom=bottom+PP6[b_val,:,num]
    axes[i].bar(e_range,PP6[b_val,:,num+1],width=.01,bottom=bottom,label='P '+str(num+2),color=color(num+1))
axes[i].text(0.1+min(axes[i].get_xlim()),0.9*max(axes[i].get_ylim()), ABC[3],fontsize=20)    
#axes[i].set_xlabel('$\epsilon$',fontsize='15')
#axes[i].set_ylabel('% of win')
axes[i].legend(bbox_to_anchor =(0.97, 0.98))
axes[i].set_xlim([0,1.25])
#axes[i].set_title('p = 6')

i,j = 0,2
i = j
Num_party=6
bottom=0
axes[i].bar(e_range,PP6_2[b_val,:,0],width=.01,label='P '+str(1),color=color(0))
for num in range(Num_party-1):
    bottom=bottom+PP6_2[b_val,:,num]
    axes[i].bar(e_range,PP6_2[b_val,:,num+1],width=.01,bottom=bottom,label='P '+str(num+2),color=color(num+1))
axes[i].text(0.1+min(axes[i].get_xlim()),0.9*max(axes[i].get_ylim()), ABC[3],fontsize=20)    
#axes[i].set_xlabel('$\epsilon$',fontsize='15')
#axes[i].set_ylabel('% of win')
axes[i].legend(bbox_to_anchor =(0.97, 0.98))
axes[i].set_xlim([0,1.25])
#axes[i].set_title('p = 6')



