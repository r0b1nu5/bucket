import pandapower as pp
import pandas as pd
import numpy as np
import os
import csv

net = pp.from_json("data-marc/network.json")

mup = pd.read_csv("data-marc/mup.csv",header=None)
sigp = pd.read_csv("data-marc/sigp.csv",header=None)
muq = pd.read_csv("data-marc/muq.csv",header=None)
sigq = pd.read_csv("data-marc/sigq.csv",header=None)

times = pd.read_csv("data-marc/times.csv",header=None)

head = ["Time",]
head.extend(range(55))

with open("data-robin/P.csv", 'w', newline='') as csvfile:
    writer = csv.writer(csvfile)
    writer.writerow(head)
with open("data-robin/Q.csv", 'w', newline='') as csvfile:
    writer = csv.writer(csvfile)
    writer.writerow(head)
with open("data-robin/V.csv", 'w', newline='') as csvfile:
    writer = csv.writer(csvfile)
    writer.writerow(head)

rho = 1.1

for i in range(86400):
    if np.mod(i,100) == 0:
        print(i)

    time = times.loc[i].values[0]
    ps = [time,]
    qs = [time,]
    vs = [time,]

    for i in range(len(mup)):
        p = np.random.normal(loc=mup.loc[i].values[0], scale=sigp.loc[i].values[0]*rho)
        q = np.random.normal(loc=muq.loc[i].values[0], scale=sigq.loc[i].values[0]*rho)

        net.load.p_mw.at[i] = 1e-6*p
        ps.append(p)
        net.load.q_mvar.at[i] = 1e-6*q
        qs.append(q)

    pp.runpp(net)

    for l in net.load.bus:
        vs.append(net.res_bus.vm_pu.loc[l])

    with open(f"data-robin/P_{rho}.csv", 'a', newline='') as csvfile:
        writer = csv.writer(csvfile)
        writer.writerow(ps)
    with open(f"data-robin/Q_{rho}.csv", 'a', newline='') as csvfile:
        writer = csv.writer(csvfile)
        writer.writerow(qs)
    with open(f"data-robin/V_{rho}.csv", 'a', newline='') as csvfile:
        writer = csv.writer(csvfile)
        writer.writerow(vs)



