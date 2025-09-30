import pandapower as pp
import pandas as pd
import numpy as np
import os
import csv

net = pp.from_json("rv/data-marc/network.json")

mup = pd.read_csv("rv/data-marc/mup.csv",header=None)
sigp = pd.read_csv("rv/data-marc/sigp.csv",header=None)
muq = pd.read_csv("rv/data-marc/muq.csv",header=None)
sigq = pd.read_csv("rv/data-marc/sigq.csv",header=None)

times = pd.read_csv("rv/data-marc/times.csv",header=None)

head = ["Time",]
head.extend(range(55))

rs = np.linspace(0,0.5,6)

for r in rs:
    with open(f"data-robin/P_{r}.csv", 'w', newline='') as csvfile:
        writer = csv.writer(csvfile)
        writer.writerow(head)
    with open(f"data-robin/Q_{r}.csv", 'w', newline='') as csvfile:
        writer = csv.writer(csvfile)
        writer.writerow(head)
    with open(f"data-robin/V_{r}.csv", 'w', newline='') as csvfile:
        writer = csv.writer(csvfile)
        writer.writerow(head)

    for i in range(10000):
        if np.mod(i,100) == 0:
            print(i)

        time = times.loc[i].values[0]
        ps = [time,]
        qs = [time,]
        vs = [time,]

        x0 = np.random.normal()
        for i in range(len(mup)):
            p = mup.loc[i].values[0] + sigp.loc[i].values[0]*((1-r)*np.random.normal() + r*x0)
            q = muq.loc[i].values[0] + sigq.loc[i].values[0]*((1-r)*np.random.normal() + r*x0)

            net.load.p_mw.at[i] = 1e-6*p
            ps.append(p)
            net.load.q_mvar.at[i] = 1e-6*q
            qs.append(q)

        pp.runpp(net)

        for l in net.load.bus:
            vs.append(net.res_bus.vm_pu.loc[l])

        with open(f"data-robin/P_{r}.csv", 'a', newline='') as csvfile:
            writer = csv.writer(csvfile)
            writer.writerow(ps)
        with open(f"data-robin/Q_{r}.csv", 'a', newline='') as csvfile:
            writer = csv.writer(csvfile)
            writer.writerow(qs)
        with open(f"data-robin/V_{r}.csv", 'a', newline='') as csvfile:
            writer = csv.writer(csvfile)
            writer.writerow(vs)



