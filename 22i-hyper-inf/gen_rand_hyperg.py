import xgi

run = '007'

n = 30
ps = [.05,.005]

H = xgi.random_hypergraph(n,ps)
while not xgi.is_connected(H):
    H = xgi.random_hypergraph(n,ps)

xgi.write_edgelist(H,"data/edgelist-n"+str(n)+"-"+run+".csv",delimiter=",")


