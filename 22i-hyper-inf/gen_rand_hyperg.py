import xgi

run = '010'

n = 1000
ps = [.0005,.000001]

H = xgi.random_hypergraph(n,ps)
while not xgi.is_connected(H):
    H = xgi.random_hypergraph(n,ps)

xgi.write_edgelist(H,"data/edgelist-n"+str(n)+"-"+run+".csv",delimiter=",")


