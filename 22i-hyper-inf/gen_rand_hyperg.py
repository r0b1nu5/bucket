import xgi

run = '013'

n = 300
ps = [.001,.00001]

H = xgi.random_hypergraph(n,ps)
#while not xgi.is_connected(H):
#    H = xgi.random_hypergraph(n,ps)

xgi.write_edgelist(H,"data/edgelist-n"+str(n)+"-"+run+".csv",delimiter=",")


