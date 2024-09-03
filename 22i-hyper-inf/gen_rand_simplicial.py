import xgi

run = '012'

n = 500
ps = [.001,.0002]

H = xgi.random_simplicial_complex(n,ps)
while not xgi.is_connected(H):
    H = xgi.random_simplicial_complex(n,ps)

xgi.write_edgelist(H,"data/edgelist-n"+str(n)+"-"+run+".csv",delimiter=",")


