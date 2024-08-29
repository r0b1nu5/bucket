import xgi

run = '002'

n = 100
ps = [.01,.001]

H = xgi.random_simplicial_complex(n,ps)
while not xgi.is_connected(H):
    H = xgi.random_simplicial_complex(n,ps)

xgi.write_edgelist(H,"data/edgelist-n"+str(n)+"-"+run+".csv",delimiter=",")


