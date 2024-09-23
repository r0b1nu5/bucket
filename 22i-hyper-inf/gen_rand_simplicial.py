import xgi

run = '014'

n = 300
ps = [.001,.00001]

H = xgi.random_simplicial_complex(n,ps)
#while not xgi.is_connected(H):
#    H = xgi.random_simplicial_complex(n,ps)

xgi.write_edgelist(H,"data/edgelist-n"+str(n)+"-"+run+".csv",delimiter=",")


