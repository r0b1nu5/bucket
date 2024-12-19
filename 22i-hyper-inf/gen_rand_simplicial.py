import xgi

run = '998'

#n = 4; ps = [.3,.2]
n = 8; ps = [.1,.05]
#n = 16; ps = [.05,.01]
#n = 32; ps = [.01,.005]
#n = 64; ps = [.005,.001]
#n = 128; ps = [.002,.0002]
#n = 256; ps = [.001,.00001]

H = xgi.random_simplicial_complex(n,ps)
#while not xgi.is_connected(H):
#    H = xgi.random_simplicial_complex(n,ps)

xgi.write_edgelist(H,"data/edgelist-n"+str(n)+"-"+run+".csv",delimiter=",")


