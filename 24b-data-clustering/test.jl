include("kmeans.jl")

n = 10000
d = 2
k = 5

#X = randn(n,d)
X = randn(n,d) + 2*rand([-1,1],n,d)
c0 = randn(k,d)
c = c0

for i in 1:10
	global g,c = my_kmeans(X,c,[1.,2.],1)
	clf()
	plot_grps(X,g)
	pause(1)
end




