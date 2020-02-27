using PyPlot, DelimitedFiles, SparseArrays

include("kuramoto.jl")

wmin = .0001
wmax = .1
nw = 30
ws = exp.(LinRange(log(1/wmax),log(1/wmin),nw))

h = .1

a0 = .2

max_iter = 100000
n_store = ceil(Int,2*pi/wmin/h)
max_iter = max(100000,n_store)

Lsp = readdlm("ntws_data/uk_w_lap_mat_sp.csv",',')
L = sparse(Lsp[:,1],Lsp[:,2],Lsp[:,3])

n = size(L)[1]

P = zeros(n)
th0 = rand(n)

j = rand(1:n)

a = zeros(n)
a[j] = a0

nhs = Array{Float64,1}()

c = 0

for w0 in ws
	global c

	c += 1
	@info "$c/$nw"
	
	w = zeros(n)
	w[j] = w0

	ths = kuramoto_sine(L,P,th0,a,w,zeros(n),n_store,max_iter,1e-8,.1,"data1")
	
	A = maximum(ths,dims=2) - minimum(ths,dims=2)

	nh = mean(2*a0./(A*w0))

	push!(nhs,nh)
end

PyPlot.plot(1 ./ws,nhs,"o")




