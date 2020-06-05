include("big_rand.jl")

NS = 10			# Number of states
N = rand(200:1000,NS)	# Population in each state
M = round.(Int64,N)	# Number of representative for each state
MU = .5*ones(NS)	# Mean opinions in each state
SI = .2*ones(NS)	# Opinions' standard deviation in each state
DE = 0*ones(NS)		# Opinions' bias in each state

eps = .2

X0 = Array{Array{Float64,1},1}()	# Natural opinions in each state
for i = 1:NS
	push!(X0,big_rand(N[i],-MU[i]+DE[i],SI[i],MU[i]+DE[i],SI[i]))
end

Y0 = Array{Array{Float64,1},1}()	# Final opinions before influence

for i in 1:NS






