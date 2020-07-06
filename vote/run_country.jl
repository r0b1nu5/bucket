include("big_rand.jl")
include("scripts_mini.jl")

NS = 10			# Number of states
N = rand(200:1000,NS)	# Population in each state
M = round.(Int64,N/50)	# Number of representative for each state
MU = .5*ones(NS)	# Mean opinions in each state
SI = .2*ones(NS)	# Opinions' standard deviation in each state
DE = 0*ones(NS)		# Opinions' bias in each state

eps = .2

X0 = Array{Array{Float64,1},1}()	# Natural opinions in each state
for i = 1:NS
	push!(X0,big_rand(N[i],-MU[i]+DE[i],SI[i],MU[i]+DE[i],SI[i]))
end

Y0 = Array{Array{Float64,1},1}()	# Final opinions before influence

Xs = Array{Array{Float64,1},1}()
state_id = Array{Int64,1}()
LDis = Array{Array{Float64,2},1}()
for i in 1:NS
	A = Float.((0 .< abs.(repeat(X0[i],1,N[i]) - repeat(X0[i]',N[i],1)) .< eps))
	d = vec(sum(A,dims=1))
	D = diagm(0 => d)
	L = D - A
	LpD = Symmetric(2*D - A)
	LDi = inv(LpD)*Diagonal(D)

	push!(Xs,LDi*X0[i])
	push!(LDis,LDi)
	state_id = [state_id;i*ones(Int64,N[i])]
end

# Single representative	

xs = vec(Xs)
oc0 = outcome(xs)
oc1 = copy(oc0)
w0 = -sign(oc0)*0.1
inf_order = mini_sort(xs,true)
c = 0

while oc0*oc1 >= 0.
	c += 1
	to_inf = inf_order[c]
	state = state_id[to_inf]
	index = (1:length(xs)).*(state_id .== state)
	xs[index] += w0*LDis[state][:,to_inf-sum(N[i] for i in 1:state-1)]
	oc1 = outcome(xs)
end

xi1 = c*0.1


# Several representatives









