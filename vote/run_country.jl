using LinearAlgebra

include("big_rand.jl")
include("scripts_mini.jl")
include("result.jl")

 #=
NS = 9				# Number of states
N = rand(200:1000,NS)		# Population in each state
M = 2*round.(Int64,N/100) .+ 1	# Number of representative for each state
MU = .5*ones(NS)		# Mean opinions in each state
SI = .2*ones(NS)		# Opinions' standard deviation in each state
DE = .1*ones(NS)		# Opinions' bias in each state

w_ref = 0.1
emi = .3
eps = emi + .01
# =#
global X0 = Array{Array{Float64,1},1}()	# Natural opinions in each state
for i in 1:NS
	global emi
	maxd = 1000.
	while maxd > emi
		global x = big_rand(N[i],-MU[i]+DE[i],SI[i],MU[i]+DE[i],SI[i])
		maxd = maximum(x[2:end]-x[1:end-1])
	end
	push!(X0,x)
end

X1 = Array{Array{Float64,1},1}()
state_id = Array{Int64,1}()
LDis = Array{Array{Float64,2},1}()
for i in 1:NS
	global state_id

	A = Float64.((0 .< abs.(repeat(X0[i],1,N[i]) - repeat(X0[i]',N[i],1)) .< eps))
	d = vec(sum(A,dims=1))
	D = diagm(0 => d)
	L = D - A
	LpD = Symmetric(2*D - A)
#	LpD = 2*D - A
	LDi = inv(LpD)*Diagonal(D)

	push!(X1,LDi*X0[i])
	push!(LDis,LDi)
	state_id = [state_id;i*ones(Int64,N[i])]
end

# ======================================================================================
# Single representative	
# ======================================================================================
# #=
x0 = Array{Float64,1}()
x1 = Array{Float64,1}()
for i in 1:NS
	global x0 = [x0;X0[i]]
	global x1 = [x1;X1[i]]
end
oc0 = outcome(x1)[1]
oc1 = copy(oc0)
w0 = -sign(oc0)*w_ref
inf_order = mini_sort(x0,(w0 < 0))
c1 = 0

while oc0*oc1 >= 0.
	c2 = 0
	while oc0*oc1 >= 0. && c2 < length(x1)
		global c1 += 1
		c2 += 1
		to_inf = inf_order[c2]
		state = state_id[to_inf]
		index = setdiff((1:length(x1)).*(state_id .== state),[0.])
		x1[index] += w0*LDis[state][:,to_inf-sum(N[1:state-1])]
		global oc1 = outcome(x1)[1]
	end
end

xi1 = c1*abs(w0)
# =#

# ======================================================================================
# Several representatives
# ======================================================================================
# #=
X2 = copy(X1)
Mp0 = Array{Int64,1}()
margMp0 = Array{Float64,1}()
Mn0 = Array{Int64,1}()
margMn0 = Array{Float64,1}()
for i in 1:NS
	xxx = result_margin(X1[i],M[i])
	push!(Mp0,xxx[1])
	push!(margMp0,xxx[2])
	push!(Mn0,xxx[3])
	push!(margMn0,xxx[4])
end

oc0 = sum(Mp0) - sum(Mn0)
oc1 = copy(oc0)
Mp1 = copy(Mp0)
Mn1 = copy(Mn0)
margMp1 = copy(margMp0)
margMn1 = copy(margMn0)
w0 = -sign(oc0)*w_ref
c1 = 0

while oc0*oc1 > 0
	if w0 > 0
		state_order = Int64.(sortslices([margMn1 1:NS],dims=1)[:,2])
	else
		state_order = Int64.(sortslices([margMp1 1:NS],dims=1)[:,2])
	end
	s = state_order[1]
	mp1 = Mp1[s]
	mp2 = copy(mp1)
	inf_order = mini_sort(X0[s],(w0 < 0))
	while mp2 == mp1
		c2 = 0
		while mp2 == mp1 && c2 < length(X2[s])
			global c1 += 1
			c2 += 1
			to_inf = inf_order[c2]
			global X2[s] += w0*LDis[s][:,to_inf]
			mp2,mn2 = result(X2[s],M[s])
		end
	end
	global Mp1[s],margMp1[s],Mn1[s],margMn1[s] = result_margin(X2[s],M[s])
	global oc1 = sum(Mp1) - sum(Mn1)
end

xi2 = c1*abs(w0)
# =#


# ======================================================================================
# Winner takes all
# ======================================================================================
# #=
X3 = copy(X1)
Mp0 = Array{Int64,1}()
margMp0 = Array{Float64,1}()
Mn0 = Array{Int64,1}()
margMn0 = Array{Float64,1}()
for i in 1:NS
	xxx = result_margin(X1[i],1)
	push!(Mp0,M[i]*xxx[1])
	push!(margMp0,xxx[2])
	push!(Mn0,M[i]*xxx[3])
	push!(margMn0,xxx[4])
end

oc0 = sum(Mp0) - sum(Mn0)
oc1 = copy(oc0)
Mp1 = copy(Mp0)
Mn1 = copy(Mn0)
margMp1 = copy(margMp0)
margMn1 = copy(margMn0)
w0 = -sign(oc0)*w_ref
c1 = 0

while oc0*oc1 > 0
	if w0 > 0
		state_order = Int64.(sortslices([M./margMn1 1:NS],dims=1,rev=true)[:,2])
	else
		state_order = Int64.(sortslices([M./margMp1 1:NS],dims=1,rev=true)[:,2])
	end
	s = state_order[1]
	mp1 = Mp1[s]
	mp2 = copy(mp1)
	inf_order = mini_sort(X0[s],(w0 < 0))
	while mp2 == mp1
		c2 = 0
		while mp2 == mp1 && c2 < length(X3[s])
			global c1 += 1
			c2 += 1
			to_inf = inf_order[c2]
			global X3[s] += w0*LDis[s][:,to_inf]
			r = result(X3[s],1)
			mp2 = M[s]*r[1]
			mn2 = M[s]*r[2]
		end
	end

	R = result_margin(X3[s],1)
	global Mp1[s] = M[s]*R[1]
	margMp1[s] = R[2]
	global Mn1[s] = M[s]*R[3]
	margMn1[s] = M[s]*R[4]
	
	global oc1 = sum(Mp1) - sum(Mn1)
end

xi3 = c1*abs(w0)
# =#







