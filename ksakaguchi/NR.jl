using LinearAlgebra, PowerModels

# Newton-Raphson solver for the power flow equations
#
## INPUT
# Y: admittance matrix
# v0: inital voltages
# θ0: initial thetas
# P0: initial active powers
# Q0: initial reactive powers
# bus_type: type of the buses (1: PQ, 2: PV, 3: slack)
#
## OUTPUT
# V: updated voltages
# T: updated thetas
# n_iter: # of iterations before convergence
function NR_solver(Y::Matrix{ComplexF64}, v0::Vector{Float64}, θ0::Vector{Float64}, P0::Vector{Float64}, Q0::Vector{Float64}, bus_type::Vector{Int64}, ϵ::Float64=1e-10, max_iter::Int64=100)
	n = size(Y)[1]
	
	# node ids whose bus type is 0
	PQ_ids = filter(i -> bus_type[i] == 1,1:n)

	#GS_solver(sp)

	# compute Y element-wise absolute values
    	Y_abs = abs.(Y)
	# compute Y element-wise angle values
    	Y_angle = angle.(Y)
	
	# node id whose bus type is 3
	slack_id = filter(i -> bus_type[i] == 3, 1:n)
	# ids from 1:n except slack_id
	ids = setdiff(1:n,slack_id)
    	
	error  = ϵ
    	n_iter = 1

	V = v0
	T = θ0

    	while(error >= ϵ && n_iter < max_iter)
		# pre-computations 1
		M1 = V*V'.*Y_abs 
		M2 = repeat(T,1,n)-repeat(T',n,1)-Y_angle
		M3 = M1.*sin.(M2)
		M4 = M1.*cos.(M2)
		V1 = diag(Y_abs).*sin.(diag(Y_angle)).*V.^2
		V2 = diag(Y_abs).*cos.(diag(Y_angle)).*V.^2

		# compute P and Q (n-dimensional vectors)
		P = vec(sum(M4,dims=2))
		Q = vec(sum(M3,dims=2))

		# pre-computations 2
		V3 = -Q-V1
		V4 = Q-V1
		V5 = P+V2
		V6 = P-V2
		
		# compare P and Q with their "known" value to get dPQ ((n-1+m)-dimensional vector)
		dP = P0[ids] - P[ids]
		dQ = Q0[PQ_ids] - Q[PQ_ids]
		dPQ = [dP;dQ]
		
		# compute the jacobian matrix
		# see Power System Analysis, Bergen/Vital, p345-347
		# L: (n-1)x(n-1), O: mxm, N: (n-1)xm M: mx(n-1), J:(n-1+m)x(n-1xm)
		L = M3[ids,ids] - diagm(0 => diag(M3[ids,ids])) + diagm(0 => V3[ids])
		O = M3[PQ_ids,PQ_ids] - diagm(0 => diag(M3[PQ_ids,PQ_ids])) + diagm(0 => V4[PQ_ids])
		N = M4 - diagm(0 => diag(M4)) + diagm(0 => V5)
		N = N[ids,PQ_ids]
		M = -M4 + diagm(0 => diag(M4)) + diagm(0 => V6)
		M = M[PQ_ids,ids]
		J = [L N;M O]
		# solve JX = dPQ
		X = J\dPQ
		

		# update V and theta
		T[ids] += X[1:n-1]
		V[PQ_ids] += X[n:end]
		error = norm(dPQ,Inf)
		n_iter += 1
    	end

	return V,T,n_iter
end

function NR_solver(nd::Dict{String,Any})
	Y,b2i,i2b = get_Y(nd)

	v0,θ0 = get_voltages(nd,i2b)

	P0,Q0 = get_powers(nd,b2i)

	bt = get_bus_types(nd,i2b)

	return NR_solver(Y,v0,θ0,P0,Q0,bt)
end

function get_Y(nd::Dict{String,Any})
	Ys = calc_admittance_matrix(nd)

	return Matrix(Ys.matrix),Ys.bus_to_idx,Ys.idx_to_bus
end

function get_voltages(nd::Dict{String,Any}, i2b::Vector{Int64})
	v0 = [nd["bus"]["$b"]["vm"] for b in i2b]
	θ0 = [nd["bus"]["$b"]["va"] for b in i2b]

	return v0,θ0
end

function get_powers(nd::Dict{String,Any}, b2i::Dict{Int64,Int64})
	P0 = zeros(n)
	Q0 = zeros(n)
	for k in keys(nd["gen"])
		P0[b2i[nd["gen"][k]["gen_bus"]]] += nd["gen"][k]["pg"]
		Q0[b2i[nd["gen"][k]["gen_bus"]]] += nd["gen"][k]["qg"]
	end
	for k in keys(nd["load"])
		P0[b2i[nd["load"][k]["load_bus"]]] -= nd["load"][k]["pd"]
		Q0[b2i[nd["load"][k]["load_bus"]]] -= nd["load"][k]["qd"]
	end
	
	return P0,Q0
end

function get_bus_types(nd::Dict{String,Any}, i2b::Vector{Int64})
	return [nd["bus"]["$i"]["bus_type"] for i in i2b]
end


