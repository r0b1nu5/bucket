using PowerModels, PyPlot, SparseArrays

include("tools.jl")

# Computes the ratio between actual angle difference and maximal angle difference for each line.
function anglediff_rate(nd::Dict{String,Any}, do_plot::Bool=true)
	θ = angle_vec(nd)
	angbnd = angdiff_bound(nd)

	return anglediff_rate(θ,angbnd,do_plot)
end

function anglediff_rate(θ::Vector{Float64}, nd::Dict{String,Any}, do_plot::Bool=true)
	angbnd = angdiff_bound(nd)

	return anglediff_rate(θ,angbnd,do_plot)
end

function anglediff_rate(θ::Array{Float64,1}, angbnd::SparseMatrixCSC{Float64,Int64}, do_plot::Bool=true)
	n = length(θ)

	I,J,V = findnz(angbnd)

	ids = setdiff((I .< J).*(1:length(I)),[0,])
	I2 = I[ids]
	J2 = J[ids]

	L = Array{Float64,1}()
	for k in 1:length(ids)
		i = I2[k]
		j = J2[k]

		δ = θ[i] - θ[j]
		push!(L,max(δ/angbnd[i,j],δ/angbnd[j,i]))
	end

	if do_plot
		c = (L .> .5) + (L .> .8) + (L .> 1.) .+ 1
		col = ["C2","C0","C1","C3"]

		figure()
		PyPlot.plot([0,length(ids)+1],[1.,1.],"--k")
		for k in 1:length(ids)
			PyPlot.bar(k,L[k],.9,0.,color=col[c[k]])
		end
		
		xticks(Array(1:length(ids)),["($(I[i]),$(J[i]))" for i in ids],rotation=90)
		ylabel("angle ratio")
	end

	return sparse([I2;J2],[J2;I2],[L;L])
end

# Computes the ratio between actual active power flow and the maximal active power flow for each line.
function flow_rate(nd::Dict{String,Any}, to_plot::Bool=true)
	θ = angle_vec(nd)
	v = volt_vec(nd)
	angbnd = angdiff_bound(nd)
	B = calc_basic_susceptance_matrix(nd)

	return flow_rate(θ,v,angbnd,B,to_plot)
end

function flow_rate(θ::Vector{Float64}, v::Vector{Float64}, nd::Dict{String,Any}, do_plot::Bool=true)
	angbnd = angdiff_bound(nd)
	B = calc_basic_susceptance_matrix(nd)

	return flow_rate(θ,v,angbnd,B,to_plot)
end

function flow_rate(θ::Vector{Float64}, v::Vector{Float64}, angbnd::SparseMatrixCSC{Float64,Int64}, B::SparseMatrixCSC{Float64,Int64}, do_plot::Bool=true)
	n = length(θ)

	I,J,V = findnz(angbnd)

	ids = setdiff((I .< J).*(1:length(I)),[0,])
	I2 = I[ids]
	J2 = J[ids]

	L = Vector{Float64}()
	for k in 1:length(ids)
		i = I2[k]
		j = J2[k]

		Bvv = B[i,j]*v[i]*v[j]

		f = Bvv*sin(θ[i] - θ[j])
		fmax = Bvv*sin(angbnd[i,j])
		fmin = Bvv*sin(angbnd[j,i])
		push!(L,max(f/fmax,f/fmin))
	end

	if do_plot
		c = (L .> .5) + (L .> .8) + (L .> 1.) .+ 1
		col = ["C2","C0","C1","C3"]

		figure()
		PyPlot.plot([0,length(ids)+1],[1.,1.],"--k")
		for k in 1:length(ids)
			PyPlot.bar(k,L[k],.9,0.,color=col[c[k]])
		end
		
		xticks(Array(1:length(ids)),["($(I[i]),$(J[i]))" for i in ids],rotation=90)
		ylabel("flow ratio")
	end

	return sparse([I2;J2],[J2;I2],[L;L])
end
	
#Computes the active power flows on each line.
function active_flow(nd::Dict{String,Any}, to_plot::Bool=true)
	θ = angle_vec(nd)
	v = volt_vec(nd)
	B = calc_basic_susceptance_matrix(nd)

	return active_flow(θ,v,B,to_plot)
end

function active_flow(θ::Vector{Float64}, v::Vector{Float64}, nd::Dict{String,Any}, do_plot::Bool=false)
	B = calc_basic_susceptance_matrix(nd)
	
	return active_flow(θ,v,B,do_plot)
end

function active_flow(θ::Vector{Float64}, v::Vector{Float64}, B::SparseMatrixCSC{Float64,Int64}, do_plot::Bool=false)
	n = length(θ)

	I,J,V = findnz(B)

	ids = setdiff((I .< J).*(1:length(I)),[0,])
	I2 = I[ids]
	J2 = J[ids]

	L = Vector{Float64}()
	for k in 1:length(ids)
		i = I2[k]
		j = J2[k]

		push!(L,B[i,j]*v[i]*v[j]*sin(θ[i]-θ[j]))
	end

	if do_plot
		figure()
		PyPlot.plot([0,length(ids)+1],[1.,1.],"--k")
		for k in 1:length(ids)
			PyPlot.bar(k,L[k],.9,0.,color="C0")
		end
		
		xticks(Array(1:length(ids)),["($(I[i]),$(J[i]))" for i in ids],rotation=90)
		ylabel("active power flow")
	end

	return sparse([I2;J2],[J2;I2],[L;-L])
end
	

function dist_ratio(Lw::SparseMatrixCSC{Float64,Int64}, do_plot::Bool=false)
	Ω = res_dist(Array(Lw))
	D = geo_dist(Lw)

	K = Ω./D

	if do_plot
		n = size(Lw)[1]
		figure()
		for i in 1:n
			k = sort(K[:,i])
			subplot(1,2,1)
			PyPlot.semilogy(1:n,k)
			subplot(1,2,2)
			PyPlot.semilogy(1:n-1,k[2:end]-k[1:end-1])
		end
	end

	return K
end

function dist_ratio(I::Array{Int64,1}, Lw::SparseMatrixCSC{Float64,Int64}, do_plot::Bool=false)
	n = size(Lw)[1]

	Ω = res_dist(Array(Lw))
	ΩI = Ω[:,I]

	DI = Array{Float64,2}(undef,n,0)
	for i in I
		DI = [DI geo_dist(i,Lw)]
	end

	KI = ΩI./DI

	if do_plot
		n = size(Lw)[1]
		figure()
		for i in length(I)
			k = sort(KI[:,i])
			subplot(1,2,1)
			PyPlot.semilogy(1:n,k,label="i = $(I[i])")
			subplot(1,2,2)
			PyPlot.semilogy(1:n-1,k[2:end]-k[1:end-1])
		end
	end

	return KI
end

function dist_ratio(i::Int64, Lw::SparseMatrixCSC{Float64,Int64}, do_plot::Bool=false)
	Ω = res_dist(Array(Lw))
	Di = geo_dist(i,Lw)

	k = Ω[:,i]./Di

	if do_plot
		kk = sort(k)
		figure()
		subplot(1,2,1)
		PyPlot.plot(1:length(k),kk)
		subplot(1,2,2)
		PyPlot.plot(1:length(k)-1,kk[2:end]-kk[1:end-1])
	end

	return k
end


function angle_vec(nd::Dict{String,Any})
	θ = Array{Float64,1}()

	i2b = calc_admittance_matrix(nd).idx_to_bus

	for i in 1:length(i2b)
		push!(θ,nd["bus"]["$(i2b[i])"]["va"])
	end

	return θ
end

function volt_vec(nd::Dict{String,Any})
	v = Array{Float64,1}()

	i2b = calc_admittance_matrix(nd).idx_to_bus

	for i in 1:length(i2b)
		push!(v,nd["bus"]["$(i2b[i])"]["vm"])
	end

	return v
end


function angdiff_bound(nd::Dict{String,Any})
	n = length(nd["bus"])

	I = Array{Int64,1}()
	J = Array{Int64,1}()
	V = Array{Float64,1}()

	b2i = calc_admittance_matrix(nd).bus_to_idx

	for (l,branch) in nd["branch"]
		f = b2i[branch["f_bus"]]
		t = b2i[branch["t_bus"]]
		if f < t
			push!(I,f)
			push!(J,t)
			push!(V,branch["angmax"])
			push!(I,t)
			push!(J,f)
			push!(V,branch["angmin"])
		else
			push!(I,f)
			push!(J,t)
			push!(V,-branch["angmin"])
			push!(I,t)
			push!(J,f)
			push!(V,-branch["angmax"])
		end
	end

	return sparse(I,J,V)
end

function get_L(nd::Dict{String,Any})
	n = length(nd["bus"])

	I,J,V = findnz(calc_admittance_matrix(nd).matrix)

	ids = setdiff((I .!= J).*(1:length(I)),[0,])
	I2 = I[ids]
	J2 = J[ids]
	V2 = ones(length(ids))

	I3 = Array(1:n)
	J3 = Array(1:n)
	V3 = vec(sum(sparse(I2,J2,V2),dims=2))

	return sparse([I2;I3],[J2;J3],[-V2;V3])
end

function get_weighted_L(nd::Dict{String,Any})
	n = length(nd["bus"])

	Lw = calc_admittance_matrix(nd).matrix
	Aw = spdiagm(0 => diag(Lw)) - Lw
	B = imag.(Aw)
	
	θ = angle_vec(nd)
	v = volt_vec(nd)

	return get_weighted_L(B,θ,v)
end


function get_weighted_L(BB::SparseMatrixCSC{Float64,Int64}, θ::Array{Float64,1}, v::Array{Float64,1})
	n = size(BB)[1]

	BB = abs.(BB - spdiagm(0 => diag(BB)))
	I,J,B = findnz(BB)

	Ri = spzeros(n,n)

	for k in 1:length(I)
		i = I[k]
		j = J[k]

		Ri[i,j] = B[k]*v[i]*v[j]*cos(θ[i] - θ[j])
	end

	D = spdiagm(0 => vec(sum(Ri,dims=1)))

	return D - Ri
end





