using PowerModels, PyPlot, SparseArrays

function load_rate(nd::Dict{String,Any}, do_plot::Bool=true)
	θ = angle_vec(nd)
	L = get_L(nd)
	angbnd = angdiff_bound(nd)

	return load_rate(θ,L,angbnd,do_plot)
end

function load_rate(θ::Array{Float64,1}, L::SparseMatrixCSC{Float64,Int64}, angbnd::SparseMatrixCSC{Float64,Int64}, do_plot::Bool=true)
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
	end

	return sparse(I2,J2,L)
end


function angle_vec(nd::Dict{String,Any})
	θ = Array{Float64,1}()

	i2b = calc_admittance_matrix(nd).idx_to_bus

	for i in 1:length(i2b)
		push!(θ,nd["bus"]["$(i2b[i])"]["va"])
	end

	return θ
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



