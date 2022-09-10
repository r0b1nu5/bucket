using PowerModels, Statistics

include("NR.jl")

#ntw_data = parse_file("data/IEEE_118_Bus.raw"; import_all=true)


function run_Nm1_PowerModels(file::String)
	ntw_data = parse_file(file*".raw"; import_all=true)
	exep = readdlm(file*"_excep.csv",',')
	exe = ["$(Int(exep[i]))" for i in 1:length(exep)]

	sol_ref = compute_ac_pf(ntw_data)
	update_data!(ntw_data,sol_ref)
	flows_ref = calc_branch_flow_ac(ntw_data)
	update_data!(ntw_data,flows_ref)

	line_ref = Dict{String,Any}()
	for kl in keys(ntw_data["branch"])
		dθ = abs(ntw_data["bus"]["$(ntw_data["branch"][kl]["f_bus"])"]["va"] - ntw_data["bus"]["$(ntw_data["branch"][kl]["t_bus"])"]["va"])
		fl = abs(ntw_data["branch"][kl]["pt"]-ntw_data["branch"][kl]["pf"])/2

		line_ref[kl] = Dict{String,Any}("dθ" => dθ, "fl" => fl)
	end


	impact = Dict{String,Any}()
	
	for kl in setdiff(keys(ntw_data["branch"]),exe)
		nd = parse_file(file; import_all=true)
		delete!(nd["branch"],kl)

@info kl
		sol_k = compute_ac_pf(nd)

		if length(sol_k) == 1
			impact[kl] = Dict{String,Any}("converged" => false)
		else
			update_data!(nd,sol_k)
			flows_k = calc_branch_flow_ac(nd)
			update_data!(nd,flows_k)

			y1 = Array{Float64,1}()
			y2 = Array{Float64,1}()
			y3 = Array{Float64,1}()
			y4 = Array{Float64,1}()

			for l in keys(nd["branch"])
				dθ = abs(nd["bus"]["$(nd["branch"][l]["f_bus"])"]["va"] - nd["bus"]["$(nd["branch"][l]["t_bus"])"]["va"])
				fl = abs(nd["branch"][l]["pt"]-nd["branch"][l]["pf"])/2

				push!(y1,dθ-line_ref[l]["dθ"])
				push!(y2,abs(dθ-line_ref[l]["dθ"])/line_ref[l]["dθ"])
				push!(y3,fl-line_ref[l]["fl"])
				push!(y4,abs(fl-line_ref[l]["fl"])/line_ref[l]["fl"])
			end

			impact[kl] = Dict{String,Any}("dθ_max" => maximum(y1),
							  "dθ_med" => median(y1),
							  "dθ_min" => minimum(y1),
							  "dθr_max" => maximum(y2),
							  "dθr_med" => median(y2),
							  "dθr_min" => minimum(y2),
							  "fl_max" => maximum(y3),
							  "fl_med" => median(y3),
							  "fl_min" => minimum(y3),
							  "flr_max" => maximum(y4),
							  "flr_med" => median(y4),
							  "flr_min" => minimum(y4),
							  "converged" => true)
		end		
	end

	return impact
end


function run_Nm1_own()
#TODO
end

## Remove line between nodes i and j.
function rm_line(ij::Tuple{Int64,Int64}, B::Array{Float64,2}, G::Array{Float64,2})

end

## Remove line #k, where numbering is given by 1:(1,2), 2:(1,3), ... n-1:(1,n), n:(2,3),...
function rm_line(k::Int64, B::Array{Float64,2}, G::Array{Float64,2})
	
end

## Extract the indices (pairs ij) and numbers (k) of the existing lines, i.e., M_ij ≠ 0.
function line_numbers(M::Array{Float64,2}, Zro::Float64=1e-8)
	n,nn = size(M)

	if n != nn
		@warn "Matrix is not square!"
	end

	m = 0
	k = 0

	ijs = Array{Tuple{Int64,Int64},1}()
	ks = Array{Int64,1}()

	for i in 1:n-1
		for j in i+1:n
			k += 1
			if abs(M[i,j]) > Zro
				m += 1
				push!(ijs,(i,j))
				push!(ks,k)
			end
		end
	end

	return ijs,ks,m
end

function ij2k(ij::Tuple{Int64,Int64}, n::Int64)
	i,j = ij

	if j <= i
		@warn "Invalid pair (j < i)!"
	end
	if i > n || j > n
		@warn "Invalid index or indices (too large)!"
	end

	nn = n - i + 1
	k = n*(n-1)/2 - nn*(nn-1)/2 + j - i
	
	return Int(k)
end

function k2ij(k::Int64, n::Int64)
	if k > n*(n-1)/2
		@warn "Invalid number (k too large)!"
	end

	N = 0
	c = n
	i = 0

	while N < k
		c -= 1
		N += c
		i += 1
	end

	x = N-n+i
	y = k-x
	j = i+y

	return (Int(i),Int(j))
end


