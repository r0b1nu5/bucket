using DataDrivenDiffEq, ModelingToolkit, LinearAlgebra, DataDrivenSparse, LinearAlgebra, PyPlot, Combinatorics

# This version infers an hypergraph with knowledge of the states 'X' and of the derivatives 'Y'.
# The maximal order of interaction we are looking at is 'd'. 
# ooi: orders of interest, i.e., orders of interaction that we want to identify.

function hyper_inf(X::Matrix{Float64}, Y::Matrix{Float64}, ooi::Vector{Int64}, dmax::Int64, thr_glob::Float64=1.)
	n,T = size(X)

	@variables x[1:n]

	problem = DirectDataDrivenProblem(X,Y,name = :HyperInference)

	prebasis = polynomial_basis([x[i] for i in 1:n],dmax)
	basis = Basis(prebasis,[x[i] for i in 1:n])

	res = solve(problem,basis,STLSQ())
	coeff = res.out[1].coefficients

	idx_o = Dict{Int64,Vector{Int64}}()
	agents_o = Dict{Int64,Vector{Vector{Int64}}}()
	inf_o = Dict{Int64,Vector{Vector{Int64}}}()
	Ainf = Dict{Int64,Any}()
	for o in ooi
		Ainf[o] = Dict{Tuple{Int64,Vector{Int64}},Float64}()
		idx_o[o],agents_o[o] = get_idx_o(o,x,prebasis)
		inf_o[o] = Vector{Int64}[]
		for i in 1:length(idx_o[o])
			id = idx_o[o][i]
			agents = agents_o[o][i]
			for a in agents
				Ainf[o][(a,agents)] = coeff[a,id]
			end
		end
	end

	return Ainf, coeff, idx_o, agents_o
end
#=
			a = coeff[agents,id]
			if sum(abs.(a) .> thr_glob) >= o
				as = sort(abs.(a),rev=true)
				da = as[1:end-1]-as[2:end]
				thr = (as[o] + as[o+1])/2
				δ = as[o] - as[o+1]
				inf = Int64.(setdiff((abs.(a) .> thr).*(1:n),[0.,]))
				push!(inf_o[o],inf)
				@info "The $o nodes most involved with $(basis[id]) are $inf (with margin δ = $δ)."
			end
		end
		if length(inf_o[o]) == 0
			@info "No $o-hyperedge has been inferred."
		end
	end
	
	return inf_o, coeff, idx_o
end
=#

# In case we want to infer only one order of hyperedge.

function hyper_inf(X::Matrix{Float64}, Y::Matrix{Float64}, ooi::Int64, dmax::Int64, thr_glob::Float64=1.)
    return hyper_inf(X,Y,[ooi,],dmax,thr_glob)
end



function get_idx_o(o::Int64, x::Symbolics.Arr{Num,1}, prebasis::Vector{Num})
	comb = combinations(1:length(x))
	idx = Int64[]
	agents = Vector{Int64}[]
	for c in comb
		mon = get_monomial(x,c)
		m,i = findmax(isequal.(prebasis,mon))
		if m > .1
			push!(idx,i)
			push!(agents,c)
		end
	end
	return idx,agents
end

function get_monomial(x::Symbolics.Arr{Num,1},c::Vector{Int64})
	mon = 1
	for i in c
		m *= x[i]
	end
	return mon
end

function get_monomial(x::Symbolics.Arr{Num,1},o::Int64)
    n = length(x)
    mon = Num[]
    comb = collect(combinations(1:n,o))

    for c in combi
	    push!(mon,get_monomial(x,c))
    end

    return mon
end


function check_inference(A2::Matrix{Float64}, A3::Array{Float64,3}, inf_o::Dict{Int64,Vector{Vector{Int64}}})
    n = size(A2)[2]
    A2t = zeros(n,n)
    A3t = zeros(n,n,n)

    if 2 in keys(inf_o)
	    e2 = [0 1;1 0.]
	    for id in inf_o[2]
	        A2t[id,id] = e2
	    end
	
	    p2 = 0
	    fp2 = 0
	    n2 = 0
	    fn2 = 0
	    for i in 1:n-1
	        for j in i+1:n
	            if A2t[i,j] > .1
	                p2 += 1
	                if A2[i,j] < .1
	                    fp2 += 1
	                end
	            else
	                n2 += 1
	                if A2[i,j] > .1
	                    fn2 += 1
	                end
	            end
	        end
	    end
	    tp2 = p2 - fp2
	    tn2 = n2 - fn2
	
	    sen2 = tp2/(tp2 + fn2)
	    spe2 = tn2/(tn2 + fp2)
	
	    @info "Detection of 2-edges: sensitivity = $sen2, specificity = $spe2."
    else
	    (sen2,spe2,tp2,fp2,tn2,fn2) = (NaN,NaN,NaN,NaN,NaN,NaN)
    end
    if 3 in keys(inf_o)
	    e3 = zeros(3,3,3)
	    e3[1,2,3] = e3[1,3,2] = e3[2,1,3] = e3[2,3,1] = e3[3,1,2] = e3[3,2,1] = 1.
	    for id in inf_o[3]
	        A3t[id,id,id] = e3
	    end
	    
	    p3 = 0
	    fp3 = 0
	    n3 = 0
	    fn3 = 0
	    for i in 1:n-2
	        for j in i+1:n-1
	            for k in j+1:n
	                if A3t[i,j,k] > .1
	                    p3 += 1
	                    if A3[i,j,k] < .1
	                        fp3 += 1
	                    end
	                else
	                    n3 += 1
	                    if A3[i,j,k] > .1
	                        fn3 += 1
	                    end
	                end
	            end
	        end
	    end
	    tp3 = p3 - fp3
	    tn3 = p3 - fn3
	
	    sen3 = tp3/(tp3 + fn3)
	    spe3 = tn3/(tn3 + fp3)
	
	    @info "Detection of 3-edges: sensitivity = $sen3, specificity = $spe3."
    else
	    (sen3,spe3,tp3,fp3,tn3,fn3) = (NaN,NaN,NaN,NaN,NaN,NaN)
    end

    return (sen2,spe2,tp2,fp2,tn2,fn2), (sen3,spe3,tp3,fp3,tn3,fn3)
end




