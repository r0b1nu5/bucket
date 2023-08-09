using LinearAlgebra, Statistics, SpecialFunctions, Distributions, DelimitedFiles, StatsBase, PyPlot

function L2B(L::Matrix{Float64})
	n = size(L)[1]

	Bout = zeros(n,0)
	Bin = zeros(n,0)
	w = Float64[]

	for i in 1:n
		for j in 1:n
			if i !== j && abs(L[i,j]) > 1e-6
				Bout = [Bout I[1:n,i]]
				Bin = [Bin I[1:n,j]]
				push!(w,-L[i,j])
			end
		end
	end

	return Bout,Bin,w
end

function winding(θ::Vector{Float64}, σ::Vector{Int64})
	Δ = θ[[σ[2:end];σ[1]]] - θ[σ]
	return round(Int64,sum(mod.(Δ .+ π,2π) .- π)/(2π))
end

function gen_cycle_undir(n::Int64)
	return diagm(0 => 2*ones(n)) - diagm(1 => ones(n-1)) - diagm(-1 => ones(n-1)) - diagm(n-1 => ones(1)) - diagm(1-n => ones(1))
end

function gen_cycle_dir(n::Int64)
	return diagm(0 => ones(n)) - diagm(1 => ones(n-1)) - diagm(1-n => ones(1))
end

function gof_normal(x::Vector{<:Any}, qmax::Int64, iter::Int64=5000)
	μ0 = mean(x)
	σ0 = std(x)
	n = length(x)

	d0 = ks_discr_normal(x,qmax,μ0,σ0)

	c = 0 

	for i in 1:iter
		y = round.(rand(Normal(μ0,σ0),n))
		μ1 = mean(y)
		σ1 = std(y)

		d1 = ks_discr_normal(y,qmax,μ1,σ1)
		
		if d1 > d0
			c +=1
		end
	end

	return c/iter
end

function ks_discr_normal(x::Vector{<:Any}, qmax::Int64, μ::Float64=0., σ::Float64=1.)
	cdft = [(.5+.5*erf((q+.5-μ)/(sqrt(2)*σ))) for q in -qmax:qmax]
	cdfe = [sum(x .<= q)/length(x) for q in -qmax:qmax]

	return norm(cdft-cdfe,Inf)
end

function gof_exp(x::Vector{<:Any}, qmax::Int64, iter::Int64=5000)
	μ0 = mean(x)

	ns = [sum(x .== q) for q in -qmax:qmax]

	β0 = mle_exp(ns,μ0,qmax)
	C0 = C_exp(β0,μ0,qmax)

	d0 = ks_exp(x,β0,μ0,qmax)

	c = 0 
	for i in 1:iter
		y = rand_exp(rand(length(x)),β0,μ0,C0,qmax)
		μ1 = mean(y)
		β1 = mle_exp([sum(y .== q) for q in -qmax:qmax],μ1,qmax)

		d1 = ks_exp(y,β1,μ1,qmax)

		if d1 > d0
			c += 1
		end
	end

	return c/iter

end

function mle_exp(ns::Vector{Int64}, μ::Float64, qmax::Int64)
	return mle_exp(Float64.(ns),μ,qmax)
end

function mle_exp(ns::Vector{Float64}, μ::Float64, qmax::Int64)
	l = sum(ns)
	
	side = 4.99
	β = 5.

	L = zeros(100)
	while side > 1e-4
		βs = LinRange(max(β-side,side/100),β+side,100)
		for i in 1:100
			L[i] = l*log(C_exp(βs[i],μ,qmax)) - βs[i]*sum(ns.*abs.((-qmax:qmax).-μ))
		end

		β = βs[findmax(L)[2]]
		side /= 10
	end

	return β
end

function C_exp(β::Float64, μ::Float64, qmax::Int64)
	return 1/sum(exp.(-β*abs.((-qmax:qmax) .- μ)))
end

function ks_exp(x::Vector{<:Any}, β::Float64, μ::Float64, qmax::Int64)
	pdft = [0.;C_exp(β,μ,qmax)*exp.(-β*abs.((-qmax:qmax) .- μ))]
	cdft = cumsum(pdft)
	cdfe = [0.;[sum(x .<= q)/length(x) for q in -qmax:qmax]]
	
	return norm(cdft-cdfe,Inf)
end

function rand_exp(yy::Vector{Float64}, β::Float64, μ::Float64, C::Float64, qmax::Int64)
	pdft = [0.;C_exp(β,μ,qmax)*exp.(-β*abs.((-qmax:qmax) .- μ))]
	cdft = cumsum(pdft)
	
	qs = [maximum((1:2*qmax+2).*(y .>= cdft)) for y in yy] .- qmax .- 1

	return qs
end

function qqplot(x::Vector{<:Any}, qmax::Int64, μ::Float64=0., σ::Float64=1.; color::String)
	cdft = [(.5+.5*erf((q+.5-μ)/(sqrt(2)*σ))) for q in -qmax:qmax]
	cdfe = [sum(x .<= q)/length(x) for q in -qmax:qmax]
	
	figure("Q-Q plot")
	PyPlot.plot(sort(cdft),sort(cdfe),".",color=color)
	xlabel("Theoretical")	
	ylabel("Empirical")
	
end

function init_Qs(n::Int64)
	qmax = floor(Int64,n/2)
	
	Qs = Dict{Int64,Dict{Int64,Int64}}()
	for qi in -qmax:qmax
		Qs[qi] = Dict{Int64,Int64}(qf => 0 for qf in -qmax:qmax)
	end
	
	return Qs
end

function save_Qs(Qs::Dict{Int64,Dict{Int64,Int64}}, file::String="./")
	spQ = zeros(Int64,0,3)

	for qi in keys(Qs)
		for qf in keys(Qs[qi])
			if Qs[qi][qf] > 0
				spQ = [spQ;[qi qf Qs[qi][qf]]]
			end
		end
	end

	writedlm(file*"Qs.csv",spQ,',')
end

function load_Qs(file::String)
	spQ = Int64.(readdlm(file,','))

	Qs = Dict{Int64,Dict{Int64,Int64}}()
	
	for l in 1:size(spQ)[1]
		qi = spQ[l,1]
		qf = spQ[l,2]
		N = spQ[l,3]

		if qi in keys(Qs)
			Qs[qi][qf] = N
		else
			Qs[qi] = Dict{Int64,Int64}(qf => N)
		end
	end

	return Qs
end

function load_Qs(folder::String, ids::Vector{Int64}, type::String, n::Int64, n_intra::Int64)
	Qs = Dict{Int64,Dict{Int64,Int64}}()

	for id in ids
		qs = load_Qs(folder*"$id-k-"*type*"-$n-$(n_intra)-Qs.csv")
		Qs = merge_Qs(Qs,qs)
	end

	return Qs
end

function merge_Qs(Q1::Dict{Int64,Dict{Int64,Int64}}, Q2::Dict{Int64,Dict{Int64,Int64}})
	Qs = Q1
	for qi in keys(Q2)
		for qf in keys(Q2[qi])
			if qi in keys(Q1)
				if qf in keys(Q1[qi])
					Qs[qi][qf] += Q2[qi][qf]
				else
					Qs[qi][qf] = Q2[qi][qf]
				end
			else
				Qs[qi] = Dict{Int64,Int64}(qf => Q2[qi][qf])
			end
		end
	end

	return Qs
end


function get_qfs(qi::Int64, Qs::Dict{Int64,Dict{Int64,Int64}})
	if qi in keys(Qs)
		qfs = sort(collect(keys(Qs[qi])))
		Ns = [Qs[qi][qf] for qf in qfs]
		return qfs,Ns
	else
		@info "No key $qi in Qs."
		return nothing
	end
end

function get_qis(qf::Int64, Qs::Dict{Int64,Dict{Int64,Int64}})
	qis = Int64[]
	Ns = Int64[]
	for qi in keys(Qs)
		if qf in keys(Qs[qi])
			push!(qis,qi)
			push!(Ns,Qs[qi][qf])
		end
	end

	if isempty(qis)
		@info "No key qi with $qf field in Qs."
		return nothing
	else
		X = sortslices([qis Ns],dims=1)
		return X[:,1],X[:,2]
	end
end

function get_stat_f(Qs)
	qis = sort(collect(keys(Qs)))

	Mf = Float64[]
	Σf = Float64[]
	Nf = Int64[]
	for qi in qis
		qfs,Ns = get_qfs(qi,Qs)
		qs = Int64[]
		for i in 1:length(qfs)
			q = qfs[i]
			qs = [qs;q*ones(Int64,Ns[i])]
		end

		push!(Mf,mean(qs))
		push!(Σf,std(qs))
		push!(Nf,length(qs))
	end

	return Mf,Σf,Nf,qis
end

function get_stat_i(Qs)
	qt = Int64[]
	for qi in keys(Qs)
		qt = [qt;collect(keys(Qs[qi]))]
	end
	qfs = sort(union(qt))

	Mi = Float64[]
	Σi = Float64[]
	Ni = Int64[]
	for qf in qfs
		qis,Ns = get_qis(qf,Qs)
		qs = Int64[]
		for i in 1:length(qis)
			q = qis[i]
			qs = [qs;q*ones(Int64,Ns[i])]
		end
		
		push!(Mi,mean(qs))
		push!(Σi,std(qs))
		push!(Ni,length(qs))
	end

	return Mi,Σi,Ni,qfs
end


function get_normal_par(Qs::Dict{Int64,Dict{Int64,Int64}})
	Mf,Σf,Nf,qis = get_stat_f(Qs)
	Mi,Σi,Ni,qfs = get_stat_i(Qs)

	# The std of our estimate of the mean scales as 1/sqrt(N) (CLT). We weight the linear regression with the inverse of the variance, i.e., with N.
	Mf_slope = sum(Mf.*qis.*Nf)/sum(Nf.*qis.^2)
	Mi_slope = sum(Mi.*qfs.*Ni)/sum(Ni.*qfs.^2)
	

	σf = mean(Σf[Nf.>1],weights(Nf[Nf.>1]))
	σi = mean(Σi[Ni.>1],weights(Ni[Ni.>1]))

	return Mf_slope,Mi_slope,σf,σi
end

function plot_Qs(Qs)
	Mf,Σf,Nf,qis = get_stat_f(Qs)
	Mi,Σi,Ni,qfs = get_stat_i(Qs)

	figure("Histograms")
	i = 0
	for qi in qis
		i += 1
		qs = sort(collect(keys(Qs[qi])))
		ns = Int64[]
		for q in qs
			push!(ns,Qs[qi][q])
		end
		if sum(ns) > 20
			subplot(1,3,1)
			PyPlot.plot(qs .- Mf[i],ns./sum(ns),"-o")
		end
	end
	subplot(1,3,1)
	xlabel("qf-μf")
	ylabel("P(qf|qi)")

	i = 0
	for qf in qfs
		i += 1
		qs = []
		for q in keys(Qs)
			if qf in keys(Qs[q])
				push!(qs,q)
			end
		end
		qs = sort(qs)
		ns = Int64[]
		for q in qs
			push!(ns,Qs[q][qf])
		end
		if sum(ns) > 20
			subplot(1,3,2)
			PyPlot.plot(qs .- Mi[i],ns./sum(ns),"-o")
		end
	end
	subplot(1,3,2)
	xlabel("qi-μi")
	ylabel("P(qi|qf)")

	xxx = get_normal_par(Qs)
	subplot(1,3,3)
	PyPlot.plot(qis,Mf,color="C0",label="μf")
	PyPlot.plot(qis,Σf,"--",color="C0",label="σf")
	PyPlot.plot(qis,xxx[1].*qis,":",color="C0")
	PyPlot.plot(qis[[1,length(qis)]],[xxx[3],xxx[3]],":",color="C0")
	PyPlot.plot(qfs,Mi,color="C1",label="μi")
	PyPlot.plot(qfs,Σi,"--",color="C1",label="σi")
	PyPlot.plot(qfs,xxx[2].*qfs,":",color="C1")
	PyPlot.plot(qfs[[1,length(qfs)]],[xxx[4],xxx[4]],":",color="C1")
	xlabel("q")
end	




