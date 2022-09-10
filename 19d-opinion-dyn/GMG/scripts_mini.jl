using Distributions, PyPlot, Statistics, LinearAlgebra, Random, SparseArrays

include("result.jl")

function influence_effort_mini(x0::Array{Float64,1}, eps::Float64, m::Int64, w0::Float64=.1)
	n = length(x0)

	A = Float64.((0 .< abs.(repeat(x0,1,n) - repeat(x0',n,1)) .< eps))
	d = vec(sum(A,dims=1))
	D = diagm(0 => d)
	L = D - A
	LpD = Symmetric(2*D - A)
	Di = diagm(0 => 1 ./ d)
	LDi = inv(LpD)*Diagonal(D)

	w = zeros(n)
	x = LDi*x0
	np,margp,nn,margn = result_margin(x,m)

	if np < nn
		x0 = -x0
		x = -x
		np,margp,nn,margn = result_margin(x,m)
	end

	xx = copy(x)

	o1 = o0
@info "Mini: computing..."
	ids = mini_sort(x0, true)
@info "Mini: computed."
	while o1 > 0.
		c = 0
		while o1 > 0. && c < n
			c += 1
			w[ids[c]] -= w0
			x -= w0*LDi[:,ids[c]]
			o1,p1,n1 = outcome(x)
		end
	end

	return sum(w), o0, xx
end

function influence_effort_mini_repr(x0::Array{Float64,1}, eps::Float64, n_repr::Int64, w0::Float64=.1)
	n = length(x0)

	A = Float64.((0 .< abs.(repeat(x0,1,n) - repeat(x0',n,1)) .< eps))
	d = vec(sum(A,dims=1))
	D = diagm(0 => d)
	L = D - A
	LpD = Symmetric(2*D - A)
	Di = diagm(0 => 1 ./ d)
	LDi = inv(LpD)*Diagonal(D)
#	LDi = inv(Di*L + diagm(0 => ones(n)))
#	LIi = inv(L + diagm(0 => ones(n)))

	w = zeros(n)
	x = LDi*(x0 + w)
	o0,p0,n0 = outcome_repr(x,n_repr)

	if o0 < 0
		x0 = -x0
		x = LDi*(x0 + w)
		o0,p0,n0 = outcome_repr(x,n_repr)
	end

	xx = copy(x)

	o1 = o0
	ids = mini_sort(x0, true)
	while o1 > 0.
		c = 0
		while o1 > 0. && c < n
			c += 1
			w[ids[c]] -= w0
			x = LDi*(x0 + w)
			o1,p1,n1 = outcome_repr(x,n_repr)
		end
	end

	return sum(w), o0, xx
end

function outcome(x::Array{Float64,1})
	return sum(sign.(x)), sum(x .> 0.), sum(x .< 0.)
end

function outcome_repr(x::Array{Float64,1}, n_repr::Int64=1)
	n_ppl = length(x)
	np = sum(x .> 0)
	nn = sum(x .< 0)

	rp = round(n_repr*np/n_ppl)
	rn = n_repr - rp

	return rp-rn, rp, rn
end

# pos2neg is true if we want to influence the outcome from positive to negative, false otherwise.

function mini_sort(x0::Array{Float64,1}, pos2neg::Bool)
	n = length(x0)

	ip = setdiff((x0 .> 0).*Array(1:n),[0,])
	xp = x0[ip]

	it = setdiff((x0 .<= 0).*Array(1:n),[0,])
	xn = x0[it]

	idp = Int.(sortslices([xp ip],dims=1)[:,2])
	idn = Int.(sortslices([abs.(xn) it],dims=1)[:,2])
	
	if pos2neg
		ids = [idp;idn]
	else
		ids = [idn;idp]
	end

	return ids
end

function loop_mini(x0::Array{Float64,1}, epss::Array{Float64,1})
	m = length(epss)

	effort = Array{Float64,1}()

	for j in 1:m
		@info "j = $j/$m"
		push!(effort,influence_effort_mini(x0,epss[j])[1])
	end
	
	return effort
end


function eps_connect(x0::Array{Float64,1}, epss::Array{Float64,1}, thr::Float64=.95)
	n = length(x0)

	dx = x0[2:end] - x0[1:end-1]
	ig = [minimum([1.;abs.(2*(setdiff((dx .> eps).*(1:n-1),[0,]) .- n/2)./n)]) for eps in epss]
	epi = minimum(setdiff((ig .> thr).*(1:length(epss)),[0,]))

	return epss[epi]
end
		
function quants(effort::Array{Float64,2})
	n,m = size(effort)
	
	q0 = Array{Float64,1}()
	q25 = Array{Float64,1}()
	q50 = Array{Float64,1}()
	q75 = Array{Float64,1}()
	q00 = Array{Float64,1}()

	for j in 1:m
		push!(q0, quantile(effort[:,j], 0.))
		push!(q25, quantile(effort[:,j], 0.25))
		push!(q50, quantile(effort[:,j], 0.5))
		push!(q75, quantile(effort[:,j], 0.75))
		push!(q00, quantile(effort[:,j], 1.))
	end

	return q0,q25,q50,q75,q00
end

function plot_quartiles(effort::Array{Float64,2}, epss::Array{Float64,1}, colo::String="C0")
	q0,q25,q50,q75,q00 = quants(effort)

	plot_quants(q0,q25,q50,q75,q00,epss,colo)
end

function plot_mean(effort::Array{Float64,2}, epss::Array{Float64,1}, colo::String="C0")
	n,m = size(effort)

	mef = [mean(effort[:,i]) for i in 1:m]
	sef = [std(effort[:,i]) for i in 1:m]

	for i in 1:m
		PyPlot.plot([epss[i],epss[i]],[mef[i]-sef[i],mef[i]+sef[i]],color=colo)
	end

	PyPlot.plot(epss,mef,"o",color=colo)
end
	
function plot_quants(q0::Array{Float64,1}, q25::Array{Float64,1}, q50::Array{Float64,1}, q75::Array{Float64,1}, q00::Array{Float64,1}, epss::Array{Float64,1}, colo::String="C0") 
#=
for i in 1:length(epss)
		PyPlot.plot(epss[i]*[1,1],[q25[i],q75[i]],color=colo)
	end

	PyPlot.plot(epss,q0,"x",color=colo)
	PyPlot.plot(epss,q00,"x",color=colo)
	PyPlot.plot(epss,q50,"o",color=colo)
=#
	PyPlot.fill([epss;epss[length(epss):-1:1]],[q0;q00[length(q00):-1:1]],color=colo,alpha=.1)
	PyPlot.fill([epss;epss[length(epss):-1:1]],[q25;q75[length(q75):-1:1]],color=colo,alpha=.3)
	PyPlot.plot(epss,q50,"-o",color=colo)
end


function plot_mini(effort::Array{Float64,1}, epss::Array{Float64,1})
	PyPlot.plot(epss,effort,color="C1",label="minimum")
end


function clusterings(L::Array{Float64,2}, x0::Array{Float64,1})
	n = length(x0)
	
	idn = setdiff((x0 .< 0.).*(1:n),[0.,])
	idp = setdiff((x0 .> 0.).*(1:n),[0.,])
	
	nn = length(idn)
	np = length(idp)
	
	A = diagm(0 => diag(L)) - L
	An = A[idn,idn]
	Ap = A[idp,idp]
	
	cn = Array{Float64,1}()
	cp = Array{Float64,1}()

	for i in idn
		nei = setdiff(Int.(A[:,i].*(1:n)),[0;idp])
		nnei = length(nei)
		c = 0
		if nnei > 1
			for j in 1:nnei-1
				for k in j+1:nnei
					if A[nei[j],nei[k]] == 1.
						c += 1
					end
				end
			end
			push!(cn,c/(nnei*(nnei-1)/2))
		end
	end

	Cn = mean(cn)

	for i in idp
		nei = setdiff(Int.(A[:,i].*(1:n)),[0;idn])
		nnei = length(nei)
		c = 0
		if nnei > 1
			for j in 1:nnei-1
				for k in j+1:nnei
					if A[nei[j],nei[k]] == 1.
						c += 1
					end
				end
			end
			push!(cp,c/(nnei*(nnei-1)/2))
		end
	end

	Cp = mean(cp)

	C0 = (sum(A)/2 - sum(An)/2 - sum(Ap)/2)/(nn*np)

	return Cn,Cp,C0
end







