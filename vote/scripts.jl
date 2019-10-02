using Distributions, PyPlot, Statistics, LinearAlgebra, Random

function influence_effort_rand(n::Int64, eps::Float64, Di::Distribution)
	x0 = rand(Di,n)
	
	A = Float64.((0 .< abs.(repeat(x0,1,n) - repeat(x0',n,1)) .< eps))
	L = diagm(0 => vec(sum(A,dims=1))) - A
	
	xr = zeros(n)
	x = consensus(L,x0,zeros(n),ones(n),ones(n))
	o0,p0,n0 = outcome(x)

	s = 1.

	if o0 > 0.
		s = -1.
	end
	
	o1 = o0
	ids = Array(1:n)
	
	while s*o1 < 0.
		ids = randperm(n)
		c = 0
		while s*o1 < 0. && c < n
			c += 1
			xr[ids[c]] += s

			x = consensus(L,x0,xr,ones(n),ones(n))
			o1,p1,n1 = outcome(x)
		end
	end

	return sum(xr), o0
end

function influence_effort_rand(n::Int64, eps::Float64, x0::Array{Float64,1})
	A = Float64.((0 .< abs.(repeat(x0,1,n) - repeat(x0',n,1)) .< eps))
	L = diagm(0 => vec(sum(A,dims=1))) - A
	
	xr = zeros(n)
	x = consensus(L,x0,zeros(n),ones(n),ones(n))
	o0,p0,n0 = outcome(x)

	s = 1.

	if o0 > 0.
		s = -1.
	end
	
	o1 = o0
	ids = Array(1:n)
cc = 0

	while s*o1 < 0.
		cc += 1
#		@info "$(o1)"
#		@info "$(xr[1])"

		ids = randperm(n)
		c = 0
		while s*o1 < 0. && c < n
			c += 1
			xr[ids[c]] += s
#=
@info "$o1"
@info "$c"
@info "$(ids[c])"
@info "$(xr[ids[c]])"
@info "-----"
=#
			x = consensus(L,x0,xr,ones(n),ones(n))
			o1,p1,n1 = outcome(x)
		end
	end

	return sum(xr), o0
end

function loops1(n_run::Int64, n::Int64, epss::Array{Float64,1}, Di::Distribution)
	m = length(epss)
	
	effort = zeros(n_run,m)

	for j in 1:m
		@info "j = $j/$m"
		for i in 1:n_run
			effort[i,j] = influence_effort_rand(n,epss[j],Di)[1]
		end
	end

	return effort
end

function loops2(n_run::Int64, n::Int64, epss::Array{Float64,1}, Di::Distribution)
	m = length(epss)

	effort = zeros(n_run,m)

	x0 = rand(Di,n)

	for j in 1:m
		@info "j = $j/$m"
		for i in 1:n_run
			effort[i,j] = influence_effort_rand(n,epss[j],x0)[1]
		end
	end

	return effort,x0
end

function loops3(n_run::Int64, n::Int64, eps::Float64, Di::Distribution)
	effort = Array{Float64,1}()
	oc = Array{Float64,1}()

	for i in 1:n_run
		@info "i = $i/$n_run"
		e,o = influence_effort_rand(n,eps,Di)
		push!(effort,e)
		push!(oc,o)
	end

	return effort,oc
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

function consensus(L::Array{Float64,2}, x0::Array{Float64,1}, xr::Array{Float64,1}, a::Array{Float64,1}, f::Array{Float64,1})
	return inv(L + diagm(0 => a))*(diagm(0 => a)*x0 + diagm(0 => f)*xr)
end

function outcome(x::Array{Float64,1})
	return sum(sign.(x)), sum(x .> 0.), sum(x .< 0.)
end

function plot_quants(q0::Array{Float64,1}, q25::Array{Float64,1}, q50::Array{Float64,1}, q75::Array{Float64,1}, q00::Array{Float64,1}, epss::Array{Float64,1}, colo::String="C0") 
	for i in 1:length(epss)
		PyPlot.plot(epss[i]*[1,1],[q25[i],q75[i]],color=colo)
	end

	PyPlot.plot(epss,q0,"x",color=colo)
	PyPlot.plot(epss,q00,"x",color=colo)
	PyPlot.plot(epss,q50,"o",color=colo)
end






