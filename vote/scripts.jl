using Distributions, PyPlot, Statistics, LinearAlgebra, Random


function influence_effort_rand(x0::Array{Float64,1}, eps::Float64, w0::Float64=1.)
	n = length(x0)

	A = Float64.((0 .< abs.(repeat(x0,1,n) - repeat(x0',n,1)) .< eps))
	D = diagm(0 => vec(sum(A,dims=1)))
	L = D - A
	Di = diagm(0 => 1 ./ diag(D))
	LDi = inv(Di*L + diagm(0 => ones(n)))
#	LIi = inv(L + diagm(0 => ones(n)))

	xr = zeros(n)
	x = LDi*(x0 + xr)
	o0,p0,n0 = outcome(x)
	
	if o0 < 0.
		x0 = -x0
		x = LDi*(x0 + xr)
		o0,p0,n0 = outcome(x)
	end
	xx = copy(x)

	o1 = o0
	ids = Array(1:n)

	while o1 > 0.
		ids = randperm(n)
		c = 0
		while o1 > 0. && c < n
			c += 1
			xr[ids[c]] -= w0
			x = LDi*(x0 + xr)
			o1,p1,n1 = outcome(x)
		end
	end

	return sum(xr), o0, xx
end

function influence_effort_fiedler(x0::Array{Float64,1}, eps::Float64, w0::Float64=1., a::Int64=2)
	n = length(x0)

	A = Float64.((0 .< abs.(repeat(x0,1,n) - repeat(x0',n,1)) .< eps))
	D = diagm(0 => vec(sum(A,dims=1)))
	L = D - A
	Di = diagm(0 => 1 ./ diag(D))
	LDi = inv(Di*L + diagm(0 => ones(n)))
#	LIi = inv(L + diagm(0 => ones(n)))

	xr = zeros(n)
	x = LDi*(x0 + xr)
	o0,p0,n0 = outcome(x)
	
	if o0 < 0.
		x0 = -x0
		x = LDi*(x0 + xr)
		o0,p0,n0 = outcome(x)
	end
	xx = copy(x)

	o1 = o0
	ids = Array(1:n)

	while o1 > 0.
		ids = fiedler_sort(x0,eps,a)
		c = 0
		while o1 > 0. && c < n
			c += 1
			xr[ids[c]] -= w0
			x = LDi*(x0 + xr)
			o1,p1,n1 = outcome(x)
		end
	end

	return sum(xr), o0, xx
end

function influence_effort_mini(x0::Array{Float64,1}, eps::Float64, w0::Float64=1.)
	n = length(x0)

	A = Float64.((0 .< abs.(repeat(x0,1,n) - repeat(x0',n,1)) .< eps))
	D = diagm(0 => vec(sum(A,dims=1)))
	L = D - A
	Di = diagm(0 => 1 ./ diag(D))
	LDi = inv(Di*L + diagm(0 => ones(n)))
#	LIi = inv(L + diagm(0 => ones(n)))

	xr = zeros(n)
	x = LDi*(x0 + xr)
	o0,p0,n0 = outcome(x)

	if o0 < 0.
		x0 = -x0
		x = LDi*(x0 + xr)
		o0,p0,n0 = outcome(x)
	end
	xx = copy(x)

	o1 = o0
	ids = Array(1:n)

	while o1 > 0.
		ids = mini_sort(x0)
		c = 0
		while o1 > 0. && c < n
			c += 1
			xr[ids[c]] -= w0
			x = LDi*(x0 + xr)
			o1,p1,n1 = outcome(x)
		end
	end

	return sum(xr), o0, xx
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

function loop_rand(n_run::Int64, x0::Array{Float64,1}, epss::Array{Float64,1})
	m = length(epss)

	effort = zeros(n_run,m)

	for j in 1:m
		@info "j = $j/$m"
		for i in 1:n_run
			effort[i,j] = influence_effort_rand(x0, epss[j])[1]
		end
	end

	return effort
end

function loop_fiedler(x0::Array{Float64,1}, epss::Array{Float64,1}, modes::Array{Int64,1}=[2,])
	m = length(epss)
	x0 = sort(x0)

	effort = zeros(length(modes),m)

	for j in 1:m
		@info "j = $j/$m"
		for i in 1:length(modes)
			effort[i,j],o0 =  influence_effort_fiedler(x0, epss[j], modes[i])
		end
	end

	return effort
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
	for i in 1:length(epss)
		PyPlot.plot(epss[i]*[1,1],[q25[i],q75[i]],color=colo)
	end

	PyPlot.plot(epss,q0,"x",color=colo)
	PyPlot.plot(epss,q00,"x",color=colo)
	PyPlot.plot(epss,q50,"o",color=colo)
end

function plot_fiedler(effort::Array{Float64,2}, epss::Array{Float64,1}, colos::Array{String,1}=["C3",], modes::Array{Int64,1}=[2,])
	for i in 1:length(modes)
		PyPlot.plot(epss,effort[i,:],color=colos[i])
	end
end

function plot_mini(effort::Array{Float64,1}, epss::Array{Float64,1})
	PyPlot.plot(epss,effort,color="C1")
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

function consensus(L::Array{Float64,2}, x0::Array{Float64,1}, xr::Array{Float64,1}, a::Array{Float64,1}, f::Array{Float64,1})
	return inv(Symmetric(L + diagm(0 => a)))*(diagm(0 => a)*x0 + diagm(0 => f)*xr)
end

function outcome(x::Array{Float64,1})
	return sum(sign.(x)), sum(x .> 0.), sum(x .< 0.)
end


function fiedler_sort(x0::Array{Float64,1}, eps::Float64, a::Int64=2)
	si = sign.(x0)
	n = length(x0)

	dx = x0[2:end] - x0[1:end-1]

	ids = Int.(setdiff((dx .> eps).*(1:n-1),[0,]))
	clusts = [[1;ids.+1] [ids;n] [ids;n]-[0;ids]]

	todo = (n+1)*ones(length(ids)+1)
	order = Array{Int64,1}()

	while maximum(todo) > 0
		m,j = findmax(min.(todo,si[clusts[:,2]].*clusts[:,3]))
		todo[j] = -(n+1)

		if clusts[j,3] == 1
			order = [order;clusts[j,1]]
		else
			nn = clusts[j,3]
			
			x = x0[clusts[j,1]:clusts[j,2]]
			A = Float64.((0 .< abs.(repeat(x,1,nn) - repeat(x',nn,1)) .< eps))
			D = diagm(0 => vec(sum(A,dims=1)))
			L = Symmetric(D - A)
			Di = diagm(0 => 1 ./ diag(D))

			ssi = sign.(x)
	
			ei = eigen(Di*L + diagm(0 => ones(n)))
			us = ei.vectors
			lsi = sortslices([real.(ei.values) 1:nn],dims=1)
			ls = lsi[:,1]
			li = Int.(lsi[:,2])
		
			l = -100.
			i = 0
			id = 0
			f = 2
		
			while l < 1 + 1e-8 && i < nn
				i += 1
				l = ls[i]
			end
			
			while f < a && i < nn
				i += 1
				f += 1
			end
		
			if i == length(ls) && ls[i] < 1e-8
				@info "WARNING: Graph is highly disconnected, probably no good Fiedler strategy."
			end
		
			l = ls[i]
			id = li[i]
			
			uf = norm.(us[:,id])

# Brute version
#			uu = sortslices([ssi.*uf 1:nn],dims=1,rev=true)
# Normalized version (in order to discard extremes
			uu = sortslices([(ssi.*uf)./(abs.((1:nn) .- nn/2)) 1:nn],dims=1,rev=true)
 			
			
			sir = sign(uf[Int(uu[1,2])])
		
			su = sir .* sign.(uf)
			
			uu = sortslices([su.*abs.(uf) 1:nn],dims=1,rev=true)
		
			test = true
		
			if uu[1,1] < 1e-8
				@info "WARNING: no positive component in the signed Fielder mode!"
				test = false
			end

			order = [order;uu[:,2] .+ clusts[j,1] .- 1]
		end
	end
	
	return Int.(order)
end

function mini_sort(x0::Array{Float64,1})
	n = length(x0)

	ip = setdiff((x0 .> 0).*Array(1:n),[0,])
	xp = x0[ip]

	it = setdiff((x0 .<= 0).*Array(1:n),[0,])
	xn = x0[it]

	idp = Int.(sortslices([xp ip],dims=1)[:,2])
	idn = Int.(sortslices([abs.(xn) it],dims=1)[:,2])

	return [idp;idn]
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







