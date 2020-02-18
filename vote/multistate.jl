include("scripts_mini.jl")
include("big_rand.jl")

# #= 
Nmin = 201
Nmax = 501

emin = 0.
emax = .5
neps = 20

nS = 9

#d = 0.
d = .5
s = .2
D = 0.

w0 = .1
# =#

function run_multistates(nS::Int64, Nmm::Tuple{Int64,Int64}, emms::Tuple{Float64,Float64,Int64}, distri::Tuple{Float64,Float64,Float64}, w0::Float64=.1)
	Ns = rand((Nmm[1]:2:Nmm[2]),nS)
	x0s = Array{Array{Float64,1},1}()

	e_connect = Array{Float64,1}()

	d,s,D = distri
	for i in 1:nS
		push!(x0s, big_rand(Ns[i],-d/2+D,s,d/2+D,s))
		xxx = sort(x0s[end])
		dx = xxx[2:end] - xxx[1:end-1]
		push!(e_connect,maximum(dx))
	end
	
	X0 = Array{Float64,1}()
	for i in 1:nS
		X0 = [X0;x0s[i]]
	end

	e1 = max(emms[1],maximum(e_connect))
	epss = LinRange(e1,emms[2],emms[3])
	
	xi1 = Array{Float64,1}()
	xi2 = Array{Float64,1}()
	xi3 = Array{Float64,1}()
	nsti = Array{Int64,1}()
	os1 = Array{Float64,2}(undef,nS,0)
	os2 = Array{Float64,2}(undef,nS,0)
	ms1 = Array{Float64,2}(undef,nS,0)
	ms2 = Array{Float64,2}(undef,nS,0)

	figure()

	for eps in epss
		@info "$(100*round((eps-e1)/(emms[2]-e1),digits=3))% (ε = $(round(eps,digits=3)))"

		x1,x2,id1,id2,o1,o2,m1,m2,dds1,dds2 = multistates_efforts(Ns,x0s,eps,w0)
		@info "Unified country"
		x3,o3,xx = influence_effort_mini_repr(X0,eps,nS,w0)

		push!(xi1,x1)
		push!(xi2,x2)
		push!(xi3,abs(x3))
		
		oo1 = zeros(nS)
		oo1[id1] = o1
		os1 = [os1 oo1]
		oo2 = zeros(nS)
		oo2[id2] = o2
		os2 = [os2 oo2]
		
		mm1 = zeros(nS)
		mm1[id1] = m1
		ms1 = [ms1 mm1]
		mm2 = zeros(nS)
		mm2[id2] = m2
		ms2 = [ms2 mm2]

		
		h1 = 0
		for i in id1
			h1 += 1
			subplot(4,1,2)
			PyPlot.plot(eps,h1,"o",color="C$(i-1)")
		end
		h2 = 0
		for i in id2
			h2 += 1
			subplot(4,1,2)
			PyPlot.plot(eps,h2,"x",color="C$(i-1)")
		end
		ylabel("nsti")
	end

	subplot(4,1,1)
	PyPlot.plot(epss,xi1,"r",label="smallest state")
	PyPlot.plot(epss,xi2,"--b",label="smallest majority")
	PyPlot.plot(epss,xi3,"-k",label="unified country")
	ylabel("ξ")
	legend()
	title("n = $(length(X0)), n_mi = $(minimum(Ns)), n_ma = $(maximum(Ns)), n_st = $nS, δ = $d, σ = $s, Δ = $D")

	for i in 1:nS
		subplot(4,1,3)
		PyPlot.plot(epss,os1[i,:],color="C$(i-1)")
		PyPlot.plot(epss,os2[i,:],"--",color="C$(i-1)")
		ylabel("outcome")

		subplot(4,1,4)
		PyPlot.plot(epss,ms1[i,:],color="C$(i-1)")
		PyPlot.plot(epss,ms2[i,:],"--",color="C$(i-1)")
		ylabel("dm")
		xlabel("ε")
	end
end


function multistates_efforts(Ns::Array{Int64,1}, x0s::Array{Array{Float64,1},1}, eps::Float64, w0::Float64=.1)
	nS = length(Ns)

	LDis = Array{Array{Float64,2},1}()
	win = Array{Int64,1}()
	os = Array{Int64,1}()
	ms = Array{Int64,1}()
	dds = Array{Array{Float64,1},1}()

	for i in 1:nS
		A = Float64.((0 .< abs.(repeat(x0s[i],1,Ns[i]) - repeat(x0s[i]',Ns[i],1)) .<= eps))
		d = vec(sum(A,dims=1))
		push!(dds,d)
		push!(ms,Int(sum(d)/2))
		D = diagm(0 => d)
		LpD = Symmetric(2*D - A)
		LDi = inv(LpD)*Diagonal(D)

		x = LDi*x0s[i]
		o0,p0,n0 = outcome(x)

		push!(LDis,LDi)
		push!(win,sign(o0))
		push!(os,o0)
	end

	nsti = ceil(Int64,abs(sum(win)/2))

	idp = Int.(setdiff((win .> 0).*(1:nS),[0,]))
	idn = Int.(setdiff((1:nS),idp))

# influencing the smallest states
	@info "Influencing the smallest states"
	if sign(sum(win)) > 0
		id1 = Int.(sortslices([Ns[idp] idp],dims=1)[:,2])
	else
		id1 = Int.(sortslices([Ns[idn] idn],dims=1)[:,2])
	end

	xi_tot_1 = 0.

	os1 = Array{Float64,1}()

	for i in 1:nsti
		j = id1[i]
		n = Ns[j]
		x0 = x0s[j]
		LDi = LDis[j]

		w = zeros(n)
		x = LDi*(x0 + w)
		o0,p0,n0 = outcome(x)

		if o0 < 0
			x0 = -x0
			x = LDi*(x0 + w)
			o0,p0,n0 = outcome(x)
		end

		push!(os1,o0)
		
		o1 = copy(o0)
#		ids = mini_sort(x0, (o0 > 0))
		ids = mini_sort(x, (o0 > 0))

		while o1 > 0.
			c = 0
			while o1 > 0. && c < n
				c += 1
				w[ids[c]] -= w0
				x = LDi*(x0 + w)
				o1,p1,n1 = outcome(x)
			end
		end

		xi_tot_1 += sum(abs.(w))
	end

# influencing the smallest majority
	@info "Influencing the smallest majorities"

	if sign(sum(win)) > 0
		id2 = Int.(sortslices([os[idp] idp],dims=1)[:,2])
	else
		id2 = Int.(sortslices([abs.(os[idn]) idn],dims=1)[:,2])
	end

	xi_tot_2 = 0.

	os2 = Array{Float64,1}()

	for i in 1:nsti
		j = id2[i]
		n = Ns[j]
		x0 = x0s[j]
		LDi = LDis[j]

		w = zeros(n)
		x = LDi*(x0 + w)
		o0,p0,n0 = outcome(x)

		if o0 < 0
			x0 = -x0
			x = LDi*(x0 + w)
			o0,p0,n0 = outcome(x)
		end

		push!(os2,o0)

		o1 = copy(o0)
#		ids = mini_sort(x0, (o0 > 0))
		ids = mini_sort(x, (o0 > 0))

		while o1 > 0.
			c = 0
			while o1 > 0. && c < n
				c += 1
				w[ids[c]] -= w0
				x = LDi*(x0 + w)
				o1,p1,n1 = outcome(x)
			end
		end

		xi_tot_2 += sum(abs.(w))
	end

	return xi_tot_1, xi_tot_2, id1[1:nsti], id2[1:nsti], os1, os2, ms[id1[1:nsti]], ms[id2[1:nsti]], dds[id1[nsti]], dds[id2[nsti]]
end






