using PyPlot, LinearAlgebra

include("big_rand.jl")
include("Leps.jl")

n = 2001

x0 = big_rand(n,-.5,.2,.5,.2)

emin = maximum(x0[2:end]-x0[1:end-1])

epss = LinRange(emin+1e-4,.8,20)

ii = Array{Int64,1}()
for f in [.05,.1,.15,.2,.3,.4,.5]
	push!(ii,findmin(abs.(x0 .- f))[2])
end

Xih = Array{Float64,2}(undef,length(ii),0)
Xs = Array{Float64,2}(undef,n,0)
nps = Array{Float64,2}(undef,length(ii),0)
nns = Array{Float64,2}(undef,length(ii),0)

for eps in epss
	global Xih,mp,mn,Xs,nps,nns
	@info "$(round(eps/maximum(epss),digits=3)*100)%"

	L,A,D = Leps(x0,eps)
	d = diag(D)

	Di = Diagonal(diagm(0 => 1 ./ d))
	
	LDi = inv(Symmetric(L + D))*Diagonal(D)

	xs = LDi*x0
	Xs = [Xs xs]

	idp = setdiff((xs .> 0.).*(1:n),[0,])
	idn = setdiff(Array(1:n),idp)

	np = Array{Float64,1}()
	nn = Array{Float64,1}()
	for i in ii
		push!(np,sum(A[i,i:end]))
		push!(nn,sum(A[i,1:i]))
	end
	nps = [nps np]
	nns = [nns nn]

	wA = diagm(0 => d[ii])*A[ii,:].*(repeat(xs[ii],1,n) - repeat(xs',length(ii),1))
	dx0 = xs[ii] - x0[ii]

	sdx = (sum(abs.(wA),dims=2) + abs.(dx0))/2

	Xih = [Xih sdx]
end

figure()
subplot(3,1,1)
for i in 1:length(ii)
	PyPlot.semilogy(epss,Xih[i,:],"o",color="C$(i-1)",label="x0 = $(round(x0[ii[i]],digits=3))")
end
legend()

subplot(3,1,2)
for i in 1:length(ii)
	PyPlot.plot(epss,Xs[ii[i],:] .- x0[ii[i]],"-o")
end

subplot(3,1,3)
for i in 1:length(ii)
	PyPlot.plot(epss,nps[i,:]./nns[i,:],"-",color="C$(i-1)")
#	PyPlot.plot(epss,nns[i,:],"--",color="C$(i-1)")
end





