using PyPlot, Statistics

include("gen-hypergraph.jl")
include("hyper-kuramoto.jl")
include("tools.jl")

n = 10
p = .2
n_iter = 1
n_α = 50
α_max = 20

figure("Condition")
figure("Condition")
PyPlot.plot([0,1.05],[0,0],"--k")
xlabel("α/α_star")
ylabel("λ2")

αs = LinRange(0,α_max,n_α)
Λ1 = zeros(n_α,0)
Λ2 = zeros(n_α,0)

θ0 = zeros(n)
xmin = 0.
xmax = 0.
ymin = 0.
ymax = 0.

for i in 1:n_iter
	A,B,E = rand_3_graph(n,p)
	W = diagm(0 => repeat(A[:,4],inner=3))
	ω = rand(n); ω .-= mean(ω)
	
	α_star = 1/maximum(inv(W)*pinv(B)*ω)
	@info "α* = $(α_star)"
	global αs = LinRange(0,α_star*1.01,n_α)

	λ1 = Float64[]
	λ2 = Float64[]
	q05 = Float64[]
	q25 = Float64[]
	q50 = Float64[]
	q75 = Float64[]
	q95 = Float64[]

	for α in αs
		@info "α/α* = $(α/α_star)"
		Θs,dΘs = hyper_k(zeros(0,3),A,α*ω,θ0,0.,0.,.01,100000,1e-3)
		if size(Θs)[2] > 1000
			v = vec(maximum(dΘs[:,end-1000:end],dims=1) - minimum(dΘs[:,end-1000:end],dims=1))
		else
			v = vec(maximum(dΘs,dims=1) - minimum(dΘs,dims=1))
		end
		push!(q05,quantile(v,.05))
		push!(q25,quantile(v,.25))
		push!(q50,quantile(v,.5))
		push!(q75,quantile(v,.75))
		push!(q95,quantile(v,.95))
			
		global θ0 = Θs[:,end]
		J = get_jac(θ0,A)
		λs = eigvals(J)
		push!(λ1,sort(real.(λs))[end])
		push!(λ2,sort(real.(λs))[end-1])
 #=
		if λ2[end] > -1e-4
			t = now()
			writedlm("data/A-$(α/α_star)-$(t).csv",A,',')
			writedlm("data/W-$(α/α_star)-$(t).csv",W,',')
			writedlm("data/E-$(α/α_star)-$(t).csv",B,',')
			writedlm("data/ω-$(α/α_star)-$(t).csv",ω,',')
		end
# =#

	end
	global Λ1 = [Λ1 λ1]
	global Λ2 = [Λ2 λ2]

	figure("Condition")
#	PyPlot.plot(αs./α_star,Λ1[:,end],":k")
	PyPlot.plot(αs./α_star,Λ2[:,end],"--o")
	global xmin = minimum([xmin;αs./α_star])
	global xmax = maximum([xmax;αs./α_star])
	global ymin = minimum([ymin;λ2])
	global ymax = maximum([ymax;λ2])

	figure("Sync?")
	PyPlot.plot(αs./α_star,q05,":",color="C$(mod(i-1,10))")
	PyPlot.plot(αs./α_star,q95,":",color="C$(mod(i-1,10))")
	PyPlot.plot(αs./α_star,q25,"-",color="C$(mod(i-1,10))")
	PyPlot.plot(αs./α_star,q75,"-",color="C$(mod(i-1,10))")
	PyPlot.plot(αs./α_star,q50,"-",color="C$(mod(i-1,10))",lw=2)
end


figure("Condition")
PyPlot.plot([1,1],[ymin,ymax],"--k")


