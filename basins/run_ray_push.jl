include("kuramoto.jl")
include("splay_states.jl")
include("cycle.jl")

n = 83
q = 5

R = 2π*rand(-1:1,n)

qmax = round(Int64,n/4)

splays = splay_states(n)

L = cycle(n)
ω = zeros(n)

nstep = 500
αs = LinRange(0,1,nstep)
qs = Vector{Int64}()
its = Vector{Int64}()

θq = splays[q]

dl = Array{Float64,1}()

for i in 1:nstep
	α = αs[i]
	
	if i%20 == 0
		@info "$(round(100*α))%"
	end

	θ,it = kuramoto(L,ω,θq + α*R)

	push!(qs,winding(θ,Array(1:n)))
	push!(its,it)
end

figure()
subplot(1,2,1)
PyPlot.plot(αs,qs)
PyPlot.plot(αs[2:end],qs[2:end]-qs[1:end-1])
xlabel("α")
ylabel("q,dq")

subplot(1,2,2)
PyPlot.hist(qs,(minimum(qs)-.25:.5:maximum(qs)+.25),density=true)
xlabel("q")


