include("kuramoto.jl")
include("splay_states.jl")

n = 83
N = 100

qmax = round(Int64,n/4)

splays = splay_states(n)

L = spdiagm(0 => 2*ones(n)) - spdiagm(1 => ones(n-1)) - spdiagm(-1 => ones(n-1)) - spdiagm(n-1 => ones(1)) - spdiagm(1-n => ones(1))
ω = zeros(n)

qs = [11,12,13]

αs = LinRange(0,20,500)

dloc = Dict{Int64,Array{Float64,1}}()

for q in qs
	θ0 = splays[q]

	dl = Array{Float64,1}()

	for i in 1:N
		x = 2*rand(n) .- 1
		ϵ = x/norm(x)

		θ = θ0
		iter = 0

		while winding(θ,Array(1:n)) == q
			iter += 1
			α = αs[iter]

			θ,it = kuramoto(L,ω,θ0 + α*ϵ)
		end

		push!(dl,αs[iter])
	end

	dloc[q] = dl
end



