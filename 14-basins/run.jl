using Statistics

include("kuramoto.jl")
include("splay_states.jl")

N = 1000
n = 83
qmax = round(Int64,n/4)

splays = splay_states(n)

L = spdiagm(0 => 2*ones(n)) - spdiagm(1 => ones(n-1)) - spdiagm(-1 => ones(n-1)) - spdiagm(n-1 => ones(1)) - spdiagm(1-n => ones(1))

ω = zeros(n)

qs = Array{Int64,1}()
θ0s_vec = Array{Array{Float64,1},1}()
θ0s = Dict{Int64,Array{Array{Float64,1},1}}(q => Array{Array{Float64,1},1}() for q in -qmax:qmax)
ds = Dict{Int64,Array{Float64,1}}(q => Array{Float64,1}() for q in -qmax:qmax)
d2s = Dict{Int64,Array{Float64,1}}(q => Array{Float64,1}() for q in -qmax:qmax)

for i in 1:N
	if i%100 == 0
		@info "==============================="
		@info "i = $i"
	end

	θ0 = 2π*rand(n)

	θs,it = kuramoto(L,ω,θ0)

	push!(qs,winding(θs,Array(1:n)))

	push!(θ0s[qs[end]],θ0)
	push!(θ0s_vec,θ0)
	
	push!(ds[qs[end]],dist(θ0,splays[qs[end]]))
	push!(d2s[qs[end]],dist2splay(θ0,qs[end]))
end








