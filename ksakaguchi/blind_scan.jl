using DelimitedFiles, PyPlot, LinearAlgebra
include("L2B.jl")
include("tools.jl")
include("ksakaguchi.jl")

n_store = 1000
n_iter = 1000
α = .5
γ = π/2

L = diagm(0 => 2*ones(10)) - diagm(1 => ones(9)) - diagm(-1 => ones(9)) - diagm(9 => ones(1)) - diagm(-9 => ones(1))
C = Array(1:10)
T = Array(1:9)
C_Σ = ones(1,10)

u = [0,]

Bu,w = L2B(L)
B = [Bu -Bu]
Bout = (abs.(B) + B)/2
Bin = Bout - B
n = size(B)[1]
W = diagm(0 => [w;w])
BoW = Bout*W

θr = Array{Float64,2}(undef,n,0)
ωr = Array{Float64,2}(undef,n,0)
for i in 1:n_store
	θt = sample_winding_cell(u,Bu,T,C_Σ)
	while cohesiveness_inc(θt,Bu) > γ
		θt = sample_winding_cell(u,Bu,T,C_Σ)
	end
	global θr = [θr θt]
	global ωr = [ωr BoW*sin.(B'*θr[:,end] .- α)]
end

R = Array{Float64,1}()
best_pair_θ = zeros(n,2)
best_pair_ω = zeros(n,2)
best_max = 0.

c = 0

while c < n_iter
	global c +=1 

	θt = sample_winding_cell(u,Bu,T,C_Σ)
	while cohesiveness_inc(θt,Bu) > γ
		θt = sample_winding_cell(u,Bu,T,C_Σ)
	end
	θ = θt
	ω = BoW*sin.(B'*θ .- α)

	r = Array{Float64,1}()

	for i in 1:n_store
		dθ = norm(mod.(θ - θr[:,i] .+ π,2π) .- π)
		dω = norm(ω - ωr[:,i])
		push!(r,dθ/dω)
	end

	rm,id = findmax(r)
	push!(R,rm)
	if rm > best_max
		global best_max = rm
		global best_pair_θ = [θ θr[:,id]]
		global best_pair_ω = [ω ωr[:,id]]
	end
end

PyPlot.hist(R,50)

θ1 = best_pair_θ[:,1]
θ2 = best_pair_θ[:,2]
ω1 = best_pair_ω[:,1]
ω2 = best_pair_ω[:,2]

x1 = ksakaguchi(L,ω1,θ2,α)
x2 = ksakaguchi(L,ω2,θ1,α)

θf1 = mod.(x1[1][:,end] .+ π,2π) .- π
θf2 = mod.(x2[1][:,end] .+ π,2π) .- π

@info "($(winding(θ1,C)),$(winding(θf1,C)))"
@info "($(winding(θ2,C)),$(winding(θf2,C)))"



#=
while c < n_store
	c += 1

	θ = 2π*rand(n)
	while cohesiveness_inc(θ,B) > γ
		θ = 2π*rand(n)
	end

	qs = winding(θ,C)

	ω = BoW*sin.(B'*θ .- α)
=#






