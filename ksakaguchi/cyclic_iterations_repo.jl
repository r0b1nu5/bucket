using DelimitedFiles, LinearAlgebra

include("toolbox_repo.jl")

# ================================================================================
"""
	iterations(Δ0::Vector{Float64}, B::Matrix{Float64}, C::Matrix{Float64}, u::Vector{Int64}, ω::Vector{Float64}, h::Union{Function,Vector{Function}}, γ::Union{Tuple{Float,Float64},Vector{Tuple{Float64,Float64}}}, δ::Float64, s::Union{Float64,Vector{Float64}}=1., max_iter::Int64=100, tol::Float64=1e-6, verb::Bool=false)

Runs the iteration scheme described in [Delabays et al. (2022)] for cyclic networks of diffusively coupled oscillators. Starts at initial conditions `Δ0`, which is projected on the afine subspace ker(B)^perp + 2π*pinv(`C`)*`u`. The script runs for at most `max_iter` iterations or until the correction is smaller than `tol`.

_INPUT_:\\
`Δ0`: Initial conditions of the iterations. Each component should ideally be bounded by the corresponding components of `γ`. \\
`B`: Incidence matrix of the (undirected) graph. \\
`C`: Cycle-edge incidence matrix associated to the cycle basis of the graph (see [Delabays et al. (2022)]. \\
`u`: Winding vector of the cell where the solution is searched. \\
`ω`: Vector of natural frequencies. \\
`h`: Vector of coupling functions over the (bidirected) edges of the graph. If a single function is given, the couplings are assumed homogenous. \\
`γ`: Vector of tuples, composed of the lower (1st comp.) and upper (2nd comp.) bounds on the domain of `h`, such that it is strictly increasing. \\
`δ`: Scaling parameter. \\
`s`: Slope of the extended coupling functions. Should have the same dimension as `h`. \\
`max_iter`: Maximal number of iterations allowed. \\
`tol`: Minimal correction allowed between two iterations. \\
`verb`: If true, enumerates the iterations.

_OUTPUT_:\\
`Δ`: Final state at the end of the iterations. \\
`Δs`: Sequence of states along the iterations. 
"""
function iterations(Δ0::Vector{Float64}, B::Matrix{Float64}, C::Matrix{Float64}, u::Vector{Int64}, ω::Vector{Float64}, h::Union{Function,Vector{Function}}, γ::Union{Tuple{Float64,Float64},Vector{Tuple{Float64,Float64}}}, δ::Float64, s::Union{Float64,Vector{Float64}}=1., max_iter::Int64=100, tol::Float64=1e-6, verb::Bool=false)
	n,m2 = size(B)
	m = Int(2*m2)

	Bb = [B -B]
	Bout = Bb.*(Bb .> 0.)
	P = cycle_proj(B)
	
	Cdag = pinv(C)

	Ld = pinv(B*B')

	Δ = B'*Ld*B*Δ0 + 2π*Cdag*u
	Δs = copy(Δ)

	t = 1
	err = 1000.
	
	while t < max_iter && err > tol
		t += 1

		if verb
			@info "iter = $t, err = $err"
		end

		Δ = Sδ(Δ,ω,B,Bout,P,Ld,δ,h,γ,s)
		err = norm(Δ - Δs[:,end])

		Δs = [Δs Δ]
	end
	
	return Δ, Δs
end


# ================================================================================
"""
	Sδ(Δ::Vector{Float64}, ω::Vector{Float64}, B::Matrix{Float64}, Bout::Matrix{Float64}, P::Matrix{Float64}, W::Matrix{Float64}, δ::Float64, h::Union{Function,Vector{Function}}, γ::Union{Tuple{Float64,Float64},Vector{Tuple{Float64,Float64}}}, s::Union{Float64,Vector{Float64}}=1.)

Iteration functions whose fixed points are solutions to the Dissipative Network Flow problem [Delabays et al. (2022)]. 

_INPUT_:\\
`Δ`: Argument of the interation function. \\
`ω`: Vector of natural frequencies. \\
`B`: Incidence matrix of the (undirected) graph. \\
`Bout`: Out-incidence matrix of the bidirected graph. \\
`P`: Cycle projection matrix (see [Delabays et al. (2022)]. \\
`W`: Weight matrix to be tuned. In [Delabays et al. (2022)], we take the pseudoinverse of the Laplacian. \\
`δ`: Scaling parameter which, if small enough, guarantees `Sδ` to be contracting. \\
`h`: Vector of the (directed) coupling functions. If a single `Function` is given, the couplings are assumed homogeneous. \\
`γ`: Vector of tuples composed of the lower (1st comp.) and upper (2nd comp.) bounds on the domain of `h`, such that it is stricly increasing. Should have the same dimension as `h`. \\
`s`: Slope of the extended coupling functions, outside of their domain. 

_OUTPUT_:\\
`Δ2`: Updated value of the state `Δ`.
"""
function Sδ(Δ::Vector{Float64}, ω::Vector{Float64}, B::Matrix{Float64}, Bout::Matrix{Float64}, P::Matrix{Float64}, W::Matrix{Float64}, δ::Float64, h::Union{Function,Vector{Function}}, γ::Union{Tuple{Float64,Float64},Vector{Tuple{Float64,Float64}}}, s::Union{Float64,Vector{Float64}}=1.)
	return Δ - δ*B'*W*(Bout*hs([Δ;-Δ],h,γ,s) - ω)
end


# ================================================================================
"""
	hs(x::Float64, h::Function, γ::Tuple{Float64,Float64}, s::Float64=1.)

Extended coupling function, extending `h` to the whole real axis. Matches `h` on [`γ`[1],`γ`[2]], is continuous, and has slope `s` outside of [`γ`[1],`γ`[2]]. 

_INPUT_:\\
`x`: Argument of the extended coupling function. \\
`h`: Coupling function to be extended. \\
`γ`: Tuple of the two bounds of the domain of `h`. \\
`s`: Slope of the extended coupling function outside of the domain of `h`. \\

_OUTPUT_:\\
`hx`: h_s(x).
"""
function hs(x::Float64, h::Function, γ::Tuple{Float64,Float64}, s::Float64=1.)
	hg1 = h(γ[1])
	hg2 = h(γ[2])
	return (s*(x - γ[1]) + hg1)*(x < γ[1]) + h(x)*(γ[1] <= x <= γ[2]) + (s*(x - γ[2]) + hg2)*(x > γ[2])
end

"""
	hs(x::Vector{Float64}, h::Vector{Function}, γ::Vector{Tuple{Float64,Float64}}, s::Union{Float64,Vector{Float64}})

Computes the value of the extended coupling function elementwise on `x`, `h`, `γ`, and `s`. 
"""
function hs(x::Vector{Float64}, h::Vector{Function}, γ::Vector{Tuple{Float64,Float64}}, s::Union{Float64,Vector{Float64}})
	if length(s) == 1
		return [hs(x[i],h[i],γ[i],s[1]) for i in 1:length(x)]
	else
		return [hs(x[i],h[i],γ[i],s[i]) for i in 1:length(x)]
	end
end

"""
	hs(x::Vector{Float64}, h::Function, γ::Tuple{Float64,Float64}, s::Float64)

Computes the value of the extended coupling function elementwise on `x`.
"""
function hs(x::Vector{Float64}, h::Function, γ::Tuple{Float64,Float64}, s::Float64)
	return [hs(x[i],h,γ,s[1]) for i in 1:length(x)]
end



