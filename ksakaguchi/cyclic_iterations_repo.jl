using DelimitedFiles, LinearAlgebra

include("toolbox_repo.jl")


function iterations(Δ0::Vector{Float64}, B::Matrix{Float64}, C::Matrix{Float64}, u::Vector{Int64}, ω::Vector{Float64}, h::Union{Function,Vector{Function}}, γ::Union{Tuple{Float64,Float64},Vector{Tuple{Float64,Float64}}}, δ::Float64, s::Union{Float64,Vector{Float64}}=1., max_iter::Int64=100, tol::Float64=1e-6, verb::Bool=false)
	n,m2 = size(B)
	m = Int(2*m2)

	if length(γ[1]) == 1
		hγ = [(h(γ[1]),h(γ[2])),]
	else
		hγ = [(h[i](γ[i][1]),h[i](γ[i][2])) for i in 1:length(γ)]
	end

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

		Δ = Sδ(Δ,ω,B,Bout,P,Ld,δ,h,γ,hγ,s)
		err = norm(Δ - Δs[:,end])

		Δs = [Δs Δ]
	end
	
	return Δ, Δs
end




function Sδ(Δ::Vector{Float64}, ω::Vector{Float64}, B::Matrix{Float64}, Bout::Matrix{Float64}, P::Matrix{Float64}, W::Matrix{Float64}, δ::Float64, h::Union{Function,Vector{Function}}, γ::Union{Tuple{Float64,Float64},Vector{Tuple{Float64,Float64}}}, hγ::Vector{Tuple{Float64,Float64}}, s::Union{Float64,Vector{Float64}}=1.)
	return Δ - δ*B'*W*(Bout*hs([Δ;-Δ],h,γ,hγ,s) - ω)
end



function hs(x::Float64, h::Function, γ::Tuple{Float64,Float64}, hγ::Tuple{Float64,Float64}, s::Float64=1.)
	return (s*(x - γ[1]) + hγ[1])*(x < γ[1]) + h(x)*(γ[1] <= x <= γ[2]) + (s*(x - γ[2]) + hγ[2])*(x > γ[2])
end

function hs(x::Vector{Float64}, h::Vector{Function}, γ::Vector{Tuple{Float64,Float64}}, hγ::Vector{Tuple{Float64,Float64}}, s::Union{Float64,Vector{Float64}})
	if length(s) == 1
		return [hs(x[i],h[i],γ[i],hγ[i],s[1]) for i in 1:length(x)]
	else
		return [hs(x[i],h[i],γ[i],hγ[i],s[i]) for i in 1:length(x)]
	end
end

function hs(x::Vector{Float64}, h::Function, γ::Tuple{Float64,Float64}, hγ::Vector{Tuple{Float64,Float64}}, s::Float64)
	return [hs(x[i],h,γ,hγ[1],s[1]) for i in 1:length(x)]
end
