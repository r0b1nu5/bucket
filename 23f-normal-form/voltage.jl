using PyPlot

function dyn_η(η0::Union{Complex{Float64},Vector{Complex{Float64}}}, F::Any, max_iter::Int64=1000, thr::Float64=1e-6, h::Float64=.01)
	iter = 0
	err = 1000.

	aν,aζ,aζb,Lκ = F

	η = η0
	ηs = η0

	while iter < max_iter && err > thr
		iter += 1
		if iter%100 == 0
			@info "iter: $iter/$max_iter"
		end

		k1 = fη_lin(η,aν,aζ,aζb,Lκ)
		k2 = fη_lin(η + k1*h/2,aν,aζ,aζb,Lκ)
		k3 = fη_lin(η + k2*h/2,aν,aζ,aζb,Lκ)
		k4 = fη_lin(η + k3*h,aν,aζ,aζb,Lκ)

		dη = (k1+2*k2+2*k3+k4)/6
		η += h*dη
		ηs = [ηs η]

		err = maximum(abs.(dη))
	end

	return ηs
end

function fη_lin(η, aν::Complex{Float64}, aζ::Complex{Float64}, aζb::Complex{Float64}, Lκ::Union{Complex{Float64},Matrix{Complex{Float64}}})
	return aν*(η + conj.(η)) - aζ*Lκ*conj.(η) - aζb*conj.(Lκ)*η
end





