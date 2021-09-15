using DelimitedFiles

include("ksakaguchi.jl")
include("tools.jl")



function freq_width(L::Array{Float64,2}, ω0::Array{Float64,1}, θ0::Array{Float64,1}, α::Float64, C::Array{Array{Int64,1},1}, verb::Bool=false, res::Float64=.0005)
	if norm(ω0) < 1e-8
		@info "The frequency vector is close to zero."
	end

	γmax = π/2 - α

	max_iter = 100000
	thresh = 1e-8

	B,w = L2B(L)

	n = length(θ0)
	
	ω = ω0 .- mean(ω0)
	ω /= norm(ω)

	q0 = winding(θ0,C)
	q = Inf
	γ = -0.2
	it = 0
	Δmax = 0.

	θ1 = Array{Float64,1}()
	θ2 = Array{Float64,1}()
	x = Array{Float64,1}()

	while q != q0 && γ < 5. && it < max_iter
		γ += .2
		x = ksakaguchi(L,γ*ω,θ0,α,false,false,.01,thresh,max_iter)
		θ1 = x[1][:,end]
		θ2 = copy(θ1)
		Δmax = maximum(abs.(dcc(B'*θ1)))
		q = winding(θ2,C)
		it = x[4]
	end

	if (q != q0 && it >= max_iter) || (q != q0 && γ >= 5.)
			return NaN,NaN
	end

	β = γ
	dβ = 1.

	while dβ > res && Δmax < γmax
		if verb
			@info "β = $β, dβ = $dβ"
		end
		
		q = q0
		it = 0
		while q == q0 && it < 100000 && Δmax < γmax
			θ1 = copy(θ2)
			β += dβ
			x = ksakaguchi(L,β*ω,θ2,α,false,false,.01,thresh,max_iter)
#			ts,x = ksakaguchi_ND(L,β*ω,θ,α,(0.,10.))
			θ2 = mod.(x[1][:,end] .+ π,2π) .- π
			Δmax = maximum(abs.(dcc(B'*θ2)))
			q = winding(θ2,C)
			it = x[4]

			if verb
				@info "q = $q, it = $it"
			end
		end
		θ2 = copy(θ1)
		Δmax = maximum(abs.(dcc(B'*θ1)))
		β -= dβ
		dβ /= 10
	end

	βint = copy(β)

	dβ = .1

	while dβ > res
		if verb
			@info "β = $β, dβ = $dβ"
		end

		q = q0
		it = 0
		while q == q0 && it < 100000
			θ1 = copy(θ2)
			β += dβ
			x = ksakaguchi(L,β*ω,θ2,α,false,false,.01,thresh,max_iter)
#			ts,x = ksakaguchi_ND(L,β*ω,θ,α,(0.,10.))
			θ2 = mod.(x[1][:,end] .+ π,2π) .- π
			q = winding(θ2,C)
			it = x[4]

			if verb
				@info "q = $q, it = $it"
			end
		end
		θ2 = copy(θ1)
		β -= dβ
		dβ /= 10
	end


	βmax = copy(β)
	fmax = mean(x[2][:,end])


	if γ > 1e-6
		β = γ
		dβ = -.2

		while abs(dβ) > res && β > 0.
			if verb
				@info "β = $β, dβ = $dβ"
			end
	
			q = q0
			it = 0
			while q == q0 && it < max_iter && β > 0.
				θ1 = copy(θ2)
				β += dβ
				x = ksakaguchi(L,β*ω,θ2,α,false,false,.01,thresh,max_iter)
	#			ts,x = ksakaguchi_ND(L,β*ω,θ,α,(0.,10.))
				θ2 = mod.(x[1][:,end] .+ π,2π) .- π
				q = winding(θ2,C)
				it = x[4]
				
				if verb
					@info "q = $q, it = $it"
				end
			end
			θ2 = copy(θ1)
			Δ1 = dcc(B'*θ1)
			β -= dβ
			dβ /= 10

			βint = max(βint,β*(maximum(abs.(Δ1)) < γmax))
		end
		
		βmin = maximum(β,0.)
		fmin = mean(x[2][:,end])
	else
		βmin = 0.
		fmin = 0. # WRONG, BUT NO NEED FOR NOW...
	end
	
	return βmax,βint,βmin,fmax,fmin
end




