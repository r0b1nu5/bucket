using Statistics, LinearAlgebra, Dates

include("ksakaguchi.jl")

function winding(θ::Array{Float64,1}, C::Array{Int64,1})
	if length(θ) < 1
		q = []
	else
		θ1 = θ[C]
		θ2 = θ[[C[2:end];C[1]]]

		dθ = θ1 - θ2

		q = round(Int,sum(mod.(dθ .+ π,2π) .- π)/(2π))
	end

	return q
end

# Computes the winding vector for cycles defined by the cycles in C.
function winding(θ::Array{Float64,1}, C::Array{Array{Int64,1},1})
	if length(C) < 1 || length(θ) < 1
		qs = []
	else
		qs = [winding(θ,C[i]) for i in 1:length(C)]
	end

	return qs
end

# Computes the max (in abs value) angle difference, using the incidence matrix.
function cohesiveness_inc(θ::Array{Float64,1}, B::Array{Float64,2})
	return maximum(abs.(mod.(transpose(B)*θ .+ π,2π) .- π))
end

# Computes the max (in abs value) angle difference, using the adjacency/Laplacian matrix.
function cohesiveness_adj(θ::Array{Float64,1}, A::Array{Float64,2})
	n = length(θ)
	
	dθ = (θ*ones(1,n) - ones(n)*θ').*(abs.(A .* (1 .- diagm(0 => ones(n)))) .> 1e-8)

	return maximum(dθ)
end

# Draws a random element of the winding cell of u
# u: Winding vector
# B: Incidence matrix (n x m)
# T: Indices of the edges composing a spanning tree
# C: cycle-edge incidence matrix (see Jafarpour et al. (2020))
# TODO Are we really sampling the whole winding cell???
function sample_winding_cell(u::Array{Int64,1}, B::Array{Float64,2}, T::Array{Int64,1}, C::Array{Float64,2}, max_attempt::Int64=1000)
	n,m = size(B)

	Bt = transpose(B)
	Ttd = pinv(transpose(B[:,T]))
	Cd = pinv(C)

	c = 0
	test = true
	θ = Array{Float64,1}()

	while c < max_attempt && test
		c += 1

		x = 2π*rand(n) .- π
		x .-= mean(x)

		y = Bt*x + 2π*Cd*u

		if norm(y,Inf) < π
			Dθ = Ttd*y[T]
			θ = θ .- mean(θ)
			test = false
		end
	end

	@info "Number of attempts: $c"
	if c == max_attempt
		@info "WARNING: No sample found."
	end

	return θ
end

function jacobian(L::Array{Float64,2}, θ::Array{Float64,1}, α::Float64)
	n = length(θ)

	A = -L.*(1 .- diagm(0 => ones(n)))

	dθ = θ*ones(1,n) - ones(n)*θ'

	J = A.*cos.(dθ .- α)
	J = J - diagm(0 => vec(sum(J,dims=2)))

	return J
end

function freq_width(L::Array{Float64,2}, ω0::Array{Float64,1}, θ0::Array{Float64,1}, α::Float64, C::Array{Int64,1}, verb::Bool=false, res::Float64=.0005)
	if norm(ω0) < 1e-8
		@info "The frequency vector is close to zero."
	end

	n = length(θ0)
	
	ω = ω0 .- mean(ω0)
	ω /= norm(ω)

	x = ksakaguchi(L,zeros(n),θ0,α,true,false,.01,1e-6)
#	ts,x = ksakaguchi_ND(L,zeros(n),θ0,α,(0.,10.))
	θ1 = x[1][:,end]
	θ = copy(θ1)
	q0 = winding(θ,C)

	β = 0.
	dβ = 1.

	while dβ > res
		if verb
			@info "β = $β, dβ = $dβ"
		end

		q = q0
		it = 0
		while q == q0 && it < 100000
			β += dβ
			x = ksakaguchi(L,β*ω,θ,α,true,false,.01,1e-6)
#			ts,x = ksakaguchi_ND(L,β*ω,θ,α,(0.,10.))
			q = winding(x[1][:,end],Array(1:n))
			it = x[4]
			
			if verb
				@info "q = $q, it = $it"
			end
		end
		θ = x[1][:,1]
		β -= dβ
		dβ /= 10
	end

	βmax = copy(β)
	fmax = mean(x[2][:,end])

	β = 0.
	dβ = 1.

	while dβ > res
		if verb
			@info "β =$β, dβ = $dβ"
		end

		q = q0
		it = 0
		while q == q0 && it < 100000
			β -= dβ
			x = ksakaguchi(L,β*ω,θ,α,true,false,.01,1e-6)
#			ts,x = ksakaguchi_ND(L,β*ω,θ,α,(0.,10.))
			q = winding(x[1][:,end],Array(1:n))
			it = x[4]

			if verb
				@info "q = $q, it = $it"
			end
		end
		θ = x[1][:,1]
		β += dβ
		dβ /= 10
	end

	βmin = copy(β)
	fmin = mean(x[2][:,end])

	return βmin,βmax,fmin,fmax
end


function cyqle(n::Int64, q::Int64=1)
	A = zeros(n,n)
	for i in 1:q
		A += diagm(i => ones(n-i)) + diagm(-i => ones(n-i)) + diagm(n-i => ones(i)) + diagm(i-n => ones(i))
	end

	D = diagm(0 => vec(sum(A,dims=2)))

	L = D - A

	return L
end

function gershgorin(M::Array{Float64,2})
	c = diag(M)

	n = length(c)

	M0 = (1 .- diagm(0 => ones(n))).*M

	rh = vec(sum(abs.(M0),dims=2))
	rv = vec(sum(abs.(M0),dims=1))

	ls = eigvals(M)

	t = LinRange(0,2pi,100)
	x1 = minimum(real.(c)-rh)-1
	x2 = maximum(real.(c)+rh)+1
	y1 = minimum(imag.(c)-rh)-1
	y2 = maximum(imag.(c)+rh)+1
	x3 = minimum(real.(c)-rv)-1
	x4 = maximum(real.(c)+rv)+1
	y3 = minimum(imag.(c)-rv)-1
	y4 = maximum(imag.(c)+rv)+1

	figure()
	subplot(1,2,1)
	PyPlot.plot([x1,x2],[0.,0.],"--k")
	PyPlot.plot([0.,0.],[y1,y2],"--k")
	axis([x1,x2,y1,y2])
	subplot(1,2,2)
	PyPlot.plot([x3,x4],[0.,0.],"--k")
	PyPlot.plot([0.,0.],[y3,y4],"--k")
	axis([x3,x4,y3,y4])

	for i in 1:n
		subplot(1,2,1)
		PyPlot.plot(rh[i]*cos.(t) .+ real(c[i]),rh[i]*sin.(t) .+ imag(c[i]),label="i = $i")
		PyPlot.plot(real(ls[i]),imag(ls[i]),"ok")
		subplot(1,2,2)
		PyPlot.plot(rv[i]*cos.(t) .+ real(c[i]),rv[i]*sin.(t) .+ imag(c[i]))
		PyPlot.plot(real(ls[i]),imag(ls[i]),"ok")
	end

	subplot(1,2,1)
	title("Row Gershgorin circles")
	legend()
	subplot(1,2,2)
	title("Column Gershgoring circles")

	return c,rh,rv
end




