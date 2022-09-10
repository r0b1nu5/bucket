using PyPlot, LinearAlgebra, DelimitedFiles

function oscillators(L::Array{Float64,2}, v0::Array{Float64,1}, vs::Array{Float64,1}, ps::Array{Float64,1}, qs::Array{Float64,1}, η::Float64, α::Float64, r::Float64=1., ℓ::Float64=.001, ω0::Float64=50., h::Float64=.01, max_iter::Int64=1000, tol::Float64=1e-6)
	nn = length(v0)
	n = Int64(nn/2)

	I2 = diagm(0 => ones(2))
	In = diagm(0 => ones(n))
	J = [0. -1.;1. 0.]
	κ = atan(ℓ/r*ω0*2π)

	LL = kron(L,I2)
	JJ = kron(In,J)
	KK = Kk(vs,ps,qs,κ)

	v1 = v0
	v2 = v0
	Vs = Array{Float64,2}(undef,nn,0)
	Vs = [Vs v1]

	err = 1000.
	iter = 0

	while err > tol && iter < max_iter
		iter += 1
		if iter%100 == 0
			@info "iter = $iter"
		end

		v1 = v2

		k1 = f(v1,vs,KK,LL,η,α)
		k2 = f(v1+h/2*k1,vs,KK,LL,η,α)
		k3 = f(v1+h/2*k2,vs,KK,LL,η,α)
		k4 = f(v1+h*k3,vs,KK,LL,η,α)

		dv = (k1 + 2*k2 + 2*k3 + k4)/6

		v2 = v1 + h*dv

		err = maximum(abs.(dv))

		Vs = [Vs v2]
	end

	return Vs
end


function rot(θ::Float64)
	return [cos(θ) -sin(θ);sin(θ) cos(θ)]
end

function Φk(v::Array{Float64,1},vs::Float64)
	return (vs - norm(v))/vs
end

function Φ(v::Array{Float64,1}, vs::Array{Float64,1})
	ϕ = [Φk(v[(2*i-1):(2*i)],vs[i]) for i in 1:length(vs)]

	return diagm(0 => kron(ϕ,[1,1]))
end

function Kk(vs::Float64, ps::Float64, qs::Float64, κ::Float64)
	return rot(κ)*[ps qs;-qs ps]/vs^2
end

function Kk(vs::Array{Float64,1}, ps::Array{Float64,1}, qs::Array{Float64,1}, κ::Float64)
	n = length(vs)

	KK = zeros(2*n,2*n)

	for i in 1:n
		KK[[2*i-1,2*i],[2*i-1,2*i]] = Kk(vs[i],ps[i],qs[i],κ)
	end

	return KK
end

function f(v::Array{Float64,1}, vs::Array{Float64,1}, KK::Array{Float64,2}, LL::Array{Float64,2}, η::Float64, α::Float64)
	Φv = Φ(v,vs)

	return η*(KK - LL)*v + α*Φv*v
end

function P(θs::Array{Float64,1}, vs::Array{Float64,1}, L::Array{Float64,2}, ℓ::Float64=.001, r::Float64=1., ω0::Float64=50.)
	b,w = L2B(L)
	B = [b -b]
	Bout = B.*(B .> 0.)
	Bin = -B.*(B .< 0.)

	dθ = B'*θs

	p = vs.^2 .*Bout*(r .- (Bout'*vs).*(Bin'*vs).*(r*cos.(dθ) - ℓ*ω0*sin.(dθ)))/(r^2 + ω0^2*ℓ^2)

	return p
end

function Q(θs::Array{Float64,1}, vs::Array{Float64,1}, L::Array{Float64,2}, ℓ::Float64=.001, r::Float64=1., ω0::Float64=50.)
	b,w = L2B(L)
	B = [b -b]
	Bout = B.*(B .> 0.)
	Bin = -B.*(B .< 0.)

	dθ = B'*θs

	q = vs.^2 .*Bout*(ω0*ℓ .- (Bout'*vs).*(Bin'*vs).*(ω0*ℓ*cos.(dθ) + r*sin.(dθ)))/(r^2 + ω0^2*ℓ^2)

	return q
end


function L2B(L::Array{Float64,2})
	n = size(L)[1]
	
	B = Array{Float64,2}(undef,n,0)
	w = Array{Float64,1}(undef,0)
	for i in 1:n-1
		for j in i+1:n
			if L[i,j] != 0.0
				ed = zeros(n)
				ed[i] = 1.0
				ed[j] = -1.0
				B = [B ed]
				push!(w,-L[i,j])
			end
		end
	end
	
	return B,w
end
