using LinearAlgebra, DelimitedFiles

include("L2B.jl")

function kuramoto(L::Matrix{Float64}, ω::Vector{Float64}, θ0::Vector{Float64}, verb::Bool=false, h::Float64=.01, iter::Int64=10000, tol::Float64=0.)
	n = length(ω)
	B,w = L2B(L)
	W = diagm(0 => w)

	θ = copy(θ0)
	θs = copy(θ0)

	t = 0
	c = 0
	err = 1000.

	while t < iter && err > tol
		t += 1

		if t%1000 == 0
			if verb
				@info "K: t = $t"
			end

			c += 1
			writedlm("temp/th_$c.csv",θs[:,1:end-1],',')
			θs = θs[:,end]
		end

		k1 = ω - B*W*sin.(B'*θ)
		k2 = ω - B*W*sin.(B'*(θ + h*k1/2))
		k3 = ω - B*W*sin.(B'*(θ + h*k2/2))
		k4 = ω - B*W*sin.(B'*(θ + h*k3))

		dθ = (k1 + 2*k2 + 2*k3 + k4)/6
		θ += h*dθ
		θs = [θs θ]

		err = maximum(abs.(dθ))
	end
	
	Θs = Matrix{Float64}(undef,n,0)
	for d in 1:c
		Θs = [Θs readdlm("temp/th_$d.csv",',')]
		rm("temp/th_$d.csv")
	end
	Θs = [Θs θs]

	return Θs
end

function complex_k1(L::Matrix{Float64}, ω::Vector{Float64}, z0::Vector{Complex{Float64}}, win::Int64=100, verb::Bool=false, h::Float64=.01, iter::Int64=10000, tol::Float64=0.)
	n = length(ω)
	B,w = L2B(L)
	W = diagm(0 => w)

	z = copy(z0)
	zs = copy(z0)

	t = 0
	c = 0
	err = 1000.

	while t < iter && err > tol
		t += 1

		if t%win == 0
			if verb
				@info "cK1: t = $t"
			end

			c += 1
			writedlm("temp/zr_$c.csv",real.(zs[:,1:end-1]),',')
			writedlm("temp/zi_$c.csv",imag.(zs[:,1:end-1]),',')
			zs = real.(zs[:,end]) .+ im*0.
			z = real.(z) .+ im*0.
		end

		k1 = ω - B*W*sin.(B'*z) + im*abs.(B)*W*cos.(B'*z)
		k2 = ω - B*W*sin.(B'*(z + h*k1/2)) + im*abs.(B)*W*cos.(B'*(z + h*k1/2))
		k3 = ω - B*W*sin.(B'*(z + h*k2/2)) + im*abs.(B)*W*cos.(B'*(z + h*k2/2))
		k4 = ω - B*W*sin.(B'*(z + h*k3)) + im*abs.(B)*W*cos.(B'*(z + h*k3))

		dz = (k1 + 2*k2 + 2*k3 + k4)/6
		z += h*dz
		zs = [zs z]

		err = maximum(abs.(dz))
	end

	Zs = Matrix{Complex{Float64}}(undef,n,0)
	for d in 1:c
		Zs = [Zs (readdlm("temp/zr_$d.csv",',') + im*readdlm("temp/zi_$d.csv",','))]
		rm("temp/zr_$d.csv")
		rm("temp/zi_$d.csv")
	end
	Zs = [Zs zs]

	return Zs
end
		


