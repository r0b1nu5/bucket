using LinearAlgebra, DelimitedFiles, Distributions

# B2: directed node-edge incidence matrix (with 1's and -1's)
# B3: directed node-facet incidence matrix (with 2's and -1's)
# τ0: correlation time of the nodal noise
# ξ0: amplitude of the  noise
# ϕ2: phase frustration (for Kuramoto-Sakaguchi)
# ϕ3: phase frustration (for 3rd-order KS)
function hyper_k_incidence(B2::Matrix{Float64}, B3::Matrix{Float64}, ω::Vector{Float64}, θ0::Vector{Float64}, τ0::Float64, ξ0::Float64=1., ϕ2::Float64=0., ϕ3::Float64=0., a2::Union{Float64,Vector{Float64}}=1., a3::Union{Float64,Vector{Float64}}=1., h::Float64=.01, max_iter::Int64=10000, tol::Float64=1e-6)
	n = length(θ0)

	B2o = B2.*(B2 .> 0.)
	B3o = B3.*(B3 .> 0.)./2

	θ = θ0
	θs = θ0
	dθs = zeros(n,1)

	iter = 0
	err = 1000.
	c = 0

	ξ = rand(Normal(0.,1.),n)

	while iter<max_iter && err>tol
		iter += 1
		if iter%100 == 0
			@info "iter: $iter"
			c += 1
			writedlm("temp/ths_$c.csv",θs[:,1:end-1],',')
			θs = θs[:,end]
			writedlm("temp/dths_$c.csv",dθs[:,2:end],',')
			dθs = dθs[:,end]
		end
		
		ξ = [cnoise(ξ[i],τ0/h) for i in 1:n]

		k1 = ω - a2.*(B2o*(sin.(B2'*θ .- ϕ2) .+ sin(ϕ2))) - a3.*(B3o*(sin.(B3'*θ .- ϕ3) .+ sin(ϕ3))) + ξ0*ξ
		k2 = ω - a2.*(B2o*(sin.(B2'*(θ + h/2*k1) .- ϕ2) .+ sin(ϕ2))) - a3.*(B3o*(sin.(B3'*(θ + h/2*k1) .- ϕ3) .+ sin(ϕ3))) + ξ0*ξ
		k3 = ω - a2.*(B2o*(sin.(B2'*(θ + h/2*k2) .- ϕ2) .+ sin(ϕ2))) - a3.*(B3o*(sin.(B3'*(θ + h/2*k2) .- ϕ3) .+ sin(ϕ3))) + ξ0*ξ
		k4 = ω - a2.*(B2o*(sin.(B2'*(θ + h*k3) .- ϕ2) .+ sin(ϕ2))) - a3.*(B3o*(sin.(B3'*(θ + h*k3) .- ϕ3) .+ sin(ϕ3))) + ξ0*ξ

		dθ = (k1 + 2*k2 + 2*k3 + k4)/6
		θ += h*dθ
		θs = [θs θ]
		dθs = [dθs dθ]

		err = maximum(abs.(dθ))
	end

	Θs = zeros(n,0)
	dΘs = zeros(n,0)
	for i in 1:c
		Θs = [Θs readdlm("temp/ths_$i.csv",',')]
		rm("temp/ths_$i.csv")
		dΘs = [dΘs readdlm("temp/dths_$i.csv",',')]
		rm("temp/dths_$i.csv")
	end
	Θs = [Θs θs]
	dΘs = [dΘs dθs[:,2:end]]

	n,T = size(dΘs)

	Jhat = (diagm(0 => ones(n)) - dΘs*dΘs'/(T*ξ0^2))/τ0

	return Θs, dΘs, Jhat
end


# τ0: correlation time of the nodal noise
# ξ0: amplitude of the  noise
function hyper_k(A2::Array{Float64,2}, A3::Union{Array{Float64,2},Array{Float64,3}}, ω::Vector{Float64}, θ0::Vector{Float64}, τ0::Float64=1., ξ0::Float64=0., ϕ2::Float64=0., ϕ3::Float64=0., h::Float64=.01, max_iter::Int64=10000, tol::Float64=1e-6)
	n = length(ω)
	θs = θ0
	θ = θ0
	dθs = zeros(n,0)

	err = 1000.
	iter = 0
	c = 0

	ξ = rand(Normal(0.,1.),n)
	
	while iter < max_iter && err > tol
		iter += 1
		
		ξ = [cnoise(ξ[i],τ0/h) for i in 1:n]

		k1 = f_kuramoto_3rd(θ,A2,A3,ω,ϕ2,ϕ3) + ξ0*ξ
		k2 = f_kuramoto_3rd(θ+h/2*k1,A2,A3,ω,ϕ2,ϕ3) + ξ0*ξ
		k3 = f_kuramoto_3rd(θ+h/2*k2,A2,A3,ω,ϕ2,ϕ3) + ξ0*ξ
		k4 = f_kuramoto_3rd(θ+h*k3,A2,A3,ω,ϕ2,ϕ3) + ξ0*ξ

		dθ = (k1 + 2*k2 + 2*k3 + k4)/6
		θ += h*dθ

		θs = [θs θ]
		dθs = [dθs dθ]

		err = maximum(abs.(dθ))
		
		if iter%100 == 0
			@info "iter: $iter"
			c += 1
			writedlm("temp/ths_$c.csv",θs[:,1:end-1],',')
			θs = θs[:,end]
			writedlm("temp/dths_$c.csv",dθs[:,1:end],',')
			dθs = zeros(n,0)
		end
	end

	Θs = zeros(n,0)
	dΘs = zeros(n,0)
	for i in 1:c
		Θs = [Θs readdlm("temp/ths_$i.csv",',')]
		rm("temp/ths_$i.csv")
		dΘs = [dΘs readdlm("temp/dths_$i.csv",',')]
		rm("temp/dths_$i.csv")
	end
	Θs = [Θs θs[:,1:end-1]]
	dΘs = [dΘs dθs]

	return Θs, dΘs
end
	

function f_kuramoto_3rd(θ::Vector{Float64}, A2l::Array{Float64,2}, A3l::Array{Float64,2}, P::Vector{Float64}, ϕ2::Float64=0., ϕ3::Float64=0.)
	n = length(θ)
	
	fθ = copy(P)
	for l in 1:size(A2l)[1]
		i,j = Int64.(A2l[l,1:2])
		a = A2l[l,3]
		fθ[i] -= a*(sin(θ[i]-θ[j]-ϕ2) + sin(ϕ2))
	end
	for l in 1:size(A3l)[1]
		i,j,k = Int64.(A3l[l,1:3])
		a = A3l[l,4]
		fθ[i] -= a*(sin(2*θ[i]-θ[j]-θ[k]-ϕ3) + sin(ϕ3))
	end

	return fθ
end


function f_kuramoto_3rd(θ::Vector{Float64}, A2::Array{Float64,2}, A3::Array{Float64,3}, P::Vector{Float64}, ϕ2::Float64=0., ϕ3::Float64=0.)
	n = length(θ)
	
	fθ = Float64[]
	for i in 1:n
		x = P[i]
		for j in 1:n
			x -= A2[i,j]*(sin(θ[i]-θ[j]-ϕ2) + sin(ϕ2))
			for k in 1:n
				x -= A3[i,j,k]*(sin(2*θ[i]-θ[j]-θ[k] - ϕ3) + sin(ϕ3))
			end
		end
		push!(fθ,x)
	end
	return fθ
end

function f_kuramoto_3rd(Θ::Matrix{Float64}, A2::Array{Float64,2}, A3::Array{Float64,3}, P::Vector{Float64}, ϕ2::Float64=0., ϕ3::Float64=0.)
	n,T = size(Θ)
	fΘ = zeros(n,0)
	for t in 1:T
		fΘ = [fΘ f_kuramoto_3rd(Θ[:,t],A2,A3,P,ϕ2,ϕ3)]
	end
	return fΘ
end
