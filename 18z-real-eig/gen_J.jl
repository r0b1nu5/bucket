using LinearAlgebra

function rand_J(n::Int64, α::Float64=.3, θmax::Float64=π/2)
	θ = θmax*rand(n)

	dθ = θ*ones(1,n) - ones(n)*θ'

	D = diagm(0 => -vec(sum(cos.(dθ .- α),dims=2)))
	cp = cos.(θ .+ α/2)
	cm = cos.(θ .- α/2)
	sp = sin.(θ .+ α/2)
	sm = sin.(θ .- α/2)
	c = cos.(θ)
	s = sin.(θ)

	J = D + cm*cp' + sm*sp'

	return J,D,cp,cm,sp,sm,c,s
end

function gen_J(θ::Array{Float64,1}, α::Float64=.3)
	n = length(θ)

	dθ = θ*ones(1,n) - ones(n)*θ'

	D = diagm(0 => -vec(sum(cos.(dθ .- α),dims=2)))
	cp = cos.(θ .+ α/2)
	cm = cos.(θ .- α/2)
	sp = sin.(θ .+ α/2)
	sm = sin.(θ .- α/2)
	c = cos.(θ)
	s = sin.(θ)

	J = D + cm*cp' + sm*sp'

	return J,D,cp,cm,sp,sm,c,s
end

