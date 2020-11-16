using Statistics, LinearAlgebra, Dates

include("ksakaguchi.jl")

function winding(th::Array{Float64,1}, C::Array{Int64,1})
	th1 = th[C]
	th2 = th[[C[2:end];C[1]]]

	dth = th2 - th1

	q = sum(mod.(dth .+ pi,2pi) .- pi)/(2pi)

	return round(Int,q)
end

function jacobian(L::Array{Float64,2}, th::Array{Float64,1}, a::Float64)
	n = length(th)

	A = -L.*(1 .- diagm(0 => ones(n)))

	dth = th*ones(1,n) - ones(n)*th'

	J = A.*cos.(dth .- a)
	J = J - diagm(0 => vec(sum(J,dims=2)))

	return J
end

function freq_width(L::Array{Float64,2}, om0::Array{Float64,1}, th0::Array{Float64,1}, a::Float64, verb::Bool=false, res::Float64=.0005)
	if norm(om0) < 1e-8
		@info "The frequency vector is close to zero."
	end

	n = length(th0)
	
	om = om0 .- mean(om0)
	om /= norm(om)

	x = ksakaguchi(L,zeros(n),th0,a,true,false,.01,1e-6)
	th1 = x[1][:,end]
	th = copy(th1)
	q0 = winding(th,Array(1:n))

	b = 0.
	db = 1.

	while db > res
		if verb
			@info "b = $b, db = $db"
		end

		q = q0
		it = 0
		while q == q0 && it < 100000
			b += db
			x = ksakaguchi(L,b*om,th,a,true,false,.01,1e-6)
			q = winding(x[1][:,end],Array(1:n))
			it = x[4]
			
			if verb
				@info "q = $q, it = $it"
			end
		end
		th = x[1][:,1]
		b -= db
		db /= 10
	end

	bmax = copy(b)
	fmax = mean(x[2][:,end])

	b = 0.
	db = 1.

	while db > res
		if verb
			@info "b =$b, db = $db"
		end

		q = q0
		it = 0
		while q == q0 && it < 100000
			b -= db
			x = ksakaguchi(L,b*om,th,a,true,false,.01,1e-6)
			q = winding(x[1][:,end],Array(1:n))
			it = x[4]

			if verb
				@info "q = $q, it = $it"
			end
		end
		th = x[1][:,1]
		b += db
		db /= 10
	end

	bmin = copy(b)
	fmin = mean(x[2][:,end])

	return bmin,bmax,fmin,fmax
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




