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




