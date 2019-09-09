using PyPlot, LinearAlgebra, Dates

function kbc(w0::Array{Float64,1}, Dth::Float64, th0::Array{Float64,1}, eps::Float64=1e-8, max_iter::Int64=Int(1e5), h::Float64=.001)
	n = length(th0)
	
	B = inc_ata(n)

	ths = Array{Float64,2}(undef,n,0)
	ths = [ths th0]
	th = th0
	
	dth = zeros(n)

	err = 1000.
	c = 0

	while err > eps && c < max_iter
		c += 1
		if c%1000 == 0
			@info "$(now()) -- iter = $c, err = $err"
		end
		
		dth = mod.(B'*th .+ pi,2pi) .- pi
		test = diagm(0 => (abs.(dth) .<= Dth))
		Bt = B*test

		k1 = w0 - Bt*sin.(dth)./n
		k2 = w0 - Bt*sin.(dth + h/2*B'*k1)./n
		k3 = w0 - Bt*sin.(dth + h/2*B'*k2)./n
		k4 = w0 - Bt*sin.(dth + h*B'*k3)./n

		thd = (k1 + 2*k2 + 2*k3 + k4)/6

		th += h*thd
		ths = [ths th]

		err = maximum(abs.(h*thd))
	end

	return ths
end

function get_lap(th::Array{Float64,1}, Dth::Float64)
	n = length(th)

	Th = repeat(th,1,n)

	dth = Th - Th'

	A = (abs.(dth) .<= Dth)
	L = diagm(0 => vec(sum(A,dims=2))) - A

	return L
end

function inc_ata(n::Int64)
	if n == 2
		B = [1;-1]
	else
		B = [ones(1,n-1) zeros(1,Int((n-1)*(n-2)/2));-diagm(0 => ones(n-1)) inc_ata(n-1)]
	end

	return B
end



