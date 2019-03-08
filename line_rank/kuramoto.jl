include("L2B.jl")

function kuramoto(L::Array{Float64,2},P::Array{Float64,1},th0::Array{Float64,1},max_iter::Int=100000,eps::Float64=1e-6,h::Float64=0.1)
	n = size(L)[1]
	
	B,w = L2B(L)
	W = diagm(0 => w)
		
	th1 = zeros(n)
	th2 = copy(th0)
	
	error = 1000.
	
	iter = 0
	
	ths = Array{Float64,2}(undef,n,0)
	ths = [ths th0]
	
	while iter < max_iter && error > eps
		iter += 1
		if iter%1000 == 0 || iter == max_iter
			@info("iter = $iter, err = $error")
		end
		
		th1 = copy(th2)
		
		k1 = P - B*W*sin.(transpose(B)*th1)
		k2 = P - B*W*sin.(transpose(B)*(th1+(h/2)*k1))
		k3 = P - B*W*sin.(transpose(B)*(th1+(h/2)*k2))
		k4 = P - B*W*sin.(transpose(B)*(th1+h*k3))
		
		dth = (k1+2*k2+2*k3+k4)./6
		
		th2 = th1 + h*dth
		
		error = maximum(abs.(dth))
		
		ths = [ths th2]
	end
	
	return ths
end





