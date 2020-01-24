using PyPlot

include("cnoise.jl")
include("L2B.jl")

function noise_inf(L::SparseMatrixCSC{Float64,Int64}, P::Array{Float64,1}, th0::Array{Float64,1}, dP0::Float64, tau0::Float64, measured_iter::Int64, max_iter::Int64, eps::Float64=1e-8, h::Float64=.1)
	B,w,Bt = L2B(L)
	W = spdiagm(0 => w)

	n = length(th0)

	iter = 0

	th1 = copy(th0)
	th2 = copy(th0)

	c = 0

	wiwj = zeros(n,n)

	xi = dP0*randn(n)

	while iter < max_iter - measured_iter
		iter += 1
		if iter%1000 == 0
			@info "$iter"
		end

		xi = [dP0*cnoise(xi[i],tau0) for i in 1:n]

		th1 = copy(th2)

		k1 = P - B*W*sin.(Bt*th1) + xi
		k2 = P - B*W*sin.(Bt*(th1+h/2*k1)) + xi
		k3 = P - B*W*sin.(Bt*(th1+h/2*k2)) + xi
		k4 = P - B*W*sin.(Bt*(th1+h*k3)) + xi

		dth = (k1+2*k2+2*k3+k4)/6
		th2 = th1 + h*dth
	end
	while iter < max_iter
		iter += 1
		if iter%1000 == 0
			@info "$iter"
		end

		xi = [dP0*cnoise(xi[i],tau0) for i in 1:n]

		th1 = copy(th2)

		k1 = P - B*W*sin.(Bt*th1) + xi
		k2 = P - B*W*sin.(Bt*(th1+h/2*k1)) + xi
		k3 = P - B*W*sin.(Bt*(th1+h/2*k2)) + xi
		k4 = P - B*W*sin.(Bt*(th1+h*k3)) + xi

		dth = (k1+2*k2+2*k3+k4)/6
		th2 = th1 + h*dth

		wiwj += dth*dth'
	end

#	wiwj ./= measured_iter*dP0^2
	
	return wiwj
#	return (-wiwj + diagm(0 => ones(n)))/tau0
end






		



