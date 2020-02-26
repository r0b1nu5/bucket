using PyPlot, DelimitedFiles, SparseArrays

include("cnoise.jl")
include("L2B.jl")

function noise_inf(L::SparseMatrixCSC{Float64,Int64}, P::Array{Float64,1}, th0::Array{Float64,1}, dP0::Float64, tau0::Float64, measured_iter::Int64, max_iter::Int64, d_meas::Int64, eps::Float64=1e-8, h::Float64=.1, history::Bool=false)
	B,w,Bt = L2B(L)
	W = spdiagm(0 => w)

	n = length(th0)

	iter = 0

	th1 = copy(th0)
	th2 = copy(th0)

	c = 0

	wiwj = zeros(n,n)
	
	dP = randn(n)
	xi = dP0*dP

	dth = zeros(n)

	while iter < max_iter - measured_iter
		iter += 1
		if iter%1000 == 0
			@info "$iter"
		end
		
		dP = [cnoise(dP[i],tau0/h) for i in 1:n]
		xi = dP0*dP
#		xi = dP0*randn(n)

		th1 = copy(th2)

		k1 = P - B*W*sin.(Bt*th1) + xi
		k2 = P - B*W*sin.(Bt*(th1+h/2*k1)) + xi
		k3 = P - B*W*sin.(Bt*(th1+h/2*k2)) + xi
		k4 = P - B*W*sin.(Bt*(th1+h*k3)) + xi

		dth = (k1+2*k2+2*k3+k4)/6
		th2 = th1 + h*dth
	end
	
	dth0 = copy(dth)
	ths = Array{Float64,2}(undef,n,0)
	dths = Array{Float64,2}(undef,n,0)
	c = 0
	cm = 0

	while iter < max_iter
		iter += 1
		if iter%1000 == 0
			@info "$iter"
		end
		
		dP = [cnoise(dP[i],tau0/h) for i in 1:n]
		xi = dP0*dP
#		xi = dP0*randn(n)

		th1 = copy(th2)

		k1 = P - B*W*sin.(Bt*th1) + xi
		k2 = P - B*W*sin.(Bt*(th1+h/2*k1)) + xi
		k3 = P - B*W*sin.(Bt*(th1+h/2*k2)) + xi
		k4 = P - B*W*sin.(Bt*(th1+h*k3)) + xi

		dth = (k1+2*k2+2*k3+k4)/6
		th2 = th1 + h*dth
		
		ddth = dth .- mean(dth)

		c += 1
		if c%d_meas == 0
			cm += 1
#			wiwj += dth0*dth'
			wiwj += ddth*ddth'
		end

		dth0 = copy(dth)
		if history
			ths = [ths th2]
			dths = [dths dth]
		end
	end

	return wiwj./cm, cm, ths, dths
end


function noise_inf_prerand(L::SparseMatrixCSC{Float64,Int64}, P::Array{Float64,1}, th0::Array{Float64,1}, dP0::Float64, dPs::Array{Float64,2}, eps::Float64=1e-8, h::Float64=.1,history::Bool=false)
	B,w,Bt = L2B(L)
	W = spdiagm(0 => w)

	n = length(th0)

	iter = 0

	th1 = copy(th0)
	th2 = copy(th0)

	c = 0

	wiwj = zeros(n,n)
	
	xis = dP0*dPs
	max_iter = size(dPs)[2]

	dth = zeros(n)

	while iter < max_iter
		iter += 1
		if iter%1000 == 0
			@info "$iter"
		end
		
		xi = xis[:,iter]
		th1 = copy(th2)

		k1 = P - B*W*sin.(Bt*th1) + xi
		k2 = P - B*W*sin.(Bt*(th1+h/2*k1)) + xi
		k3 = P - B*W*sin.(Bt*(th1+h/2*k2)) + xi
		k4 = P - B*W*sin.(Bt*(th1+h*k3)) + xi

		dth = (k1+2*k2+2*k3+k4)/6
		th2 = th1 + h*dth
	end
	
	return iter, th2, dth
end




		



