using LinearAlgebra, PyPlot

include("tools.jl")

# DIRECTED VERSION OF THE KURAMOTO MODEL.

function kuramoto(L::Array{Float64,2}, om::Array{Float64,1}, th0::Array{Float64,1}, h::Float64=.01, eps::Float64=1e-6, max_iter::Int64=Int(1e6), doplot::Bool=true)
	n = length(th0)

	B,S,T,w = L2B_dir(L)
	W = diagm(0 => w)
	Bt = Array(B')

	err = 1000.
	iter = 0

	th1 = th0
	th2 = th0
	ths = Array{Float64,2}(undef,n,0)
	ths = [ths th0]
	dths = Array{Float64,2}(undef,n,0)

	while err > eps && iter < max_iter
		iter += 1
		if iter%100 == 0
			@info "iter: $iter, err: $err"
		end
		
		k1 = om - S*W*sin.(Bt*th1)
		k2 = om - S*W*sin.(Bt*(th1 + h/2*k1))
		k3 = om - S*W*sin.(Bt*(th1 + h/2*k2))
		k4 = om - S*W*sin.(Bt*(th1 + h*k3))

		dth = (k1 + 2*k2 + 2*k3 + k4)/6

		th2 = th1 + h*dth

		err = maximum(abs.(dth*ones(1,n) - ones(n)*dth'))

		th1 = copy(th2)
		ths = [ths th1]
		dths = [dths dth]
	end

	if doplot
		figure()
		for i in 1:n
			subplot(1,2,1)
			PyPlot.plot(ths[i,:])
			subplot(1,2,2)
			PyPlot.plot(dths[i,:])
		end
	end

	return ths,dths,iter
end






