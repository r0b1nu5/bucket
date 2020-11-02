using SparseArrays,Dates, LinearAlgebra, PyPlot

include("L2B.jl")
include("../get_rgb.jl")

## INPUT
# L: graph Laplacian
# P: natural frequencies
# th0: initial angles
# store_history: if "true", then returns all angles' time evolution
# verb: displays info along the simulation
function kuramoto(L::Array{Float64,2},P::Array{Float64,1},th0::Array{Float64,1},store_history::Bool=false,verb::Bool=true,max_iter::Int=10000,eps::Float64=1e-6,h::Float64=0.1)
	n = size(L)[1]
	
	B,w = L2B(L)
	W = diagm(0 => w)
		
	th1 = zeros(n)
	th2 = copy(th0)
	
	error = 1000.
	
	iter = 0
	
	ths = Array{Float64,2}(undef,n,0)
	ths = [ths th0]
	dths = Array{Float64,2}(undef,n,0)
	dth = Array{Float64,1}()

	while iter < max_iter && error > eps
		iter += 1
		if (iter%1000 == 0 || iter == max_iter) && verb
			@info "$(now()) -- iter = $iter, err = $error"
		end
		
		th1 = copy(th2)
		
		k1 = P - B*W*sin.(transpose(B)*th1)
		k2 = P - B*W*sin.(transpose(B)*(th1+(h/2)*k1))
		k3 = P - B*W*sin.(transpose(B)*(th1+(h/2)*k2))
		k4 = P - B*W*sin.(transpose(B)*(th1+h*k3))
		
		dth = (k1+2*k2+2*k3+k4)./6
		
		th2 = th1 + h*dth
		
		error = maximum(abs.(dth))
		
		if store_history
			ths = [ths th2]
			dths = [dths dth]
		end
	end
	
	if !store_history
		ths = th2
		dths = dth
	end
	return ths,dths,iter
end



function kuramoto_long(L::Array{Float64,2},P::Array{Float64,1},th0::Array{Float64,1},verb::Bool=true,max_iter::Int=10000000,eps::Float64=1e-6,h::Float64=0.1)
	n = size(L)[1]
	
	B,w = L2B(L)
	W = diagm(0 => w)
		
	th1 = zeros(n)
	th2 = copy(th0)
	
	error = 1000.
	
	iter = 0
	
	ths = Array{Float64,2}(undef,n,0)
	ths = [ths th0]
	dths = Array{Float64,2}(undef,n,0)
	
	while iter < max_iter && error > eps
		iter += 1
		if iter%10000 == 0 || iter == max_iter
			@info "$(now()) -- iter = $iter, err = $error"
		end
		
		th1 = copy(th2)
		
		k1 = P - B*W*sin.(transpose(B)*th1)
		k2 = P - B*W*sin.(transpose(B)*(th1+(h/2)*k1))
		k3 = P - B*W*sin.(transpose(B)*(th1+(h/2)*k2))
		k4 = P - B*W*sin.(transpose(B)*(th1+h*k3))
		
		dth = (k1+2*k2+2*k3+k4)./6
		
		th2 = mod.(th1 + h*dth .+ pi,2pi) .- pi
		
		error = maximum(abs.(dth))
		
		if th1[1]*th2[1] < 0
			ths = [ths th2]
			dths = [dths dth]
		end
	end
	
	return ths,dths,iter
end


function kuramoto_poincaremap(L::Array{Float64,2}, P::Array{Float64,1}, th0::Array{Float64,1}, col::String="C0", verb::Bool=true, max_iter::Int=10000, eps::Float64=1e-6, h::Float64=0.1)
	n = size(L)[1]
	
	B,w = L2B(L)
	W = diagm(0 => w)
		
	th1 = zeros(n)
	th2 = copy(th0)
	
	error = 1000.
	
	iter = 0
	
	ths = Array{Float64,2}(undef,n,0)
	ths = [ths th0]
	dths = Array{Float64,2}(undef,n,0)
	dth = Array{Float64,1}()

	figure("poincare")
	PyPlot.plot([-pi/2,-pi/2],[-pi,pi],"--k")
	PyPlot.plot([pi/2,pi/2],[-pi,pi],"--k")
	PyPlot.plot([-pi,pi],[-pi/2,-pi/2],"--k")
	PyPlot.plot([-pi,pi],[pi/2,pi/2],"--k")

	while iter < max_iter && error > eps
		iter += 1
		if (iter%1000 == 0 || iter == max_iter) && verb
			@info "$(now()) -- iter = $iter, err = $error"
		end
		
		th1 = copy(th2)
		
		k1 = P - B*W*sin.(transpose(B)*th1)
		k2 = P - B*W*sin.(transpose(B)*(th1+(h/2)*k1))
		k3 = P - B*W*sin.(transpose(B)*(th1+(h/2)*k2))
		k4 = P - B*W*sin.(transpose(B)*(th1+h*k3))
		
		dth = (k1+2*k2+2*k3+k4)./6
		
		th2 = th1 + h*dth
		
		error = maximum(abs.(dth))

		if (mod(th2[3]+pi,2pi)-pi)*(mod(th1[3]+pi,2pi)-pi) < 0 && th2[3] > 0
			PyPlot.plot(mod(th2[1]+pi,2pi)-pi,mod(th2[2]+pi,2pi)-pi,".",color=get_rgb(Float64(iter),0.,Float64(max_iter)))
		end
		
	end
end

