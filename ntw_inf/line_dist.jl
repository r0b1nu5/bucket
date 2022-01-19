using LinearAlgebra, DelimitedFiles

include("L2B.jl")
include("cnoise.jl")

###################### KURAMOTO DYNAMICS ##############################

function kuramoto_sine(L::Array{Float64,2}, P::Array{Float64,1}, th0::Array{Float64,1}, l::Int64, a0::Float64, w0::Float64, p0::Float64, store::Bool=false, max_iter::Int64=10000, eps::Float64=1e-8, h::Float64=.01)
	B,w = L2B(L)
	W = diagm(0 => w)
	Bt = transpose(B)

	dW = zeros(size(W))
	dW[l,l] = 1.

	n = length(th0)

	err = 1000.
	iter = 0

	th1 = copy(th0)
	th2 = copy(th0)
	if store
		ths = Array{Float64,2}(undef,n,0)
		ths = [ths th0]
	else
		ths = Array{Float64,1}()
	end
	
	c = 0

	while err > eps && iter < max_iter
		iter += 1
		if iter%1000 == 0
			@info "$iter"
		end

		xi = a0.*sin.(w0.*h.*iter + p0)

		th1 = copy(th2)
		
		k1 = P - B*(W + xi*dW)*sin.(Bt*th1)
		k2 = P - B*(W + xi*dW)*sin.(Bt*(th1 + h/2*k1))
		k3 = P - B*(W + xi*dW)*sin.(Bt*(th1 + h/2*k2))
		k4 = P - B*(W + xi*dW)*sin.(Bt*(th1 + h*k3))

		dth = (k1+2*k2+2*k3+k4)/6

		th2 = th1 + h*dth

		if store && iter%1000 == 0
			c += 1
			writedlm("data1/ths_$c.csv",ths[:,1:end],',')
			ths = Array{Float64,2}(undef,n,0)
			ths = [ths th2]
		elseif store
			ths = [ths th2]
		else
			ths = copy(th2)
		end
		
		err = maximum(abs.(dth))
	end

	Ths = Array{Float64,2}(undef,n,0)
	for i in 1:c
		Ths = [Ths readdlm("data1/ths_$i.csv",',')]
		rm("data1/ths_$i.csv")
	end
	Ths = [Ths ths]

	return Ths,iter
end

function kuramoto_sine(ntw::String, L::SparseMatrixCSC{Float64,Int64}, P::Array{Float64,1}, th0::Array{Float64,1}, l::Int64, a0::Float64, w0::Float64, p0::Float64, store::Bool=false, max_iter::Int64=10000, eps::Float64=1e-8, h::Float64=.01)
	B,w,Bt = L2B(L)
	W = spdiagm(0 => w)
	n,m = size(B)

	dW = spzeros(m,m)
	dW[l,l] = 1.

	err = 1000.
	iter = 0

	th1 = copy(th0)
	th2 = copy(th0)
	if store
		ths = Array{Float64,2}(undef,n,0)
		ths = [ths th0]
	else
		ths = Array{Float64,1}()
	end
	
	c = 0

	while err > eps && iter < max_iter
		iter += 1
		if iter%1000 == 0
			@info ntw*"$iter"
		end

		xi = a0.*sin.(w0.*h.*iter + p0)

		th1 = copy(th2)
		
		k1 = P - B*(W + xi*dW)*sin.(Bt*th1)
		k2 = P - B*(W + xi*dW)*sin.(Bt*(th1 + h/2*k1))
		k3 = P - B*(W + xi*dW)*sin.(Bt*(th1 + h/2*k2))
		k4 = P - B*(W + xi*dW)*sin.(Bt*(th1 + h*k3))

		dth = (k1+2*k2+2*k3+k4)/6

		th2 = th1 + h*dth

		if store && iter%1000 == 0
			c += 1
			writedlm("data1/"*ntw*"_ths_$c.csv",ths[:,1:end],',')
			ths = Array{Float64,2}(undef,n,0)
			ths = [ths th2]
		elseif store
			ths = [ths th2]
		else
			ths = copy(th2)
		end
		
		err = maximum(abs.(dth))
	end

	Ths = Array{Float64,2}(undef,n,0)
	for i in 1:c
		Ths = [Ths readdlm("data1/"*ntw*"_ths_$i.csv",',')]
		rm("data1/"*ntw*"_ths_$i.csv")
	end
	Ths = [Ths ths]

	return Ths,iter
end

function kuramoto_step(L::Array{Float64,2}, P::Array{Float64,1}, th0::Array{Float64,1}, l::Int64, a0::Float64, T0::Float64, Tf::Float64, store::Bool=false, max_iter::Int64=10000, eps::Float64=1e-8, h::Float64=.01)
	B,w = L2B(L)
	W = diagm(0 => w)
	Bt = transpose(B)

	dW = zeros(size(W))
	dW[l,l] = 1.

	n = length(th0)

	err = 1000.
	iter = 0

	th1 = copy(th0)
	th2 = copy(th0)
	if store
		ths = Array{Float64,2}(undef,n,0)
		ths = [ths th0]
	else
		ths = Array{Float64,1}()
	end
	
	c = 0

	while err > eps && iter < max_iter
		iter += 1
		if iter%1000 == 0
			@info "$iter"
		end

		xi = a0.*(iter*h > T0).*(iter*h < Tf)

		th1 = copy(th2)
		
		k1 = P - B*(W + xi*dW)*sin.(Bt*th1)
		k2 = P - B*(W + xi*dW)*sin.(Bt*(th1 + h/2*k1))
		k3 = P - B*(W + xi*dW)*sin.(Bt*(th1 + h/2*k2))
		k4 = P - B*(W + xi*dW)*sin.(Bt*(th1 + h*k3))

		dth = (k1+2*k2+2*k3+k4)/6

		th2 = th1 + h*dth

		if store && iter%1000 == 0
			c += 1
			writedlm("data1/ths_$c.csv",ths[:,1:end],',')
			ths = Array{Float64,2}(undef,n,0)
			ths = [ths th2]
		elseif store
			ths = [ths th2]
		else
			ths = copy(th2)
		end
		
		err = maximum(abs.(dth))
	end

	Ths = Array{Float64,2}(undef,n,0)
	for i in 1:c
		Ths = [Ths readdlm("data1/ths_$i.csv",',')]
		rm("data1/ths_$i.csv")
	end
	Ths = [Ths ths]

	return Ths
end



function kuramoto_ramp(L::Array{Float64,2}, P::Array{Float64,1}, th0::Array{Float64,1}, l::Int64, a0::Float64, Ti::Float64=1., Tf::Float64=2., store::Bool=false, max_iter::Int64=10000, eps::Float64=1e-8, h::Float64=.01)
	B,w = L2B(L)
	W = diagm(0 => w)
	Bt = transpose(B)

	dW = zeros(size(W))
	dW[l,l] = 1.

	n = length(th0)

	err = 1000.
	iter = 0

	th1 = copy(th0)
	th2 = copy(th0)
	if store
		ths = Array{Float64,2}(undef,n,0)
		ths = [ths th0]
	else
		ths = Array{Float64,1}()
	end
	
	c = 0

	while err > eps && iter < max_iter
		iter += 1
		if iter%1000 == 0
			@info "$iter"
		end

		xi = a0*(iter*h - Ti)/(Tf - Ti)*(Ti < iter*h <= Tf) + a0*(iter*h > Tf)

		th1 = copy(th2)
		
		k1 = P - B*(W + xi*dW)*sin.(Bt*th1)
		k2 = P - B*(W + xi*dW)*sin.(Bt*(th1 + h/2*k1))
		k3 = P - B*(W + xi*dW)*sin.(Bt*(th1 + h/2*k2))
		k4 = P - B*(W + xi*dW)*sin.(Bt*(th1 + h*k3))

		dth = (k1+2*k2+2*k3+k4)/6

		th2 = th1 + h*dth

		if store && iter%1000 == 0
			c += 1
			writedlm("data1/ths_$c.csv",ths[:,1:end],',')
			ths = Array{Float64,2}(undef,n,0)
			ths = [ths th2]
		elseif store
			ths = [ths th2]
		else
			ths = copy(th2)
		end
		
		err = maximum(abs.(dth))
	end

	Ths = Array{Float64,2}(undef,n,0)
	for i in 1:c
		Ths = [Ths readdlm("data1/ths_$i.csv",',')]
		rm("data1/ths_$i.csv")
	end
	Ths = [Ths ths]

	return Ths

end

function kuramoto_noise(L::Array{Float64,2}, P::Array{Float64,1}, th0::Array{Float64,1}, l::Int64, a0::Float64, tau::Float64, store::Bool=false, max_iter::Int64=10000, eps::Float64=1e-8, h::Float64=.01)
	B,w = L2B(L)
	W = diagm(0 => w)
	Bt = transpose(B)

	dW = zeros(size(W))
	dW[l,l] = 1.

	n = length(th0)

	err = 1000.
	iter = 0

	th1 = copy(th0)
	th2 = copy(th0)
	if store
		ths = Array{Float64,2}(undef,n,0)
		ths = [ths th0]
	else
		ths = Array{Float64,1}()
	end
	
	c = 0
	xi = randn()

	while err > eps && iter < max_iter
		iter += 1
		if iter%1000 == 0
			@info "$iter"
		end

		xi = cnoise(xi,tau/h)

		th1 = copy(th2)
		
		k1 = P - B*(W + a0*xi*dW)*sin.(Bt*th1)
		k2 = P - B*(W + a0*xi*dW)*sin.(Bt*(th1 + h/2*k1))
		k3 = P - B*(W + a0*xi*dW)*sin.(Bt*(th1 + h/2*k2))
		k4 = P - B*(W + a0*xi*dW)*sin.(Bt*(th1 + h*k3))

		dth = (k1+2*k2+2*k3+k4)/6

		th2 = th1 + h*dth

		if store && iter%1000 == 0
			c += 1
			writedlm("data1/ths_$c.csv",ths[:,1:end],',')
			ths = Array{Float64,2}(undef,n,0)
			ths = [ths th2]
		elseif store
			ths = [ths th2]
		else
			ths = copy(th2)
		end
		
		err = maximum(abs.(dth))
	end

	Ths = Array{Float64,2}(undef,n,0)
	for i in 1:c
		Ths = [Ths readdlm("data1/ths_$i.csv",',')]
		rm("data1/ths_$i.csv")
	end
	Ths = [Ths ths]

	return Ths
end

###################### CONSENSUS DYNAMICS ##############################

function linear_sine(L::Array{Float64,2}, P::Array{Float64,1}, th0::Array{Float64,1}, ij::Tuple{Int64,Int64}, a0::Float64, w0::Float64, p0::Float64, store::Bool=false, max_iter::Int64=10000, eps::Float64=1e-8, h::Float64=.01)
	n = length(th0)
	l = [ij[1],ij[2]]

	dL = zeros(n,n)
	dL[l,l] = [1. -1.;-1. 1.]

	err = 1000.
	iter = 0

	th1 = copy(th0)
	th2 = copy(th0)
	if store
		ths = Array{Float64,2}(undef,n,0)
		ths = [ths th0]
	else
		ths = Array{Float64,1}()
	end
	
	c = 0

	while err > eps && iter < max_iter
		iter += 1
		if iter%1000 == 0
			@info "$iter"
		end

		xi = a0.*sin.(w0.*h.*iter + p0)

		th1 = copy(th2)
		
		k1 = P - (L + xi*dL)*th1
		k2 = P - (L + xi*dL)*(th1 + h/2*k1)
		k3 = P - (L + xi*dL)*(th1 + h/2*k2)
		k4 = P - (L + xi*dL)*(th1 + h*k3)

		dth = (k1+2*k2+2*k3+k4)/6

		th2 = th1 + h*dth

		if store && iter%1000 == 0
			c += 1
			writedlm("data1/ths_$c.csv",ths[:,1:end],',')
			ths = Array{Float64,2}(undef,n,0)
			ths = [ths th2]
		elseif store
			ths = [ths th2]
		else
			ths = copy(th2)
		end
		
		err = maximum(abs.(dth))
	end

	Ths = Array{Float64,2}(undef,n,0)
	for i in 1:c
		Ths = [Ths readdlm("data1/ths_$i.csv",',')]
		rm("data1/ths_$i.csv")
	end
	Ths = [Ths ths]

	return Ths
end


function linear_step(L::Array{Float64,2}, P::Array{Float64,1}, th0::Array{Float64,1}, ij::Tuple{Int64,Int64}, a0::Float64, Ti::Float64, Tf::Float64, store::Bool=false, max_iter::Int64=10000, eps::Float64=1e-8, h::Float64=.01)
	n = length(th0)
	l = [ij[1],ij[2]]

	dL = zeros(n,n)
	dL[l,l] = [1. -1.;-1. 1.]

	err = 1000.
	iter = 0

	th1 = copy(th0)
	th2 = copy(th0)
	if store
		ths = Array{Float64,2}(undef,n,0)
		ths = [ths th0]
	else
		ths = Array{Float64,1}()
	end
	
	c = 0

	while err > eps && iter < max_iter
		iter += 1
		if iter%1000 == 0
			@info "$iter"
		end

		xi = a0.*(iter*h > Ti).*(iter*h < Tf)

		th1 = copy(th2)
		
		k1 = P - (L + xi*dL)*th1
		k2 = P - (L + xi*dL)*(th1 + h/2*k1)
		k3 = P - (L + xi*dL)*(th1 + h/2*k2)
		k4 = P - (L + xi*dL)*(th1 + h*k3)

		dth = (k1+2*k2+2*k3+k4)/6

		th2 = th1 + h*dth

		if store && iter%1000 == 0
			c += 1
			writedlm("data1/ths_$c.csv",ths[:,1:end],',')
			ths = Array{Float64,2}(undef,n,0)
			ths = [ths th2]
		elseif store
			ths = [ths th2]
		else
			ths = copy(th2)
		end
		
		err = maximum(abs.(dth))
	end

	Ths = Array{Float64,2}(undef,n,0)
	for i in 1:c
		Ths = [Ths readdlm("data1/ths_$i.csv",',')]
		rm("data1/ths_$i.csv")
	end
	Ths = [Ths ths]

	return Ths
end



function linear_ramp(L::Array{Float64,2}, P::Array{Float64,1}, th0::Array{Float64,1}, ij::Tuple{Int64,Int64}, a0::Float64, Ti::Float64=1., Tf::Float64=2., store::Bool=false, max_iter::Int64=10000, eps::Float64=1e-8, h::Float64=.01)
	n = length(th0)
	l = [ij[1],ij[2]]

	dL = zeros(n,n)
	dL[l,l] = [1. -1.;-1. 1.]

	err = 1000.
	iter = 0

	th1 = copy(th0)
	th2 = copy(th0)
	if store
		ths = Array{Float64,2}(undef,n,0)
		ths = [ths th0]
	else
		ths = Array{Float64,1}()
	end
	
	c = 0

	while err > eps && iter < max_iter
		iter += 1
		if iter%1000 == 0
			@info "$iter"
		end

		xi = a0*(iter*h - Ti)/(Tf - Ti)*(Ti < iter*h <= Tf) + a0*(iter*h > Tf)

		th1 = copy(th2)
		
		k1 = P - (L + xi*dL)*th1
		k2 = P - (L + xi*dL)*(th1 + h/2*k1)
		k3 = P - (L + xi*dL)*(th1 + h/2*k2)
		k4 = P - (L + xi*dL)*(th1 + h*k3)

		dth = (k1+2*k2+2*k3+k4)/6

		th2 = th1 + h*dth

		if store && iter%1000 == 0
			c += 1
			writedlm("data1/ths_$c.csv",ths[:,1:end],',')
			ths = Array{Float64,2}(undef,n,0)
			ths = [ths th2]
		elseif store
			ths = [ths th2]
		else
			ths = copy(th2)
		end
		
		err = maximum(abs.(dth))
	end

	Ths = Array{Float64,2}(undef,n,0)
	for i in 1:c
		Ths = [Ths readdlm("data1/ths_$i.csv",',')]
		rm("data1/ths_$i.csv")
	end
	Ths = [Ths ths]

	return Ths

end

function linear_noise(L::Array{Float64,2}, P::Array{Float64,1}, th0::Array{Float64,1}, ij::Tuple{Int64,Int64}, a0::Float64, tau::Float64, store::Bool=false, max_iter::Int64=10000, eps::Float64=1e-8, h::Float64=.01)
	n = length(th0)
	l = [ij[1],ij[2]]

	dL = zeros(n,n)
	dL[l,l] = [1. -1.;-1. 1.]

	err = 1000.
	iter = 0

	th1 = copy(th0)
	th2 = copy(th0)
	if store
		ths = Array{Float64,2}(undef,n,0)
		ths = [ths th0]
	else
		ths = Array{Float64,1}()
	end
	
	c = 0
	xi = randn()

	while err > eps && iter < max_iter
		iter += 1
		if iter%1000 == 0
			@info "$iter"
		end

		xi = cnoise(xi,tau/h)

		th1 = copy(th2)
		
		k1 = P - (L + a0*xi*dL)*th1
		k2 = P - (L + a0*xi*dL)*(th1 + h/2*k1)
		k3 = P - (L + a0*xi*dL)*(th1 + h/2*k2)
		k4 = P - (L + a0*xi*dL)*(th1 + h*k3)

		dth = (k1+2*k2+2*k3+k4)/6

		th2 = th1 + h*dth

		if store && iter%1000 == 0
			c += 1
			writedlm("data1/ths_$c.csv",ths[:,1:end],',')
			ths = Array{Float64,2}(undef,n,0)
			ths = [ths th2]
		elseif store
			ths = [ths th2]
		else
			ths = copy(th2)
		end
		
		err = maximum(abs.(dth))
	end

	Ths = Array{Float64,2}(undef,n,0)
	for i in 1:c
		Ths = [Ths readdlm("data1/ths_$i.csv",',')]
		rm("data1/ths_$i.csv")
	end
	Ths = [Ths ths]

	return Ths
end


###################### OPINION DYNAMICS ##############################

function opinion_sine(L::Array{Float64,2}, P::Array{Float64,1}, th0::Array{Float64,1}, l::Int64, a0::Float64, w0::Float64, p0::Float64, store::Bool=false, max_iter::Int64=10000, eps::Float64=1e-8, h::Float64=.01)
	n = length(th0)
	
	B,S,T,w = L2B_dir(L)
	W = diagm(0 => w)
	Bt = Array(B')

	dW = zeros(size(W))
	dW[l,l] = 1.

	err = 1000.
	iter = 0

	th1 = copy(th0)
	th2 = copy(th0)
	if store
		ths = Array{Float64,2}(undef,n,0)
		ths = [ths th0]
	else
		ths = Array{Float64,1}()
	end
	
	c = 0

	while err > eps && iter < max_iter
		iter += 1
		if iter%1000 == 0
			@info "$iter"
		end

		xi = a0.*sin.(w0.*h.*iter + p0)

		th1 = copy(th2)
		
		k1 = P - S*(W + xi*dW)*tanh.(Bt*th1)
		k2 = P - S*(W + xi*dW)*tanh.(Bt*(th1 + h/2*k1))
		k3 = P - S*(W + xi*dW)*tanh.(Bt*(th1 + h/2*k2))
		k4 = P - S*(W + xi*dW)*tanh.(Bt*(th1 + h*k3))

		dth = (k1+2*k2+2*k3+k4)/6

		th2 = th1 + h*dth

		if store && iter%1000 == 0
			c += 1
			writedlm("data1/ths_$c.csv",ths[:,1:end],',')
			ths = Array{Float64,2}(undef,n,0)
			ths = [ths th2]
		elseif store
			ths = [ths th2]
		else
			ths = copy(th2)
		end
		
		err = maximum(abs.(dth))
	end

	Ths = Array{Float64,2}(undef,n,0)
	for i in 1:c
		Ths = [Ths readdlm("data1/ths_$i.csv",',')]
		rm("data1/ths_$i.csv")
	end
	Ths = [Ths ths]

	return Ths
end


function opinion_step(L::Array{Float64,2}, P::Array{Float64,1}, th0::Array{Float64,1}, l::Int64, a0::Float64, T0::Float64=1., store::Bool=false, max_iter::Int64=10000, eps::Float64=1e-8, h::Float64=.01)
	n = length(th0)
	
	B,S,T,w = L2B_dir(L)
	W = diagm(0 => w)
	Bt = transpose(B)

	dW = zeros(size(W))
	dW[l,l] = 1.
	
	err = 1000.
	iter = 0

	th1 = copy(th0)
	th2 = copy(th0)
	if store
		ths = Array{Float64,2}(undef,n,0)
		ths = [ths th0]
	else
		ths = Array{Float64,1}()
	end
	
	c = 0

	while err > eps && iter < max_iter
		iter += 1
		if iter%1000 == 0
			@info "$iter"
		end

		xi = a0.*(iter*h > T0)

		th1 = copy(th2)
		
		k1 = P - S*(W + xi*dW)*tanh.(Bt*th1)
		k2 = P - S*(W + xi*dW)*tanh.(Bt*(th1 + h/2*k1))
		k3 = P - S*(W + xi*dW)*tanh.(Bt*(th1 + h/2*k2))
		k4 = P - S*(W + xi*dW)*tanh.(Bt*(th1 + h*k3))

		dth = (k1+2*k2+2*k3+k4)/6

		th2 = th1 + h*dth

		if store && iter%1000 == 0
			c += 1
			writedlm("data1/ths_$c.csv",ths[:,1:end],',')
			ths = Array{Float64,2}(undef,n,0)
			ths = [ths th2]
		elseif store
			ths = [ths th2]
		else
			ths = copy(th2)
		end
		
		err = maximum(abs.(dth))
	end

	Ths = Array{Float64,2}(undef,n,0)
	for i in 1:c
		Ths = [Ths readdlm("data1/ths_$i.csv",',')]
		rm("data1/ths_$i.csv")
	end
	Ths = [Ths ths]

	return Ths
end


function opinion_ramp(L::Array{Float64,2}, P::Array{Float64,1}, th0::Array{Float64,1}, l::Int64, a0::Float64, Ti::Float64=1., Tf::Float64=2., store::Bool=false, max_iter::Int64=10000, eps::Float64=1e-8, h::Float64=.01)
	n = length(th0)
	B,S,T,w = L2B_dir(L)
	W = diagm(0 => w)
	Bt = transpose(B)

	dW = zeros(size(W))
	dW[l,l] = 1.
	
	err = 1000.
	iter = 0

	th1 = copy(th0)
	th2 = copy(th0)
	if store
		ths = Array{Float64,2}(undef,n,0)
		ths = [ths th0]
	else
		ths = Array{Float64,1}()
	end
	
	c = 0

	while err > eps && iter < max_iter
		iter += 1
		if iter%1000 == 0
			@info "$iter"
		end

		xi = a0*(iter*h - Ti)/(Tf - Ti)*(Ti < iter*h <= Tf) + a0*(iter*h > Tf)

		th1 = copy(th2)
		
		k1 = P - S*(W + xi*dW)*tanh.(Bt*th1)
		k2 = P - S*(W + xi*dW)*tanh.(Bt*(th1 + h/2*k1))
		k3 = P - S*(W + xi*dW)*tanh.(Bt*(th1 + h/2*k2))
		k4 = P - S*(W + xi*dW)*tanh.(Bt*(th1 + h*k3))

		dth = (k1+2*k2+2*k3+k4)/6

		th2 = th1 + h*dth

		if store && iter%1000 == 0
			c += 1
			writedlm("data1/ths_$c.csv",ths[:,1:end],',')
			ths = Array{Float64,2}(undef,n,0)
			ths = [ths th2]
		elseif store
			ths = [ths th2]
		else
			ths = copy(th2)
		end
		
		err = maximum(abs.(dth))
	end

	Ths = Array{Float64,2}(undef,n,0)
	for i in 1:c
		Ths = [Ths readdlm("data1/ths_$i.csv",',')]
		rm("data1/ths_$i.csv")
	end
	Ths = [Ths ths]

	return Ths

end

function opinion_noise(L::Array{Float64,2}, P::Array{Float64,1}, th0::Array{Float64,1}, l::Int64, a0::Float64, tau::Float64, store::Bool=false, max_iter::Int64=10000, eps::Float64=1e-8, h::Float64=.01)
	n = length(th0)
	B,S,T,w = L2B_dir(L)
	W = diagm(0 => w)
	Bt = transpose(B)

	dW = zeros(size(W))
	dW[l,l] = 1.
	
	err = 1000.
	iter = 0

	th1 = copy(th0)
	th2 = copy(th0)
	if store
		ths = Array{Float64,2}(undef,n,0)
		ths = [ths th0]
	else
		ths = Array{Float64,1}()
	end
	
	c = 0
	xi = randn()

	while err > eps && iter < max_iter
		iter += 1
		if iter%1000 == 0
			@info "$iter"
		end

		xi = cnoise(xi,tau/h)

		th1 = copy(th2)
		
		k1 = P - S*(W + a0*xi*dW)*tanh.(Bt*th1)
		k2 = P - S*(W + a0*xi*dW)*tanh.(Bt*(th1 + h/2*k1))
		k3 = P - S*(W + a0*xi*dW)*tanh.(Bt*(th1 + h/2*k2))
		k4 = P - S*(W + a0*xi*dW)*tanh.(Bt*(th1 + h*k3))

		dth = (k1+2*k2+2*k3+k4)/6

		th2 = th1 + h*dth

		if store && iter%1000 == 0
			c += 1
			writedlm("data1/ths_$c.csv",ths[:,1:end],',')
			ths = Array{Float64,2}(undef,n,0)
			ths = [ths th2]
		elseif store
			ths = [ths th2]
		else
			ths = copy(th2)
		end
		
		err = maximum(abs.(dth))
	end

	Ths = Array{Float64,2}(undef,n,0)
	for i in 1:c
		Ths = [Ths readdlm("data1/ths_$i.csv",',')]
		rm("data1/ths_$i.csv")
	end
	Ths = [Ths ths]

	return Ths
end

####################### HIGHER-ORDER KURAMOTO ####################

function kuramoto_q(L::Array{Float64,2}, qs::Array{Float64,1}, Ks::Array{Float64,1}, om::Array{Float64,1}, th0::Array{Float64,1}, ls::Array{Int64,1}, a0::Array{Float64,1}, Ti::Array{Float64,1}, Tf::Array{Float64,1}, sig::Float64=0., store::Bool=true, max_iter::Int64=10000, eps::Float64=1e-8, h::Float64=.01)
	n = length(th0)
	Q = length(qs)
	P = length(ls)

	B,w = L2B(L)
	m = length(w)
	W = diagm(0 => w)
	Bt = transpose(B)
	
	dWs = Array{Array{Float64,2},1}()
	for p in 1:P
		dWt = zeros(size(W))
		dWt[ls[p],ls[p]] = 1.
		push!(dWs,dWt)
	end
	
	omega = om .-= mean(om)
		
	th1 = copy(th0)
	th2 = copy(th0)
	if store
		ths = Array{Float64,2}(undef,n,0)
		ths = [ths th0]
		dhs = Array{Float64,2}(undef,n,0)
		dhs = [dhs zeros(n)]
		ds = Array{Float64,2}(undef,n,0)
		ds = [ds diag(L)]
	else
		ths = Array{Float64,1}
		dhs = Array{Float64,1}
		ds = Array{Float64,1}
	end
	err = 1000.
	c = 0
	iter = 0

	while err > eps && iter < max_iter
		if iter%1000 == 0
			@info "c = $c, err = $err"
		end
		iter += 1
		
		xis = [a0[i]*(iter*h - Ti[i])/(Tf[i] - Ti[i])*(Ti[i] < iter*h <= Tf[i]) + a0[i]*(iter*h > Tf[i]) for i in 1:P]
		dW = zeros(m,m)
		for p in 1:P
			dW += xis[p]*dWs[p]
		end

		Lt = B*(W-dW)*Bt

		eta = sig*randn(n)

		th1 = copy(th2)
		
		k1 = omega + eta - B*(W - dW)*sin.(Bt*th1*qs')*Ks
		k2 = omega + eta - B*(W - dW)*sin.(Bt*(th1+h/2*k1)*qs')*Ks
		k3 = omega + eta - B*(W - dW)*sin.(Bt*(th1+h/2*k2)*qs')*Ks
		k4 = omega + eta - B*(W - dW)*sin.(Bt*(th1+h*k3)*qs')*Ks
		
		dth = (k1+2*k2+2*k3+k4)/6
		
		th2 = th1 + h*dth
		err = maximum(abs.(dth))
		
		if store && iter%1000 == 0
			c += 1
			writedlm("data1/ths_$c.csv",ths[:,1:end],',')
			ths = Array{Float64,2}(undef,n,0)
			ths = [ths th2]
			writedlm("data1/dhs_$c.csv",dhs[:,1:end],',')
			dhs = Array{Float64,2}(undef,n,0)
			dhs = [dhs dth]
			writedlm("data1/ds_$c.csv",ds[:,1:end],',')
			ds = Array{Float64,2}(undef,n,0)
			ds = [ds diag(Lt)]
		elseif store
			ths = [ths th2]
			dhs = [dhs dth]
			ds = [ds diag(Lt)]
		else
			ths = copy(th2)
			dhs = copy(dth)
			ds = diag(Lt)
		end
	end

	Ths = Array{Float64,2}(undef,n,0)
	Dhs = Array{Float64,2}(undef,n,0)
	Ds = Array{Float64,2}(undef,n,0)
	for i in 1:c
		Ths = [Ths readdlm("data1/ths_$i.csv",',')]
		rm("data1/ths_$i.csv")
		Dhs = [Dhs readdlm("data1/dhs_$i.csv",',')]
		rm("data1/dhs_$i.csv")
		Ds = [Ds readdlm("data1/ds_$i.csv",',')]
		rm("data1/ds_$i.csv")
	end
	Ths = [Ths ths]
	Dhs = [Dhs dhs]
	Ds = [Ds ds]

	return Ths,Dhs,Ds
end


function kuramoto_q(L::SparseMatrixCSC{Float64,Int64}, qs::Array{Float64,1}, Ks::Array{Float64,1}, om::Array{Float64,1}, th0::Array{Float64,1}, ls::Array{Int64,1}, a0::Array{Float64,1}, Ti::Array{Float64,1}, Tf::Array{Float64,1}, sig::Float64=0., store::Bool=true, max_iter::Int64=10000, eps::Float64=1e-8, h::Float64=.01)
	n = length(th0)
	Q = length(qs)
	P = length(ls)

	B,w,Bt = L2B(L)
	m = length(w)
	W = spdiagm(0 => w)
	
	dWs = Array{SparseMatrixCSC{Float64,Int64},1}()
	for p in 1:P
		dWt = spzeros(length(w),length(w))
		dWt[ls[p],ls[p]] = 1.
		push!(dWs,dWt)
	end
	
	omega = om .-= mean(om)
		
	th1 = copy(th0)
	th2 = copy(th0)
	if store
		ths = Array{Float64,2}(undef,n,0)
		ths = [ths th0]
		dhs = Array{Float64,2}(undef,n,0)
		dhs = [dhs zeros(n)]
		ds = Array{Float64,2}(undef,n,0)
		ds = [ds diag(L)]
	else
		ths = Array{Float64,1}
		dhs = Array{Float64,1}
		ds = Array{Float64,1}
	end
	err = 1000.
	c = 0
	iter = 0

	io_t = open("data1/th.csv","a")
	io_dt = open("data1/dh.csv","a")
	io_prt = open("data1/prt.csv","a")

	while err > eps && iter < max_iter
		if iter%1000 == 0
			@info "iter = $iter, err = $err"
		end
		iter += 1
		
		xis = [a0[i]*(iter*h - Ti[i])/(Tf[i] - Ti[i])*(Ti[i] < iter*h <= Tf[i]) + a0[i]*(iter*h > Tf[i]) for i in 1:P]
		dW = spzeros(m,m)
		for p in 1:P
			dW += xis[p]*dWs[p]
		end

		Lt = B*(W-dW)*Bt

		eta = sig*randn(n)

		th1 = copy(th2)
		
		k1 = omega + eta - B*(W - dW)*sin.(Bt*th1*qs')*Ks
		k2 = omega + eta - B*(W - dW)*sin.(Bt*(th1+h/2*k1)*qs')*Ks
		k3 = omega + eta - B*(W - dW)*sin.(Bt*(th1+h/2*k2)*qs')*Ks
		k4 = omega + eta - B*(W - dW)*sin.(Bt*(th1+h*k3)*qs')*Ks
		
		dth = (k1+2*k2+2*k3+k4)/6
		
		th2 = th1 + h*dth
		err = maximum(abs.(dth))
		
		writedlm(io_t,th2',',')
		writedlm(io_dt,dth',',')
		writedlm(io_prt,xis',',')
		
#=
		if store && iter%100 == 0
			c += 1
			writedlm("data1/ths_$c.csv",ths[:,1:end],',')
			ths = Array{Float64,2}(undef,n,0)
			ths = [ths th2]
			writedlm("data1/dhs_$c.csv",dhs[:,1:end],',')
			dhs = Array{Float64,2}(undef,n,0)
			dhs = [dhs dth]
			writedlm("data1/ds_$c.csv",ds[:,1:end],',')
			ds = Array{Float64,2}(undef,n,0)
			ds = [ds diag(Lt)]
		elseif store
			ths = [ths th2]
			dhs = [dhs dth]
			ds = [ds diag(Lt)]
		else
			ths = copy(th2)
			dhs = copy(dth)
			ds = diag(Lt)
		end

=#
	end
#=
	Ths = Array{Float64,2}(undef,n,0)
	Dhs = Array{Float64,2}(undef,n,0)
	Ds = Array{Float64,2}(undef,n,0)
	for i in 1:c
		Ths = [Ths readdlm("data1/ths_$i.csv",',')]
		rm("data1/ths_$i.csv")
		Dhs = [Dhs readdlm("data1/dhs_$i.csv",',')]
		rm("data1/dhs_$i.csv")
		Ds = [Ds readdlm("data1/ds_$i.csv",',')]
		rm("data1/ds_$i.csv")
	end
	Ths = [Ths ths]
	Dhs = [Dhs dhs]
	Ds = [Ds ds]
=#
	close(io_t)
	close(io_dt)
	close(io_prt)

	Ths = Array(readdlm("data1/th.csv",',')')
	Dhs = Array(readdlm("data1/dh.csv",',')')
	Ps = Array(readdlm("data1/prt.csv",',')')

	rm("data1/th.csv")
	rm("data1/dh.csv")
	rm("data1/prt.csv")

	return Ths,Dhs,Ps
end


function linear_noise(L::SparseMatrixCSC{Float64,Int64}, om::Array{Float64,1}, x0::Array{Float64,1}, ls::Array{Int64,1}, a0::Array{Float64,1}, tau::Array{Float64,1}, sig::Float64=0., store::Bool=true, max_iter::Int64=10000, eps::Float64=1e-8, h::Float64=.01)
	n = length(th0)
	P = length(ls)

	B,w,Bt = L2B(L)
	m = length(w)
	W = spdiagm(0 => w)
	
	dLs = Array{SparseMatrixCSC{Float64,Int64},1}()
	for p in 1:P
		dWt = spzeros(length(w),length(w))
		dWt[ls[p],ls[p]] = 1.
		push!(dLs,B*dWt*Bt)
	end
	
	omega = om .-= mean(om)
		
	x1 = copy(x0)
	x2 = copy(x0)
	if store
		xs = Array{Float64,2}(undef,n,0)
		xs = [xs x0]
		dxs = Array{Float64,2}(undef,n,0)
		dxs = [dxs zeros(n)]
	else
		xs = Array{Float64,1}
		dxs = Array{Float64,1}
	end
	err = 1000.
	c = 0
	iter = 0

	xis = randn(P)

	io_x = open("data1/xs.csv","a")
	io_dx = open("data1/dxs.csv","a")
	io_prt = open("data1/prt.csv","a")

	while err > eps && iter < max_iter
		if iter%1000 == 0
			@info "iter = $iter, err = $err"
		end
		iter += 1
		
		xis = [cnoise(xis[p],tau[p]/h) for p in 1:P]
		dL = spzeros(n,n)
		for p in 1:P
			dL += a0[p]*xis[p]*dLs[p]
		end

		Lt = L + dL

		eta = sig*randn(n)

		x1 = copy(x2)
		
		k1 = omega + eta - (L + dL)*x1
		k2 = omega + eta - (L + dL)*(x1+h/2*k1)
		k3 = omega + eta - (L + dL)*(x1+h/2*k2)
		k4 = omega + eta - (L + dL)*(x1+h*k3)
		
		dx = (k1+2*k2+2*k3+k4)/6
		
		x2 = x1 + h*dx
		err = maximum(abs.(dx))

		writedlm(io_x,x2',',')
		writedlm(io_dx,dx',',')
		writedlm(io_prt,xis',',')

#=
		if store && iter%100 == 0
			c += 1
			writedlm("data1/xs_$c.csv",xs[:,1:end],',')
			xs = Array{Float64,2}(undef,n,0)
			xs = [xs x2]
			writedlm("data1/dxs_$c.csv",dxs[:,1:end],',')
			dxs = Array{Float64,2}(undef,n,0)
			dxs = [dxs dx]
		elseif store
			xs = [xs x2]
			dxs = [dxs dx]
		else
			xs = copy(x2)
			dxs = copy(dx)
		end
=#
end
#=
	Xs = Array{Float64,2}(undef,n,0)
	dXs = Array{Float64,2}(undef,n,0)
	for i in 1:c
		if i%10 == 0
			@info "Retrieve c = $i"
		end
		Xs = [Xs readdlm("data1/xs_$i.csv",',')]
		rm("data1/xs_$i.csv")
		dXs = [dXs readdlm("data1/dxs_$i.csv",',')]
		rm("data1/dxs_$i.csv")
	end
	Xs = [Xs xs]
	dXs = [dXs dxs]

=#
	close(io_x)
	close(io_dx)
	close(io_prt)

	Xs = Array(readdlm("data1/xs.csv",',')')
	dXs = Array(readdlm("data1/dxs.csv",',')')
	Ps = Array(readdlm("data1/prt.csv",',')')

	rm("data1/xs.csv")
	rm("data1/dxs.csv")
	rm("data1/prt.csv")

	return Xs,dXs,Ps
end



