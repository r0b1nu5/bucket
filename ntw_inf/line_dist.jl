using LinearAlgebra

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

	return Ths
end


function kuramoto_step(L::Array{Float64,2}, P::Array{Float64,1}, th0::Array{Float64,1}, l::Int64, a0::Float64, T0::Float64=1., store::Bool=false, max_iter::Int64=10000, eps::Float64=1e-8, h::Float64=.01)
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

		xi = a0.*(iter*h > T0)

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


function linear_step(L::Array{Float64,2}, P::Array{Float64,1}, th0::Array{Float64,1}, ij::Tuple{Int64,Int64}, a0::Float64, T0::Float64=1., store::Bool=false, max_iter::Int64=10000, eps::Float64=1e-8, h::Float64=.01)
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

		xi = a0.*(iter*h > T0)

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



