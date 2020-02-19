using LinearAlgebra, Distributions, DelimitedFiles, SparseArrays

include("L2B.jl")

#TODO

function kuramoto_corr_noise(L::Array{Float64,2}, P::Array{Float64,1}, th0::Array{Float64,1}, dP0::Float64, tau0::Float64, store::Bool=false, max_iter::Int64=100000, eps::Float64=1e-8, h::Float64=.1)
	B,w = L2B(L)
	W = diagm(0 => w)
	Bt = transpose(B)

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

	while err > eps && iter < max_iter
		iter += 1
		if iter%1000 == 0
			@info "$iter"
		end

		xi = dP0^2*rand(Normal(0.,1.),n)

		th1 = copy(th2)
		
		k1 = P - B*W*sin.(Bt*th1) + xi
		k2 = P - B*W*sin.(Bt*(th1+h/2*k1)) + xi
		k3 = P - B*W*sin.(Bt*(th1+h/2*k2)) + xi
		k4 = P - B*W*sin.(Bt*(th1+h*k3)) + xi

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
	end

	return Ths
end
function kuramoto_white_noise(L::Array{Float64,2}, P::Array{Float64,1}, th0::Array{Float64,1}, dP0::Float64, store::Bool=false, max_iter::Int64=100000, eps::Float64=1e-8, h::Float64=.1)
	B,w = L2B(L)
	W = diagm(0 => w)
	Bt = transpose(B)

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

		xi = dP0^2*rand(Normal(0.,1.),n)

		th1 = copy(th2)
		
		k1 = P - B*W*sin.(Bt*th1) + xi
		k2 = P - B*W*sin.(Bt*(th1+h/2*k1)) + xi
		k3 = P - B*W*sin.(Bt*(th1+h/2*k2)) + xi
		k4 = P - B*W*sin.(Bt*(th1+h*k3)) + xi

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
	end

	return Ths
end

# max_iter: maximum number of iterations
# store: number of iterations to save (at the end of the time series)

function kuramoto_sine(L::Array{Float64,2}, P::Array{Float64,1}, th0::Array{Float64,1}, a0::Array{Float64,1}, w0::Array{Float64,1}, p0::Array{Float64,1}, store::Int64=1, max_iter::Int64=100000, eps::Float64=1e-8, h::Float64=.1)
	B,w = L2B(L)
	W = diagm(0 => w)
	Bt = transpose(B)

	n = length(th0)

	err = 1000.
	iter = 0

	th1 = copy(th0)
	th2 = copy(th0)
	if store > 1
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

		xi = a0.*cos.(w0.*h.*iter + p0)

		th1 = copy(th2)
		
		k1 = P - B*W*sin.(Bt*th1) + xi
		k2 = P - B*W*sin.(Bt*(th1+h/2*k1)) + xi
		k3 = P - B*W*sin.(Bt*(th1+h/2*k2)) + xi
		k4 = P - B*W*sin.(Bt*(th1+h*k3)) + xi

		dth = (k1+2*k2+2*k3+k4)/6

		th2 = th1 + h*dth

		if iter > (max_iter - store) && store > 1 && iter%1000 == 0
			c += 1
			writedlm("data1/ths_$c.csv",ths[:,1:end],',')
			ths = Array{Float64,2}(undef,n,0)
			ths = [ths th2]
		elseif store > 1 && iter > (max_iter - store)
			ths = [ths th2]
		else
			ths = copy(th2)
		end

		err = maximum(abs.(dth))
	end

	Ths = Array{Float64,2}(undef,n,0)
	for i in 1:c
		Ths = [Ths readdlm("data1/ths_$i.csv",',')]
	end

	return Ths
end

function kuramoto_sine(L::SparseMatrixCSC{Float64,Int64}, P::Array{Float64,1}, th0::Array{Float64,1}, a0::Array{Float64,1}, w0::Array{Float64,1}, p0::Array{Float64,1}, store::Int64=1, max_iter::Int64=100000, eps::Float64=1e-8, h::Float64=.1, file_name::String="data1")
	B,W,Bt = L2B(L)
	W = spdiagm(0 => W)

	n = length(th0)

	iter = store - max_iter

	th1 = copy(th0)
	th2 = copy(th0)
	if store > 1
		ths = Array{Float64,2}(undef,n,0)
	else
		ths = Array{Float64,1}()
	end
	
	c = 0

	while iter < 0
		iter += 1
		if iter%10000 == 0
			@info "$iter"
		end

		xi = a0.*sin.(w0.*h.*iter .+ p0)

		th1 = copy(th2)
		
		k1 = P - B*W*sin.(Bt*th1) + xi
		k2 = P - B*W*sin.(Bt*(th1+h/2*k1)) + xi
		k3 = P - B*W*sin.(Bt*(th1+h/2*k2)) + xi
		k4 = P - B*W*sin.(Bt*(th1+h*k3)) + xi

		dth = (k1+2*k2+2*k3+k4)/6

		th2 = th1 + h*dth
	end

	ths = [ths th2]
	
	while iter < store
		iter += 1
		if iter%1000 == 0
			@info "$iter"
		end
		

		xi = a0.*sin.(w0.*h.*iter + p0)

		th1 = copy(th2)

		k1 = P - B*W*sin.(Bt*th1) + xi
		k2 = P - B*W*sin.(Bt*(th1+h/2*k1)) + xi
		k3 = P - B*W*sin.(Bt*(th1+h/2*k2)) + xi
		k4 = P - B*W*sin.(Bt*(th1+h*k3)) + xi

		dth = (k1+2*k2+2*k3+k4)/6

		th2 = th1 + h*dth

		if store > 1 && iter%1000 == 0
			c += 1
			writedlm(file_name*"/ths_$c.csv",ths[:,1:end],',')
			ths = Array{Float64,2}(undef,n,0)
			ths = [ths th2]
		elseif store > 1
			ths = [ths th2]
		else
			ths = copy(th2)
		end
	end

	Ths = Array{Float64,2}(undef,n,0)
	for i in 1:c
		Ths = [Ths readdlm(file_name*"/ths_$i.csv",',')]
	end
	Ths = [Ths ths]

	return Ths
end


function kuramoto_step(L::Array{Float64,2}, P::Array{Float64,1}, th0::Array{Float64,1}, a0::Array{Float64,1}, T::Int64, store::Bool=false, max_iter::Int64=100000, eps::Float64=1e-8, h::Float64=.1)
	B,w = L2B(L)
	W = diagm(0 => w)
	Bt = transpose(B)

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

	while err > eps && iter < max_iter
		iter += 1

		if iter > T
			xi = a0
		else
			x = zeros(n)
		end

		th1 = copy(th2)
		
		k1 = P - B*W*sin.(Bt*th1) + xi
		k2 = P - B*W*sin.(Bt*(th1+h/2*k1)) + xi
		k3 = P - B*W*sin.(Bt*(th1+h/2*k2)) + xi
		k4 = P - B*W*sin.(Bt*(th1+h*k3)) + xi

		dth = (k1+2*k2+2*k3+k4)/6

		th2 = th1 + h*dth

		if store
			ths = [ths th2]
		else
			ths = copy(th2)
		end

		err = maximum(abs.(dth))
	end

	return ths
end




