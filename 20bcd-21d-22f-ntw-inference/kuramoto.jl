using LinearAlgebra, Distributions, DelimitedFiles, SparseArrays

include("L2B.jl")

#TODO

function kuramoto(L::Array{Float64,2}, P::Array{Float64,1}, θ0::Array{Float64,1}, store::Bool=false, max_iter::Int64=100000, ε::Float64=1e-8, h::Float64=.01)
	B,w = L2B(L)
	W = diagm(0 => w)
	Bt = transpose(B)

	n = length(θ0)

	err = 1000.
	iter = 0

	θ1 = copy(θ0)
	θ2 = copy(θ0)
	if store
		θs = Array{Float64,2}(undef,n,0)
		θs = [θs θ0]
	else
		θs = Array{Float64,1}()
	end

	while err > ε && iter < max_iter
		iter += 1
		if iter%1000 == 0
			@info "$iter"
		end
		
		θ1 = copy(θ2)
		
		k1 = P - B*W*sin.(Bt*θ1)
		k2 = P - B*W*sin.(Bt*(θ1+h/2*k1))
		k3 = P - B*W*sin.(Bt*(θ1+h/2*k2))
		k4 = P - B*W*sin.(Bt*(θ1+h*k3))

		dθ = (k1+2*k2+2*k3+k4)/6

		θ2 = θ1 + h*dθ

		if store && iter%1000 == 0
			c += 1
			writedlm("data1/ths_$c.csv",θs[:,1:end],',')
			θs = Array{Float64,2}(undef,n,0)
			θs = [θs θ2]
		elseif store
			θs = [θs θ2]
		else
			θs = copy(θ2)
		end

		err = maximum(abs.(dθ))
	end

	Θs = Array{Float64,2}(undef,n,0)
	for i in 1:c
		Θs = [Θs readdlm("data1/ths_$i.csv",',')]
	end

	return Θs
end

function f_kuramoto(θ::Vector{Float64}, B::Matrix{Float64}, w::Vector{Float64}, P::Vector{Float64})
	return P - B*diagm(0 => w)*sin.(transpose(B)*θ)
end

function f_kuramoto(Θ::Matrix{Float64}, B::Matrix{Float64}, w::Vector{Float64}, P::Vector{Float64})
	n,T = size(Θ)
	fΘ = zeros(n,0)
	for t in 1:T
		fΘ = [fΘ f_kuramoto(Θ[:,t],B,w,P)]
	end
	return fΘ
end

function f_kuramoto_3rd(θ::Vector{Float64}, A2::Array{Float64,2}, A3::Array{Float64,3}, P::Vector{Float64})
	n = length(θ)
	fθ = Float64[]
	for i in 1:n
		x = P[i]
		for j in 1:n
			x -= A2[i,j]*sin(θ[i]-θ[j])
			for k in 1:n
				x -= A3[i,j,k]*sin(2*θ[i]-θ[j]-θ[k])
			end
		end
		push!(fθ,x)
	end
	return fθ
end

function f_kuramoto_3rd(Θ::Matrix{Float64}, A2::Array{Float64,2}, A3::Array{Float64,3}, P::Vector{Float64})
	n,T = size(Θ)
	fΘ = zeros(n,0)
	for t in 1:T
		fΘ = [fΘ f_kuramoto_3rd(Θ[:,t],A2,A3,P)]
	end
	return fΘ
end



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

	xi = dP0^2*rand(Normal(0.,1.),n)

	while err > eps && iter < max_iter
		iter += 1
		if iter%1000 == 0
			@info "$iter"
		end
		
		xi = dP0^2*cnoise(xi,tau0)

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

function kuramoto_white_noise(L::Array{Float64,2}, P::Array{Float64,1}, th0::Array{Float64,1}, dP0::Float64, store::Bool=false, max_iter::Int64=100000, eps::Float64=1e-8, h::Float64=.01)
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
		dths = Matrix{Float64}(undef,n,0)
	else
		ths = Array{Float64,1}()
		dths = Vector{Float64}()
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
		dths = [dths dth]

		if store && iter%1000 == 0
			c += 1
			writedlm("data1/ths_$c.csv",ths[:,1:end],',')
			ths = Array{Float64,2}(undef,n,0)
			ths = [ths th2]
			writedlm("data1/dths_$c.csv",dths[:,1:end],',')
			dths = Matrix{Float64}(undef,n,0)
		elseif store
			ths = [ths th2]
		else
			ths = copy(th2)
			dths = copy(dth)
		end

		err = maximum(abs.(dth))
	end

	if store
		Ths = Array{Float64,2}(undef,n,0)
		dThs = Matrix{Float64}(undef,n,0)
		for i in 1:c
			Ths = [Ths readdlm("data1/ths_$i.csv",',')]
			dThs = [dThs readdlm("data1/dths_$i.csv",',')]
		end
	else
		Ths = ths
		dThs = dths
	end

	return Ths,dThs
end

# max_iter: maximum number of iterations
# store: number of iterations to save (at the end of the time series)

function kuramoto_sine(L::Array{Float64,2}, P::Array{Float64,1}, th0::Array{Float64,1}, a0::Array{Float64,1}, w0::Array{Float64,1}, p0::Array{Float64,1}, store::Int64=1, max_iter::Int64=100000, eps::Float64=1e-8, h::Float64=.01)
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

		xi = a0.*sin.(w0.*h.*iter + p0)

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
	Ths = [Ths ths]

	return Ths
end

### NODAL DISTURBANCE!!! 

function kuramoto_sine(L::SparseMatrixCSC{Float64,Int64}, P::Vector{Float64}, th0::Vector{Float64}, i::Int64, a0::Float64, w0::Float64, p0::Float64, store::Int64=1, max_iter::Int64=100000, eps::Float64=1e-8, h::Float64=.1,file_name::String="data1")
	B,w,Bt = L2B(L)
	n,m = size(B)
	a = zeros(m)
	a[l] = a0
	w = zeros(m)
	w[l] = w0
	p = zeros(m)
	p[l] = p0

	return kuramoto_sine(L,P,th0,a,w,p,store,max_iter,eps,h,file_name)
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




