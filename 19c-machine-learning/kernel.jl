using LinearAlgebra

# ========================== POLYNOMIAL KERNEL ==================================

# training_data: Tuple containing the array of inputs, xt, and the array of outpus, yt.
function polynomial_kernel_training(training_data::Tuple{Array{Float64,2},Array{Float64,2}},degree::Int64,gamma::Float64=.1)
	xt = training_data[1]
	yt = training_data[2]
	
	nx = size(xt)[1]
	ny = size(yt)[1]
	
	@info "Computing Gram matrix..."
	phi = polynomial_basis(xt,degree)
	
	K = transpose(phi)*phi
		
	@info "Optimizing coefficients..."
	c = yt*inv(K + gamma*diagm(0 => ones(size(K)[1])))
	
	return c,phi
end

# xs: inputs for the prediction, i.e., the number of columns of xs is the time horizon
# xt: training inputs

function polynomial_kernel_prediction(xs::Array{Float64,2},phit::Array{Float64,2},c::Array{Float64,2},degree::Int64)
	@info "Predicting future values..."
	
	T = size(xs)[2]
	
	phis = polynomial_basis(xs,degree)
	
	ys = c*transpose(phit)*phis
	
	return ys
end


function polynomial_basis(x::Array{Float64,1},degree::Int64)
	n = length(x)
	
	p = [1.,]
	d = 0
	while d < degree
		d += 1
		q = copy(p)
		for i in 1:n
			p = [p;x[i]*q]
		end
	end
	
	return p
end

function polynomial_basis(xs::Array{Float64,2},degree::Int64)
	n = size(xs)[1]
	
	p = ones(1,size(xs)[2])
	d = 0
	while d < degree
		d += 1
		q = copy(p)
		for i in 1:n
			p = [p;repeat(xs[[i,],:],size(p)[1]).*p]
		end
	end
	
	return p
end


# ============================= GAUSSIAN KERNEL =========================

# training_data: Tuple containing the array of inputs, xt, and the array of outpus, yt.
function gaussian_kernel_training(training_data::Tuple{Array{Float64,2},Array{Float64,2}},rho::Float64,gamma::Float64=.1)
	xt = training_data[1]
	yt = training_data[2]
	
	nx = size(xt)[1]
	ny = size(yt)[1]
	N = size(xt)[2]
	
	@info "Computing Gram matrix..."
	K = zeros(N,N)
	for i in 1:N-1
		K[i,i] = 1.
		for j in i:N
			K[i,j] = gaussian_k(xt[:,i],xt[:,j])
			K[j,i] = gaussian_k(xt[:,i],xt[:,j])
		end
	end
			
	@info "Optimizing coefficients..."
	c = yt*inv(K + gamma*diagm(0 => ones(size(K)[1])))
	
	return c
end


# xs: inputs for the prediction, i.e., the number of columns of xs is the time horizon
# xt: training inputs

function gaussian_kernel_prediction(xs::Array{Float64,2},xt::Array{Float64,2},c::Array{Float64,2},rho::Float64)
	@info "Predicting future values..."
	
	T = size(xs)[2]
	
	ys = Array{Float64,2}(undef,size(c)[1],0)
	for t in 1:T
		k = gaussian_k(xt,xs[:,t])
		y = c*k
		ys = [ys y]
	end
	
	return ys
end


function gaussian_k(x1::Array{Float64,1},x2::Array{Float64,1},rho::Float64=1.)
	return exp(-rho*sum((x1-x2).^2))
end
	
function gaussian_k(x1::Array{Float64,2},x2::Array{Float64,1},rho::Float64=1.)
	return exp.(-rho*(vec(sum((x1-repeat(x2,1,size(x1)[2])).^2,dims=1))))
end



