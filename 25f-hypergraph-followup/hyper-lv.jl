using LinearAlgebra, DelimitedFiles, Distributions

function hyper_lv(A2::Array{Float64,2}, 
		 A3::Union{Array{Float64,2},Array{Float64,3}}, 
		 r::Vector{Float64}, 
		 l::Vector{Float64}, 
		 x0::Vector{Float64}, 
		 h::Float64=.01, 
		 max_iter::Int64=10000, 
		 tol::Float64=1e-6,
		 verb::Bool = false)

	n = length(r)
	xs = x0
	x = x0
	dxs = zeros(n,0)

	err = 1000.
	iter = 0
	c = 0

	while iter < max_iter && err > tol
		iter += 1
		
		k1 = f_lv_3rd(x,A2,A3,r,l)
		k2 = f_lv_3rd(x+h/2*k1,A2,A3,r,l)
		k3 = f_lv_3rd(x+h/2*k2,A2,A3,r,l)
		k4 = f_lv_3rd(x+h*k3,A2,A3,r,l)

		dx = (k1 + 2*k2 + 2*k3 + k4)/6
                x += h*dx

		xs = [xs x]
		dxs = [dxs dx]

		err = abs(maximum(dx)-minimum(dx))
		
		if iter%100 == 0
			if verb
				@info "iter: $iter"
			end
			c += 1
			writedlm("temp/xs_$c.csv",xs[:,1:end-1],',')
			xs = xs[:,end]
			writedlm("temp/dxs_$c.csv",dxs[:,1:end],',')
			dxs = zeros(n,0)
		end
	end
#	@info "Total iter: $iter"

	Xs = zeros(n,0)
	dXs = zeros(n,0)
	for i in 1:c
		Xs = [Xs readdlm("temp/xs_$i.csv",',')]
		rm("temp/xs_$i.csv")
		dXs = [dXs readdlm("temp/dxs_$i.csv",',')]
		rm("temp/dxs_$i.csv")
	end
	Xs = [Xs xs[:,1:end-1]]
	dXs = [dXs dxs]

	return Xs, dXs, iter
end
	
function hyper_lv_gaussian_noise(A2::Array{Float64,2}, 
				 A3::Union{Array{Float64,2},Array{Float64,3}}, 
				 r::Vector{Float64}, 
				 l::Vector{Float64}, 
				 x0::Vector{Float64}, 
				 ξ0::Float64,
				 δt::Float64, # Time between two random shocks
				 h::Float64=.01, 
				 max_iter::Int64=10000, 
				 tol::Float64=1e-6,
				 verb::Bool = false)

	n = length(r)
	xs = x0
	x = x0
	dxs = zeros(n,0)
        Δ = floor(Int64,δt/h)

	err = 1000.
	iter = 0
	c = 0

	while iter < max_iter && err > tol
		iter += 1
		
		k1 = f_lv_3rd(x,A2,A3,r,l)
		k2 = f_lv_3rd(x+h/2*k1,A2,A3,r,l)
		k3 = f_lv_3rd(x+h/2*k2,A2,A3,r,l)
		k4 = f_lv_3rd(x+h*k3,A2,A3,r,l)

                dx = (k1 + 2*k2 + 2*k3 + k4)/6 + ξ0*randn(n)*(mod(iter,Δ) == 0)
                x += h*dx
		
		xs = [xs x]
		dxs = [dxs dx]

		err = abs(maximum(dx)-minimum(dx))
		
		if iter%100 == 0
			if verb
				@info "iter: $iter"
			end
			c += 1
			writedlm("temp/xs_$c.csv",xs[:,1:end-1],',')
			xs = xs[:,end]
			writedlm("temp/dxs_$c.csv",dxs[:,1:end],',')
			dxs = zeros(n,0)
		end
	end

	Xs = zeros(n,0)
	dXs = zeros(n,0)
	for i in 1:c
		Xs = [Xs readdlm("temp/xs_$i.csv",',')]
		rm("temp/xs_$i.csv")
		dXs = [dXs readdlm("temp/dxs_$i.csv",',')]
		rm("temp/dxs_$i.csv")
	end
	Xs = [Xs xs[:,1:end-1]]
	dXs = [dXs dxs]

	return Xs, dXs, iter
end

function hyper_lv_drooped_gaussian_noise(A2::Array{Float64,2}, 
					 A3::Union{Array{Float64,2},Array{Float64,3}}, 
					 r::Vector{Float64}, 
					 l::Vector{Float64}, 
					 x0::Vector{Float64}, 
					 b::Vector{Float64},
					 xstar::Vector{Float64},
					 ξ0::Float64,
					 δt::Float64,
					 h::Float64=.01, 
					 max_iter::Int64=10000, 
					 tol::Float64=1e-6,
					 verb::Bool = false)

	n = length(r)
	xs = x0
	x = x0
	dxs = zeros(n,0)
        Δ = floor(Int64,δt/h)

	err = 1000.
	iter = 0
	c = 0

	while iter < max_iter && err > tol
		iter += 1
		
		k1 = f_lv_3rd_droop(x,A2,A3,r,l,b,xstar)
		k2 = f_lv_3rd_droop(x+h/2*k1,A2,A3,r,l,b,xstar)
		k3 = f_lv_3rd_droop(x+h/2*k2,A2,A3,r,l,b,xstar)
		k4 = f_lv_3rd_droop(x+h*k3,A2,A3,r,l,b,xstar)

		dx = (k1 + 2*k2 + 2*k3 + k4)/6 + ξ0*randn(n)*(mod(iter,Δ) == 0)

                x += h*dx

		xs = [xs x]
		dxs = [dxs dx]

		err = abs(maximum(dx)-minimum(dx))
		
		if iter%100 == 0
			if verb
				@info "iter: $iter"
			end
			c += 1
			writedlm("temp/xs_$c.csv",xs[:,1:end-1],',')
			xs = xs[:,end]
			writedlm("temp/dxs_$c.csv",dxs[:,1:end],',')
			dxs = zeros(n,0)
		end
	end
#	@info "Total iter: $iter"

	Xs = zeros(n,0)
	dXs = zeros(n,0)
	for i in 1:c
		Xs = [Xs readdlm("temp/xs_$i.csv",',')]
		rm("temp/xs_$i.csv")
		dXs = [dXs readdlm("temp/dxs_$i.csv",',')]
		rm("temp/dxs_$i.csv")
	end
	Xs = [Xs xs[:,1:end-1]]
	dXs = [dXs dxs]

	return Xs, dXs, iter
end

function hyper_lv_damped_gaussian_noise(A2::Array{Float64,2}, 
					A3::Union{Array{Float64,2},Array{Float64,3}}, 
					r::Vector{Float64},
					l::Vector{Float64},
					x0::Vector{Float64}, 
					d::Vector{Float64}, 
					ξ0::Float64,
					h::Float64=.01, 
					max_iter::Int64=10000, 
					tol::Float64=1e-6,
					verb::Bool = false)

	n = length(r)
	xs = x0
	x = x0
	dxs = zeros(n,0)

	err = 1000.
	iter = 0
	c = 0

	while iter < max_iter && err > tol
		iter += 1
		
		k1 = f_lv_3rd(x,A2,A3,r,l)
		k2 = f_lv_3rd(x+h/2*k1,A2,A3,r,l)
		k3 = f_lv_3rd(x+h/2*k2,A2,A3,r,l)
		k4 = f_lv_3rd(x+h*k3,A2,A3,r,l)

		dx = (k1 + 2*k2 + 2*k3 + k4)/6 + ξ0*randn(n)
		x += h*dx.*(1 .- d)


		xs = [xs x]
		dxs = [dxs dx]

		err = abs(maximum(dx)-minimum(dx))
		
		if iter%100 == 0
			if verb
				@info "iter: $iter"
			end
			c += 1
			writedlm("temp/xs_$c.csv",xs[:,1:end-1],',')
			xs = xs[:,end]
			writedlm("temp/dxs_$c.csv",dxs[:,1:end],',')
			dxs = zeros(n,0)
		end
	end
#	@info "Total iter: $iter"

	Xs = zeros(n,0)
	dXs = zeros(n,0)
	for i in 1:c
		Xs = [Xs readdlm("temp/xs_$i.csv",',')]
		rm("temp/xs_$i.csv")
		dXs = [dXs readdlm("temp/dxs_$i.csv",',')]
		rm("temp/dxs_$i.csv")
	end
	Xs = [Xs xs[:,1:end-1]]
	dXs = [dXs dxs]

	return Xs, dXs, iter
end

function f_lv_3rd(x::Vector{Float64}, A2l::Array{Float64,2}, A3l::Array{Float64,2}, r::Vector{Float64}, l::Vector{Float64})
	n = length(x)

	fx = r.*x.*(1 .- x./l)
	
	for l in 1:size(A2l)[1]
		i,j = Int64.(A2l[l,1:2])
		a = A2l[l,3]
		x[i] += a*x[i]*x[j]
	end
	for l in 1:size(A3l)[1]
		i,j,k = Int64.(A3l[l,1:3])
		a = A3l[l,4]
		fx[i] += a*x[i]*x[j]*x[k]
	end

	return fx
end


function f_lv_3rd(x::Vector{Float64}, A2::Array{Float64,2}, A3::Array{Float64,3}, r::Vector{Float64}, l::Vector{Float64})
	n = length(x)
	
	fx = Float64[]
	for i in 1:n
		y = r[i]*x[i]*(1 - x[i]/l[i])
		for j in 1:n
			y += A2[i,j]*x[i]*x[j]
			for k in 1:n
				y += A3[i,j,k]*x[i]*x[j]*x[k]
			end
		end
		push!(fx,y)
	end
	return fx
end

function f_lv_3rd(X::Matrix{Float64}, A2::Array{Float64,2}, A3::Array{Float64,3}, r::Vector{Float64}, l::Vector{Float64})
	n,T = size(X)
	fX = zeros(n,0)
	for t in 1:T
		fX = [fX f_lv_3rd(X[:,t],A2,A3,r,l)]
	end
	return fX
end


function f_lv_3rd_droop(x::Vector{Float64}, A2l::Array{Float64,2}, A3l::Array{Float64,2}, r::Vector{Float64}, l::Vector{Float64}, b::Vector{Float64}, xstar::Vector{Float64})
	n = length(x)
	
	fx = r.*x.*(1 .- r./l)
	for l in 1:size(A2l)[1]
		i,j = Int64.(A2l[l,1:2])
		a = A2l[l,3]
		fx[i] += a*x[i]*x[j]
	end
	for l in 1:size(A3l)[1]
		i,j,k = Int64.(A3l[l,1:3])
		a = A3l[l,4]
		fx[i] += a*x[i]*x[j]*x[k]
	end

	return fx - b.*(x - xstar)
end
