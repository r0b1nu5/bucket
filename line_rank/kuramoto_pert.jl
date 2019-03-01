using LinearAlgebra,Statistics

#= 
Possible perturbations:
	box: "box" TODO
	white noise: "noise_w" TODO
	colored noise: ("noise_c",[tau0,dP0])

L: Laplacian matrix of the graph, possibly weighted, diag terms are positive, off-diag terms are negative.
=#

function kuramoto_pert(L::Array{Float64,2},P::Array{Float64,1},th0::Array{Float64,1},pert,eps::Float64=1e-8,h::Float64=.1,max_iter::Int=100000)
	n = size(L)[1]
	
	B,w = L2B(L)
	W = diagm(0 => w)
	
	err = 1000.0
	
	th1 = zeros(n)
	th2 = copy(th0)
	
	while err > eps
		th1 = copy(th2)
		
		k1 = P - B*W*sin.(transpose(B)*th1)
		k2 = P - B*W*sin.(transpose(B)*(th1+(h/2)*k1))
		k3 = P - B*W*sin.(transpose(B)*(th1+(h/2)*k2))
		k4 = P - B*W*sin.(transpose(B)*(th1+h*k3))
		
		thd = (k1+2*k2+2*k3+k4)/6
		
		th2 = th1 + h*thd
		
		err = norm(th2-th1,Inf)
	end
	
	th = mod.(th2 .+ pi,2pi) .- pi
	
	if pert[1] == "noise_c"
		iter = 0
		tau0 = pert[2][1]
		dP0 = pert[2][2]
		dP = zeros(n)
		
		th1 = zeros(n)
		th2 = copy(th)
		
		C1 = 0.0
		C2 = 0.0
		
		while iter < max_iter
			iter += 1
			if iter%10000 == 0
				@info("iter = $iter")
			end
			
			dP = [cnoise(dP[i],tau0,dP0) for i in 1:n]
			th1 = copy(th2)
			
			k1 = P + dP - B*W*sin.(transpose(B)*th1)
			k2 = P + dP - B*W*sin.(transpose(B)*(th1+(h/2)*k1))
			k3 = P + dP - B*W*sin.(transpose(B)*(th1+(h/2)*k2))
			k4 = P + dP - B*W*sin.(transpose(B)*(th1+h*k3))
			
			thd = (k1+2*k2+2*k3+k4)./6
				
			th2 = th1 + h*thd
			
			dth = th2 - th
			C1 += sum((dth .- mean(dth)).^2)
			C2 += sum((thd .- mean(thd)).^2)
		end
	end
	
	C1 /= iter
	C2 /= iter
	
	return C1,C2
end


function cnoise(dP::Float64,tau0::Float64,dP0::Float64)
	f = exp(-1/tau0)
	g = randn()
	
	return f*dP + sqrt(1-f^2)*g
end

function L2B(L::Array{Float64,2})
	n = size(L)[1]
	
	B = Array{Float64,2}(undef,n,0)
	w = Array{Float64,1}(undef,0)
	for i in 1:n-1
		for j in i+1:n
			if L[i,j] != 0.0
				ed = zeros(n)
				ed[i] = 1.0
				ed[j] = -1.0
				B = [B ed]
				push!(w,-L[i,j])
			end
		end
	end
	
	return B,w
end
