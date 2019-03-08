using LinearAlgebra,Statistics

#= 
Possible perturbations:
	box: ("box",[dP0,T0,id]), dP0: amplitude of the box, T0: duration of the box (number of iterations), id: index of the perturbed node
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
	ths = Array{Float64,2}(undef,n,0)
		
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
			if iter%10000 == 0 || iter == max_iter
				@info("iter = $iter")
			end
			
			dP = [cnoise(dP[i],tau0,dP0) for i in 1:n]
			th1 = copy(th2)
			ths = [ths th1]
			
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
	
	C1 /= iter
	C2 /= iter
	
	
	end
	
	if pert[1] == "box"
		iter = 0
		dP0 = pert[2][1]
		T0 = pert[2][2]
		id = Int(pert[2][3])
		dP = zeros(n)
		error = 1000.
		
		th1 = zeros(n)
		th2 = copy(th)
		
		C1 = 0.0
		C2 = 0.0
		
		while error > eps
			iter += 1
			if iter%10000 == 0
				@info("iter = $iter")
			end
			
			dP[id] = dP0*(iter <= T0)
			th1 = copy(th2)
			ths = [ths th1]
				
			k1 = P + dP - B*W*sin.(transpose(B)*th1)
			k2 = P + dP - B*W*sin.(transpose(B)*(th1+(h/2)*k1))
			k3 = P + dP - B*W*sin.(transpose(B)*(th1+(h/2)*k2))
			k4 = P + dP - B*W*sin.(transpose(B)*(th1+h*k3))
			
			thd = (k1+2*k2+2*k3+k4)./6
			
			th2 = th1 + h*thd
			
			dth = th2 - th
			C1 += sum((dth .- mean(dth)).^2)
			C2 += sum((thd .- mean(thd)).^2)
			
			error = maximum(abs.(thd))
			if error < eps
				@info("iter = $iter")
			end
		end
	end
	
	return C1,C2,ths
end


function cnoise(dP::Float64,tau0::Float64,dP0::Float64)
	f = exp(-1/tau0)
	g = randn()
	
	return f*dP + sqrt(1-f^2)*g
end

