using Statistics,PyPlot

function kuramoto_q(q::Int64,om::Array{Float64,1},K::Float64,th0::Array{Float64,1})
	n = length(th0)
	h = .001
		
	B = Array{Float64,2}(undef,n,0)
	for i in 1:n-1
		for j in i+1:n
			l = zeros(n)
			l[[i,j]] = [1.,-1.]
			B = [B l]
		end
	end
	Bt = transpose(B)
	
	omega = om .-= mean(om)
		
	th1 = copy(th0)
	th2 = copy(th0)
	ths = Array{Float64,2}(undef,n,0)
	ths = [ths th0]
	err = 1000.
	c = 0

	while err > 1e-6 && c < 1e6
		if c%1000 == 0
			@info "c = $c, err = $err"
		end
		c += 1
		
		th1 = copy(th2)
		
		k1 = omega - K*B*sin.(q*Bt*th1)
		k2 = omega - K*B*sin.(q*Bt*(th1+h/2*k1))
		k3 = omega - K*B*sin.(q*Bt*(th1+h/2*k2))
		k4 = omega - K*B*sin.(q*Bt*(th1+h*k3))
		
		dth = (k1+2*k2+2*k3+k4)/6
		
		th2 = th1 + h*dth
		ths = [ths th2]
		
		err = maximum(abs.(dth))
	end
	
	return ths
end
	

function plot_q(n::Int64,q::Int64)
	th0 = 2*pi/q*rand(n)
	omega = 4*rand(n)
	
	ths00 = mod.(kuramoto_q(q,omega,1.,th0),2*pi)
	ths01 = mod.(kuramoto_q(q,omega,1.,th0+2*pi/q*[i%2 for i in 1:n]),2*pi)
	ths10 = mod.(kuramoto_q(q,omega,1.,th0+2*pi/q*[i%2 for i in 0:n-1]),2*pi)
	ths11 = mod.(kuramoto_q(q,omega,1.,th0+2*pi/q*ones(n)),2*pi)
	pss = mod.(kuramoto_q(1,q*omega,Float64(q),q*th0),2*pi)
	
	for i in 1:round(Int,n/2)
		figure()
		PyPlot.plot(ths00[2*i-1,:],ths00[2*i,:])
		PyPlot.plot(ths00[2*i-1,1],ths00[2*i,1],"ok")
		PyPlot.plot(ths01[2*i-1,:],ths01[2*i,:])
		PyPlot.plot(ths10[2*i-1,:],ths10[2*i,:])
		PyPlot.plot(ths11[2*i-1,:],ths11[2*i,:])
		PyPlot.plot(pss[2*i-1,:],pss[2*i,:])
		PyPlot.plot(pss[2*i-1,:]/2,pss[2*i,:]/2,":")
	end
	
	return ths00,ths01,ths10,ths11,pss
end
