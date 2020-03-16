using Statistics, PyPlot

function sir(N::Int64, b::Float64, g::Float64, sir0::Array{Float64,1}, plot::Bool=false, eps::Float64=1e-8, h::Float64=.01)
	if abs(sum(sir0) - 1) > 1e-6
		@info "Sum of S+I+R is not 1."
	end

	sir1 = sir0
	sir2 = sir0

	sirs = Array{Float64,2}(undef,3,0)
	sirs = [sirs sir0]

	dsir = 1000.

	while maximum(abs.(dsir)) > eps
		sir1 = copy(sir2)
		
		k1 = f_sir(sir1,b,g)
		k2 = f_sir(sir1 + h/2*k1,b,g)
		k3 = f_sir(sir1 + h/2*k2,b,g)
		k4 = f_sir(sir1 + h*k3,b,g)

		dsir = (k1 + 2*k2 + 2*k3 + k4)/6

		sir2 = sir1 + h*dsir

		sirs = [sirs sir2]
	end

	if plot
		PyPlot.plot(h*(1:size(sirs)[2]), sirs[1,:], color="C0")
		PyPlot.plot(h*(1:size(sirs)[2]), sirs[2,:], color="C3")
		PyPlot.plot(h*(1:size(sirs)[2]), sirs[3,:], color="C2")
	end

	return sirs
end

function sir(N::Int64, bs::Array{Float64,1}, g::Float64, sir0::Array{Float64,1}, plot::Bool=false, eps::Float64=1e-8, h::Float64=.01)
	if abs(sum(sir0) - 1) > 1e-6
		@info "Sum of S+I+R is not 1."
	end

	sir1 = sir0
	sir2 = sir0

	sirs = Array{Float64,2}(undef,3,0)
	sirs = [sirs sir0]

	dsir = 1000.
	iter = 0

	while maximum(abs.(dsir)) > eps && iter < length(bs)
		iter += 1
		b = bs[iter]
		sir1 = copy(sir2)
		
		k1 = f_sir(sir1,b,g)
		k2 = f_sir(sir1 + h/2*k1,b,g)
		k3 = f_sir(sir1 + h/2*k2,b,g)
		k4 = f_sir(sir1 + h*k3,b,g)

		dsir = (k1 + 2*k2 + 2*k3 + k4)/6

		sir2 = sir1 + h*dsir

		sirs = [sirs sir2]

#		@info "$(maximum(abs.(dsir))) > $eps"
#		@info "$iter < $(length(b))"
	end

	if plot
		subplot(2,1,1)
#		PyPlot.plot(h*(1:size(sirs)[2]), sirs[1,:], color="C0")
		PyPlot.plot(h*(1:size(sirs)[2]), sirs[2,:], color="C3")
#		PyPlot.plot(h*(1:size(sirs)[2]), sirs[3,:], color="C2")

		subplot(2,1,2)
		PyPlot.plot(h*(1:size(sirs)[2]), bs[1:size(sirs)[2]], color="C1")
	end

	return sirs
end

function sir(N::Int64, bs_data::Tuple{Float64,Float64,Int64,Int64,Int64}, g::Float64, sir0::Array{Float64,1}, plot::Bool=false, eps::Float64=1e-8, h::Float64=.01)
	if abs(sum(sir0) - 1) > 1e-6
		@info "Sum of S+I+R is not 1."
	end

	bi,bf,ti,tt,tf = bs_data
	bs = bs_gen(bi,bf,ti,tt,tf)

	sir1 = sir0
	sir2 = sir0

	sirs = Array{Float64,2}(undef,3,0)
	sirs = [sirs sir0]

	ii = sir0[2]
	double = [1,]

	dsir = 1000.
	iter = 0

	while maximum(abs.(dsir)) > eps && iter < length(bs)
		iter += 1
		b = bs[iter]
		sir1 = copy(sir2)
		
		k1 = f_sir(sir1,b,g)
		k2 = f_sir(sir1 + h/2*k1,b,g)
		k3 = f_sir(sir1 + h/2*k2,b,g)
		k4 = f_sir(sir1 + h*k3,b,g)

		dsir = (k1 + 2*k2 + 2*k3 + k4)/6

		sir2 = sir1 + h*dsir

		sirs = [sirs sir2]

		if sir2[2] > 2*ii
			push!(double,iter)
			ii = sir2[2]
		end

#		@info "$(maximum(abs.(dsir))) > $eps"
#		@info "$iter < $(length(b))"
	end

	if plot
		subplot(3,1,1)
#		PyPlot.plot(h*[ti,ti], [0.,1.], "--k")
#		PyPlot.plot(h*[ti+tt,ti+tt], [0.,1.], "--k")
#		for d in double
#			PyPlot.plot(h*[d,d], [0.,1.], ":k")
#		end
#		PyPlot.plot(h*(1:size(sirs)[2]), sirs[1,:], color="C0")
		PyPlot.plot(h*(1:size(sirs)[2]), sirs[2,:])
#		PyPlot.plot(h*(1:size(sirs)[2]), sirs[3,:], color="C2")

		subplot(3,1,2)
		PyPlot.semilogy(h*(1:size(sirs)[2]), sirs[2,:])
#		PyPlot.plot(h*(1:size(sirs)[2]), bs[1:size(sirs)[2]])
		
		subplot(3,1,3)
		PyPlot.plot(h*(1:size(sirs)[2]), sirs[2,:]+sirs[3,:])
	end

	return sirs
end

function f_sir(sir,b,g)
	s,i,r = sir

	return [-b*s*i,
		b*s*i - g*i,
		g*i]
end

function bs_gen(bi::Float64, bf::Float64, ti::Int64, tt::Int64, tf::Int64)
	return [bi*ones(ti);
		(bf-bi)/tt*(ti+1:ti+tt) .+ bi .- (bf-bi)*ti/tt;
		bf*ones(tf)]
end






