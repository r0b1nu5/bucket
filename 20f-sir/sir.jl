using Statistics, PyPlot

# R0: avg number of people infected by an infected (~2.7 [\pm 0.1])
# d: contagious time (~5-6 [2-14], as of March 17)


function sir(N::Int64, R0::Float64, d::Float64, sir0::Array{Float64,1}, fignum::Int64=0, eps::Float64=1e-8, h::Float64=.01)
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
		
		k1 = f_sir(sir1,R0,d)
		k2 = f_sir(sir1 + h/2*k1,R0,d)
		k3 = f_sir(sir1 + h/2*k2,R0,d)
		k4 = f_sir(sir1 + h*k3,R0,d)

		dsir = (k1 + 2*k2 + 2*k3 + k4)/6

		sir2 = sir1 + h*dsir

		sirs = [sirs sir2]
	end

	if fignum > 0
		figure(fignum)
		PyPlot.plot(h*(1:size(sirs)[2]), sirs[1,:], color="C0")
		PyPlot.plot(h*(1:size(sirs)[2]), sirs[2,:], color="C3")
		PyPlot.plot(h*(1:size(sirs)[2]), sirs[3,:], color="C2")
	end

	return sirs
end

function sir(N::Int64, R0s::Array{Float64,1}, d::Float64, sir0::Array{Float64,1}, fignum::Int64=0, eps::Float64=1e-8, h::Float64=.01)
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
		R0 = R0s[iter]
		sir1 = copy(sir2)
		
		k1 = f_sir(sir1,R0,d)
		k2 = f_sir(sir1 + h/2*k1,R0,d)
		k3 = f_sir(sir1 + h/2*k2,R0,d)
		k4 = f_sir(sir1 + h*k3,R0,d)

		dsir = (k1 + 2*k2 + 2*k3 + k4)/6

		sir2 = sir1 + h*dsir

		sirs = [sirs sir2]

#		@info "$(maximum(abs.(dsir))) > $eps"
#		@info "$iter < $(length(b))"
	end

	if fignum > 0
		figure(fignum)

		subplot(2,1,1)
#		PyPlot.plot(h*(1:size(sirs)[2]), sirs[1,:], color="C0")
		PyPlot.plot(h*(1:size(sirs)[2]), sirs[2,:], color="C3")
#		PyPlot.plot(h*(1:size(sirs)[2]), sirs[3,:], color="C2")

		subplot(2,1,2)
		PyPlot.plot(h*(1:size(sirs)[2]), bs[1:size(sirs)[2]], color="C1")
	end

	return sirs
end

function sir(N::Int64, R0s_data::Tuple{Float64,Float64,Int64,Int64,Int64}, d::Float64, sir0::Array{Float64,1}, start_day::Int64, fignum::Int64=0, eps::Float64=1e-8, h::Float64=.01)
	if abs(sum(sir0) - 1) > 1e-6
		@info "Sum of S+I+R is not 1."
	end

	R0i,R0f,ti,tt,tf = R0s_data
	R0s = R0s_gen(R0i,R0f,ti,tt,tf)

	sir1 = sir0
	sir2 = sir0

	sirs = Array{Float64,2}(undef,3,0)
	sirs = [sirs sir0]

	ii = sir0[2]
	double = [1,]

	dsir = 1000.
	iter = 0

	while maximum(abs.(dsir)) > eps && iter < length(R0s)
		iter += 1
		R0 = R0s[iter]
		sir1 = copy(sir2)
		
		k1 = f_sir(sir1,R0,d)
		k2 = f_sir(sir1 + h/2*k1,R0,d)
		k3 = f_sir(sir1 + h/2*k2,R0,d)
		k4 = f_sir(sir1 + h*k3,R0,d)

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

	reti = floor(Int64,start_day/h)
	sris = Array{Float64,2}(undef,3,0)
	sris = [sris sir0]
	ris1 = copy(sir0)
	ris2 = copy(sir0)

	while reti > 0
		reti -= 1
		R0 = R0i

		ris1 = copy(ris2)

		k1 = -f_sir(ris1,R0,d)
		k2 = -f_sir(ris1 + h/2*k1,R0,d)
		k3 = -f_sir(ris1 + h/2*k2,R0,d)
		k4 = -f_sir(ris1 + h*k3,R0,d)

		dris = (k1 + 2*k2 + 2*k3 + k4)/6

		ris2 = ris1 + h*dris

		sris = [ris2 sris]
	end

	if fignum > 0
		figure(fignum)

		@info "$(size(sris)[2])"
		@info "$(size(sirs)[2])"

		subplot(2,1,1)
#		PyPlot.plot(h*[ti,ti], [0.,1.], "--k")
#		PyPlot.plot(h*[ti+tt,ti+tt], [0.,1.], "--k")
#		for d in double
#			PyPlot.plot(h*[d,d], [0.,1.], ":k")
#		end
#		PyPlot.plot(h*(1:size(sirs)[2]), sirs[1,:], color="C0")
		PyPlot.plot([h*(1:size(sris)[2]-1); h*(1:size(sirs)[2]) .+ start_day], [sris[2,1:end-1];sirs[2,:]],label="R0 = $(R0s[end]), day 1 = $(ti/100), day d = $((ti+tt)/100.)")
#		PyPlot.plot(h*(1:size(sirs)[2]), sirs[3,:], color="C2")

		subplot(2,1,2)
		PyPlot.semilogy([h*(1:size(sris)[2]-1); h*(1:size(sirs)[2]) .+ start_day], [sris[2,1:end-1];sirs[2,:]])
#		PyPlot.plot(h*(1:size(sirs)[2]), bs[1:size(sirs)[2]])
		
#		subplot(3,1,3)
#		PyPlot.plot([h*(1:size(sris)[2]-1); h*(1:size(sirs)[2]) .+ start_day], [sris[2,1:end-1]+sris[3,1:end-1]; sirs[2,:]+sirs[3,:]])
	end

	return sirs, sris
end

function f_sir(sir::Array{Float64,1},R0::Float64,d::Float64)
	s,i,r = sir

	return [-R0/d*s*i,
		R0/d*s*i - i/d,
		i/d]
end

function R0s_gen(R0i::Float64, R0f::Float64, ti::Int64, tt::Int64, tf::Int64)
	return [R0i*ones(ti);
		(R0f-R0i)/tt*(ti+1:ti+tt) .+ R0i .- (R0f-R0i)*ti/tt;
		R0f*ones(tf)]
end

# First step of the time series is February 29th, 2020 (last day with 0 cases in CH).

function plot_data(data_set::String, fignum::Int64, N::Int64=8500000)
	I_ref = vec(readdlm(data_set,','))
	i_ref = I_ref./N

	figure(fignum)

	subplot(2,1,1)
	PyPlot.plot(0:length(i_ref)-1, i_ref, "ok")
	
	subplot(2,1,2)
	PyPlot.plot(0:length(i_ref)-1, i_ref, "ok")
end

function plot_bound(fignum::Int64=0,bound::Float64=.05)
	figure(fignum)

	subplot(2,1,1)
	PyPlot.plot([0,350],[bound,bound],":k")
	subplot(2,1,2)
	PyPlot.plot([0,350],[bound,bound],":k")
end







