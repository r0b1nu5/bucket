using JuMP, Ipopt, PyPlot, Statistics

include("locate_forced.jl")

function test_w_forcing_l0(id::Int64, L::Array{Float64,2}, d::Array{Float64,1}, a0::Array{Float64,1}, f0::Array{Float64,1}, p0::Array{Float64,1}, Xs::Array{Float64,2}, dt::Float64, l::Float64=.01)
	nn,T = size(Xs)
	n = Int.(nn/2)

	system_id = Model(with_optimizer(Ipopt.Optimizer))
        
	@variable(system_id, Lm[i = 1:n, j = 1:n])
	@variable(system_id, dm[i = 1:n])
	@variable(system_id, a >= 0)
	@variable(system_id, f >= 0)
	@variable(system_id, 0 <= phi <= 2pi)
        
	for i in 1:n-1
                for j in i+1:n
                        @constraint(system_id, Lm[i,j] <= 0)
                        @constraint(system_id, Lm[j,i] <= 0)
                end
        end
        
        for i in 1:n
                @constraint(system_id, sum(Lm[i,:]) == 0)
        end

	@NLexpression(system_id, err[i = [1:id-1;id+1:n], t = 1:T-1], Xs[n+i,t+1] - Xs[n+i,t] - dt * (sum(-Lm[i,k]*Xs[k,t] for k = 1:n) - dm[i]*Xs[n+i,t]))
	@NLexpression(system_id, erri[t = 1:T-1], Xs[n+id,t+1] - Xs[n+id,t] - dt * (sum(-Lm[id,k]*Xs[k,t] for k = 1:n) - dm[id]*Xs[n+id,t] + a*cos(2*pi*f*dt*t + phi)))
#	@NLexpression(system_id, err[i = 1:n, t = 1:T-1], Xs[n+i,t+1] - Xs[n+i,t] - dt * (sum(-Lm[i,k]*Xs[k,t] for k = 1:n) - dm[i]*Xs[n+i,t]))

@NLobjective(system_id, Min, sum(err[i,t]^2 for i = [1:id-1;id+1:n] for t = 1:T-1) + sum(erri[t]^2 for t = 1:T-1)) 
        
        optimize!(system_id)

	reLm = abs.(L - value.(Lm))./abs.(L)
	redm = abs.(d - value.(dm))./abs.(d)
	rea = abs.(a0 .- value.(a))./abs.(a0)
	ref = abs.(f0 .- value.(f))./abs.(f0)
	rep = abs.(p0 .- value.(phi))./abs.(p0)

	aeLm = abs.(L - value.(Lm))
	aedm = abs.(d - value.(dm))
	aea = abs.(a0 .- value.(a))
	aef = abs.(f0 .- value.(f))
	aep = abs.(p0 .- value.(phi))
	
	figure()
	subplot(211)
	PyPlot.plot(1:n^2, vec(reLm), "o", label="M^{-1}*L")
	PyPlot.plot(n^2 .+ (1:n), redm, "o", label="d/m")
	PyPlot.plot(n^2+n .+ (1:n), rea, "o", label="a0")
	PyPlot.plot(n^2+2*n .+ (1:n), ref, "o", label="f0")
	PyPlot.plot(n^2+3*n .+ (1:n), rep, "o", label="phi0")
	title("Relative error")

	subplot(212)
	PyPlot.plot(1:n^2, vec(aeLm), "o", label="M^{-1}*L")
	PyPlot.plot(n^2 .+ (1:n), aedm, "o", label="d/m")
	PyPlot.plot(n^2+n .+ (1:n), aea, "o", label="a0")
	PyPlot.plot(n^2+2*n .+ (1:n), aef, "o", label="f0")
	PyPlot.plot(n^2+3*n .+ (1:n), aep, "o", label="phi0")
	title("Absolute error")
	legend()

        return (Lm = value.(Lm), dm = value.(dm), a = value.(a), f = value.(f), phi = value.(phi))
end

function test_w_forcing(L::Array{Float64,2}, d::Array{Float64,1}, a0::Array{Float64,1}, f0::Array{Float64,1}, p0::Array{Float64,1}, Xs::Array{Float64,2}, dt::Float64, l1::Float64=0., l2::Float64=0., mu::Float64=1e-1, bp::Float64=1e-1)
	nn,T = size(Xs)
	n = Int.(nn/2)

	system_id = Model(with_optimizer(Ipopt.Optimizer, mu_init = mu, bound_push = bp))
  	
	Ah,x = system_identification_correl(Xs,dt)
	L0 = -Ah[n+1:2*n,1:n]/dt

	@variable(system_id, Lm[i = 1:n, j = 1:n])
	set_start_value.(Lm,L0)
	@variable(system_id, dm[i = 1:n])
	@variable(system_id, a[i = 1:n])
	@variable(system_id, f[i = 1:n])
	@variable(system_id, phi[i = 1:n])
        
	#=
	for i in 1:n-1
                for j in i+1:n
                        @constraint(system_id, Lm[i,j] <= 0)
                        @constraint(system_id, Lm[j,i] <= 0)
                end
        end
        
        for i in 1:n
                @constraint(system_id, sum(Lm[i,:]) == 0)
        end
=#

	@NLexpression(system_id, err[i = 1:n, t = 1:T-1], Xs[n+i,t+1] - Xs[n+i,t] - dt * (sum(-Lm[i,k]*Xs[k,t] for k = 1:n) - dm[i]*Xs[n+i,t] + a[i]*cos(2*pi*f[i]*t*dt + phi[i])))
#	@NLexpression(system_id, err[i = 1:n, t = 1:T-1], Xs[n+i,t+1] - Xs[n+i,t] - dt * (sum(-Lm[i,k]*Xs[k,t] for k = 1:n) - dm[i]*Xs[n+i,t]))

@NLobjective(system_id, Min, sum(err[i,t]^2 for i = 1:n for t = 1:T-1) - l1 * sum(Lm[i,j] + Lm[j,i] for i = 1:n-1 for j = i+1:n) + l2 * sum(a[i] for i in 1:n)) 
        
        optimize!(system_id)
#=
	reLm = abs.(L - value.(Lm))./abs.(L)
	redm = abs.(d - value.(dm))./abs.(d)
	rea = abs.(a0 - value.(a))./abs.(a0)
	ref = abs.(f0 - value.(f))./abs.(f0)
	rep = abs.(p0 - value.(phi))./abs.(p0)

	aeLm = abs.(L - value.(Lm))
	aedm = abs.(d - value.(dm))
	aea = abs.(a0 - value.(a))
	aef = abs.(f0 - value.(f))
	aep = abs.(p0 - value.(phi))
	
	figure()
	subplot(211)
	PyPlot.plot(1:n^2, vec(reLm), "o", label="M^{-1}*L")
	PyPlot.plot(n^2 .+ (1:n), redm, "o", label="d/m")
	PyPlot.plot(n^2+n .+ (1:n), rea, "o", label="a0")
	PyPlot.plot(n^2+2*n .+ (1:n), ref, "o", label="f0")
	PyPlot.plot(n^2+3*n .+ (1:n), rep, "o", label="phi0")
	title("Relative error")

	subplot(212)
	PyPlot.plot(1:n^2, vec(aeLm), "o", label="M^{-1}*L")
	PyPlot.plot(n^2 .+ (1:n), aedm, "o", label="d/m")
	PyPlot.plot(n^2+n .+ (1:n), aea, "o", label="a0")
	PyPlot.plot(n^2+2*n .+ (1:n), aef, "o", label="f0")
	PyPlot.plot(n^2+3*n .+ (1:n), aep, "o", label="phi0")
	title("Absolute error")
	legend()
=#
        return (Lm = value.(Lm), dm = value.(dm), a = value.(a), f = value.(f), phi = value.(phi))
end

function test_w_forcing2(L::Array{Float64,2}, d::Array{Float64,1}, a0::Array{Float64,1}, f0::Array{Float64,1}, p0::Array{Float64,1}, Xs::Array{Float64,2}, Ablind::Array{Float64,2}, dt::Float64, mu::Float64=1e-1, bp::Float64=1e-1)
	nn,T = size(Xs)
	n = Int.(nn/2)
	
	Ah,x = system_identification_correl(Xs,dt)
	Li = -Ah[n+1:2*n,1:n]/dt

	df = rand(Normal(0,f0[1]),n,T)

	ids = Array{Tuple{Int64,Int64},1}()
	idx = Array{Int64,1}()
	li = zeros(n^2)
	for i in 1:n-1
		for j in i+1:n
			if Ablind[i,j] > 1e-8
				push!(ids,(i,j))
				push!(idx,Int(n*(i-1) + j))
				push!(ids,(j,i))
				push!(idx,Int(n*(j-1) + i))
				li[n*(i-1) + j] = Li[i,j]
				li[n*(j-1) + i] = Li[j,i]
			end
		end
	end

	system_id = Model(with_optimizer(Ipopt.Optimizer, mu_init = mu, bound_push = bp))
        
	@variable(system_id, Lm[i = idx] <= 0)
	for ii in idx
		set_start_value(Lm[ii],li[ii])
	end
	@variable(system_id, dm[i = 1:n])
	@variable(system_id, a[i = 1:n])
	@variable(system_id, f[i = 1:n])
	@variable(system_id, phi[i = 1:n])
        
	ks = Array{Array{Int64,1},1}()
	for i in 1:n
		push!(ks,Array{Int64,1}())
		for k in 1:n
			if length(intersect(ids,[(i,k),(k,i)])) > 0
				push!(ks[end],k)
			end
		end
	end
	@NLexpression(system_id, LX[i = 1:n, t = 1:T-1], sum(Lm[Int(n*(i-1) + k)]*(Xs[k,t] - Xs[i,t]) for k = ks[i]))
	
	@NLexpression(system_id, err[i = 1:n, t = 1:T-1], Xs[n+i,t+1] - Xs[n+i,t] - dt * (-LX[i,t] - dm[i]*Xs[n+i,t] + a[i]*cos(2*pi*(1+df[i,t])*f[i]*t*dt + phi[i])))

	@NLobjective(system_id, Min, 1/T*sum(err[i,t]^2 for i = 1:n for t = 1:T-1)) 
        
        optimize!(system_id)
	
	Lf = zeros(n,n)
	for i in 1:n-1
		for j in i+1:n
			if length(intersect(ids,[(i,j),(j,i)])) > 0
				Lf[i,j] = value(Lm[Int(n*(i-1) + j)])
				Lf[i,i] -= Lf[i,j]
				Lf[j,i] = value(Lm[Int(n*(j-1) + i)])
				Lf[j,j] -= Lf[j,i]
			end
		end
	end
#=
	reLm = abs.(L - Lf)./abs.(L)
	redm = abs.(d - value.(dm))./abs.(d)
	rea = abs.(a0 - value.(a))./abs.(a0)
	ref = abs.(f0 - value.(f))./abs.(f0)
	rep = abs.(p0 - value.(phi))./abs.(p0)

	aeLm = abs.(L - Lf)
	aedm = abs.(d - value.(dm))
	aea = abs.(a0 - value.(a))
	aef = abs.(f0 - value.(f))
	aep = abs.(p0 - value.(phi))
	
	figure()
	subplot(211)
	PyPlot.plot(1:n^2, vec(reLm), "o", label="M^{-1}*L")
	PyPlot.plot(n^2 .+ (1:n), redm, "o", label="d/m")
	PyPlot.plot(n^2+n .+ (1:n), rea, "o", label="a0")
	PyPlot.plot(n^2+2*n .+ (1:n), ref, "o", label="f0")
	PyPlot.plot(n^2+3*n .+ (1:n), rep, "o", label="phi0")
	title("Relative error")

	subplot(212)
	PyPlot.plot(1:n^2, vec(aeLm), "o", label="M^{-1}*L")
	PyPlot.plot(n^2 .+ (1:n), aedm, "o", label="d/m")
	PyPlot.plot(n^2+n .+ (1:n), aea, "o", label="a0")
	PyPlot.plot(n^2+2*n .+ (1:n), aef, "o", label="f0")
	PyPlot.plot(n^2+3*n .+ (1:n), aep, "o", label="phi0")
	title("Absolute error")
	legend()
=#
        return (Lm = Lf, dm = value.(dm), a = value.(a), f = value.(f), phi = value.(phi))
end

function find_freq_autocorr(Xs,dt,plot::Bool=false)
	nn,T = size(Xs)
	n = Int(nn/2)
	
	ac = Array{Float64,2}(undef,nn,0)
	leng = floor(Int,T/10)
	for i in 1:leng
		XX = Xs[:,1:end-i].*Xs[:,i+1:end]
		ac = [ac vec(sum(XX,dims=2)./(T-i))]
	end
	@info "Autocorrelations computed."

	if plot
		figure()
		for i in 1:n
			PyPlot.plot((1:floor(Int,T/10)).*dt,ac[i,:],"--",color="C$(i-1)")
			PyPlot.plot((1:floor(Int,T/10)).*dt,ac[n+i,:],color="C$(i-1)")
		end
	end

	fh = Array{Float64,1}()
	for i in 1:n
		@info "$i/$(2*n)"

		tops = Array{Int64,1}()
		bots = Array{Int64,1}()

		for j in 2:leng-1
			if (ac[i,j-1] > ac[i,j]) && (ac[i,j+1] > ac[i,j])
				push!(bots,j)
			elseif (ac[i,j-1] < ac[i,j]) && (ac[i,j+1] < ac[i,j])
				push!(tops,j)
			end
		end
		dtops = tops[2:end] - tops[1:end-1]
		dbots = bots[2:end] - bots[1:end-1]
		del = median([dtops;dbots])
		
		push!(fh,1/(del*dt))
	end

	return fh,ac
end


function test_w_forcing3(L::Array{Float64,2}, d::Array{Float64,1}, aa::Array{Float64,1}, ff::Array{Float64,1}, pp::Array{Float64,1}, Xs::Array{Float64,2}, dt::Float64, Dt::Int64, l1::Float64=0., l2::Float64=0., mu::Float64=1e-1, bp::Float64=1e-1)
	nn,T = size(Xs)
	n = Int.(nn/2)
	n_bin = floor(Int,T/Dt)

	system_id = Model(with_optimizer(Ipopt.Optimizer, mu_init = mu, bound_push = bp))
  	
	Ah,x = system_identification_correl(Xs,dt)
	L0 = -Ah[n+1:2*n,1:n]/dt
	
	f0 = mean(find_freq_autocorr(Xs,dt)[1])
	
	a0 = maximum(Xs[n+1:2*n,:])*sqrt(n)/f0

	@variable(system_id, Lm[i = 1:n, j = 1:n])
	set_start_value.(Lm,L0)
	@variable(system_id, dm[i = 1:n])
	@variable(system_id, a[i = 1:n])
	set_start_value.(a,a0*ones(n))
	@variable(system_id, f[i = 1:n])
	set_start_value.(f,f0*ones(n))
	@variable(system_id, phi[i = 1:n, j = 1:n_bin])
        
# #=
	for i in 1:n-1
                for j in i+1:n
                        @constraint(system_id, Lm[i,j] <= 1e-5)
                        @constraint(system_id, Lm[j,i] <= 1e-5)
                end
        end
# =#
# #=
        for i in 1:n
                @constraint(system_id, sum(Lm[i,:]) == 0)
        end
# =#

	@NLexpression(system_id, err[i = 1:n, t = 1:Dt-1, j = 1:n_bin], Xs[n+i,(j-1)*Dt + t + 1] - Xs[n+i,(j-1)*Dt + t] - dt * (sum(-Lm[i,k]*Xs[k,(j-1)*Dt + t] for k = 1:n) - dm[i]*Xs[n+i,(j-1)*Dt + t] + a[i]*cos(2*pi*f[i]*((j-1)*Dt + t)*dt + phi[i,j])))
#	@NLexpression(system_id, err[i = 1:n, t = 1:T-1], Xs[n+i,t+1] - Xs[n+i,t] - dt * (sum(-Lm[i,k]*Xs[k,t] for k = 1:n) - dm[i]*Xs[n+i,t]))

#@NLobjective(system_id, Min, 1/T*sum(err[i,t,j]^2 for i = 1:n for t = 1:Dt-1 for j = 1:n_bin) + l1 * sum(abs(Lm[i,j]) + abs(Lm[j,i]) for i = 1:n-1 for j = i+1:n) + l2 * sum(abs(a[i]) for i in 1:n)) 
@NLobjective(system_id, Min, 1/T*sum(err[i,t,j]^2 for i = 1:n for t = 1:Dt-1 for j = 1:n_bin)) 
        
        optimize!(system_id)
        
	return (Lm = value.(Lm), dm = value.(dm), a = value.(a), f = value.(f), phi = value.(phi))
end

function test_w_forcing3_l0(id::Int64, L::Array{Float64,2}, d::Array{Float64,1}, a0::Array{Float64,1}, f0::Array{Float64,1}, p0::Array{Float64,1}, Xs::Array{Float64,2}, dt::Float64, Dt::Int64, mu::Float64=1e-1, bp::Float64=1e-1)
	nn,T = size(Xs)
	n = Int.(nn/2)
	n_bin = floor(Int,T/Dt)

	system_id = Model(with_optimizer(Ipopt.Optimizer, mu_init = mu, bound_push = bp))
        
	Ah,x = system_identification_correl(Xs,dt)
	L0 = -Ah[n+1:2*n]/dt

	f0 = median(find_freq_autocorr(Xs,dt))

	a0 = maximum(abs.(Xs[n+1:2*n,:])*sqrt(n)/f0)

	@variable(system_id, Lm[i = 1:n, j = 1:n])
	set_start_value.(Lm,L0)
	@variable(system_id, dm[i = 1:n])
	@variable(system_id, a)
	set_start_value(a,a0)
	@variable(system_id, f)
	set_start_value(f,f0)
	@variable(system_id, phi[j = 1:n_bin])
        
	for i in 1:n-1
                for j in i+1:n
                        @constraint(system_id, Lm[i,j] <= 0)
                        @constraint(system_id, Lm[j,i] <= 0)
                end
        end
        
        for i in 1:n
                @constraint(system_id, sum(Lm[i,:]) == 0)
        end

	@NLexpression(system_id, err[i = [1:id-1;id+1:n], t = 1:Dt-1, j = 1:n_bin], Xs[n+i,(j-1)*Dt + t + 1] - Xs[n+i,(j-1)*Dt + t] - dt * (sum(-Lm[i,k]*Xs[k,(j-1)*Dt + t] for k = 1:n) - dm[i]*Xs[n+i,(j-1)*Dt + t]))
	@NLexpression(system_id, erri[t = 1:Dt-1, j = 1:n_bin], Xs[n+id,(j-1)*Dt + t + 1] - Xs[n+id,(j-1)*Dt + t] - dt * (sum(-Lm[id,k]*Xs[k,(j-1)*Dt + t] for k = 1:n) - dm[id]*Xs[n+id,(j-1)*Dt + t] + a*cos(2*pi*f*dt*((j-1)*Dt + t) + phi[j])))
#	@NLexpression(system_id, err[i = 1:n, t = 1:T-1], Xs[n+i,t+1] - Xs[n+i,t] - dt * (sum(-Lm[i,k]*Xs[k,t] for k = 1:n) - dm[i]*Xs[n+i,t]))

	@NLobjective(system_id, Min, 1/T*sum(err[i,t,j]^2 for i = [1:id-1;id+1:n] for t = 1:Dt-1 for j = 1:n_bin) + sum(erri[t,j]^2 for t = 1:Dt-1 for j = 1:n_bin)) 
        
        optimize!(system_id)

#=
	reLm = abs.(L - value.(Lm))./abs.(L)
	redm = abs.(d - value.(dm))./abs.(d)
	rea = abs.(a0 .- value.(a))./abs.(a0)
	ref = abs.(f0 .- value.(f))./abs.(f0)
	rep = abs.(p0 .- value.(phi))./abs.(p0)

	aeLm = abs.(L - value.(Lm))
	aedm = abs.(d - value.(dm))
	aea = abs.(a0 .- value.(a))
	aef = abs.(f0 .- value.(f))
	aep = abs.(p0 .- value.(phi))
	
	figure()
	subplot(211)
	PyPlot.plot(1:n^2, vec(reLm), "o", label="M^{-1}*L")
	PyPlot.plot(n^2 .+ (1:n), redm, "o", label="d/m")
	PyPlot.plot(n^2+n .+ (1:n), rea, "o", label="a0")
	PyPlot.plot(n^2+2*n .+ (1:n), ref, "o", label="f0")
	PyPlot.plot(n^2+3*n .+ (1:n), rep, "o", label="phi0")
	title("Relative error")

	subplot(212)
	PyPlot.plot(1:n^2, vec(aeLm), "o", label="M^{-1}*L")
	PyPlot.plot(n^2 .+ (1:n), aedm, "o", label="d/m")
	PyPlot.plot(n^2+n .+ (1:n), aea, "o", label="a0")
	PyPlot.plot(n^2+2*n .+ (1:n), aef, "o", label="f0")
	PyPlot.plot(n^2+3*n .+ (1:n), aep, "o", label="phi0")
	title("Absolute error")
	legend()
=#

        return (Lm = value.(Lm), dm = value.(dm), a = value.(a), f = value.(f), phi = value.(phi))
end
function test_w_forcing_warm_start(L::Array{Float64,2}, d::Array{Float64,1}, a0::Array{Float64,1}, f0::Array{Float64,1}, p0::Array{Float64,1}, Xs::Array{Float64,2}, dt::Float64, L0::Array{Float64,2}, d0::Array{Float64,1}, l1::Float64=0., l2::Float64=0., mu::Float64=1e-1)
	nn,T = size(Xs)
	n = Int.(nn/2)

	system_id = Model(with_optimizer(Ipopt.Optimizer, mu_init = mu, warm_start_init_point="yes", bound_push=1e-5))
 
	@variable(system_id, Lm[i = 1:n, j = 1:n])
	set_start_value.(Lm, L0)
	@variable(system_id, dm[i = 1:n])
	set_start_value.(dm, d0)
	@variable(system_id, a[i = 1:n] >= 0)
#	@variable(system_id, a[i = 1:n])
#	@constraint(system_id, 1000*a[i = 1:n] >= 0)
	set_start_value.(a,a0)
	@variable(system_id, f[i = 1:n] >= 0)
#	@variable(system_id, f[i = 1:n])
#	@constraint(system_id, 1000*f[i = 1:n] >= 0)
	set_start_value.(f,f0)
	@variable(system_id, 0 <= phi[i = 1:n] <= 2pi)
#	@variable(system_id, phi[i = 1:n])
	set_start_value.(phi,p0)

	e = zeros(n,T-1)
	for i in 1:n
		for t in 1:T-1
			e[i,t] = Xs[n+i,t+1] - Xs[n+i,t] - dt * (sum(-L0[i,k]*Xs[k,t] for k = 1:n) - d0[i]*Xs[n+i,t] + a0[i]*cos(2*pi*f0[i]*t*dt + p0[i]))
		end
	end

	ee = sum(e[i,t]^2 for i = 1:n for t = 1:T-1) - l1 * sum(L0[i,j] + L0[j,i] for i = 1:n-1 for j = i+1:n) + l2 * sum(a0[i] for i = 1:n)
	@info "error = $ee, max. $(maximum(abs.(e))), min. $(minimum(abs.(e)))"
	

	for i in 1:n-1
                for j in i+1:n
                        @constraint(system_id, Lm[i,j] <= 0)
                        @constraint(system_id, Lm[j,i] <= 0)
                end
        end
        
        for i in 1:n
                @constraint(system_id, sum(Lm[i,:]) == 0)
        end

	@NLexpression(system_id, err[i = 1:n, t = 1:T-1], (Xs[n+i,t+1] - Xs[n+i,t] - dt * (sum(-Lm[i,k]*Xs[k,t] for k = 1:n) - dm[i]*Xs[n+i,t] + a[i]*cos(2*pi*f[i]*t*dt + phi[i]))))
#	@NLexpression(system_id, err[i = 1:n, t = 1:T-1], Xs[n+i,t+1] - Xs[n+i,t] - dt * (sum(-Lm[i,k]*Xs[k,t] for k = 1:n) - dm[i]*Xs[n+i,t]))

@NLobjective(system_id, Min, sum(err[i,t]^2 for i = 1:n for t = 1:T-1) - l1 * sum(Lm[i,j] + Lm[j,i] for i = 1:n-1 for j = i+1:n) + l2 * sum(a[i] for i in 1:n)) 
        
        optimize!(system_id)

	reLm = abs.(L - value.(Lm))./abs.(L)
	redm = abs.(d - value.(dm))./abs.(d)
	rea = abs.(a0 - value.(a))./abs.(a0)
	ref = abs.(f0 - value.(f))./abs.(f0)
	rep = abs.(p0 - value.(phi))./abs.(p0)

	aeLm = abs.(L - value.(Lm))
	aedm = abs.(d - value.(dm))
	aea = abs.(a0 - value.(a))
	aef = abs.(f0 - value.(f))
	aep = abs.(p0 - value.(phi))
	
	figure()
	subplot(211)
	PyPlot.plot(1:n^2, vec(reLm), "o", label="M^{-1}*L")
	PyPlot.plot(n^2 .+ (1:n), redm, "o", label="d/m")
	PyPlot.plot(n^2+n .+ (1:n), rea, "o", label="a0")
	PyPlot.plot(n^2+2*n .+ (1:n), ref, "o", label="f0")
	PyPlot.plot(n^2+3*n .+ (1:n), rep, "o", label="phi0")
	title("Relative error")

	subplot(212)
	PyPlot.plot(1:n^2, vec(aeLm), "o", label="M^{-1}*L")
	PyPlot.plot(n^2 .+ (1:n), aedm, "o", label="d/m")
	PyPlot.plot(n^2+n .+ (1:n), aea, "o", label="a0")
	PyPlot.plot(n^2+2*n .+ (1:n), aef, "o", label="f0")
	PyPlot.plot(n^2+3*n .+ (1:n), aep, "o", label="phi0")
	title("Absolute error")
	legend()

        return (Lm = value.(Lm), dm = value.(dm), a = value.(a), f = value.(f), phi = value.(phi))
end




function test_wo_forcing(L::Array{Float64,2}, d::Array{Float64,1}, Xs::Array{Float64,2}, dt::Float64, l::Float64=.01)
	nn,T = size(Xs)
	n = Int.(nn/2)

	system_id = Model(with_optimizer(Ipopt.Optimizer))
        
	@variable(system_id, Lm[i = 1:n, j = 1:n])
	@variable(system_id, dm[i = 1:n])
#	@variable(system_id, a[i = 1:n] >= 0)
#	@variable(system_id, f[i = 1:n] >= 0)
#	@variable(system_id, 0 <= phi[i = 1:n] <= 2pi)
        
	for i in 1:n-1
                for j in i+1:n
                        @constraint(system_id, Lm[i,j] <= 0)
                        @constraint(system_id, Lm[j,i] <= 0)
                end
        end
        
        for i in 1:n
                @constraint(system_id, sum(Lm[i,:]) == 0)
        end

#	@NLexpression(system_id, err[i = 1:n, t = 1:T-1], Xs[n+i,t+1] - Xs[n+i,t] - dt * (sum(-Lm[i,k]*Xs[k,t] for k = 1:n) - dm[i]*Xs[n+i,t] + a[i]*cos(2*pi*f[i]*t*dt + phi[i])))
	@NLexpression(system_id, err[i = 1:n, t = 1:T-1], Xs[n+i,t+1] - Xs[n+i,t] - dt * (sum(-Lm[i,k]*Xs[k,t] for k = 1:n) - dm[i]*Xs[n+i,t]))

	@NLobjective(system_id, Min, sum(err[i,t]^2 for i = 1:n for t = 1:T-1) - l * sum(Lm[i,j] + Lm[j,i] for i = 1:n-1 for j = i+1:n)) 
        
        optimize!(system_id)

	reLm = abs.(L - value.(Lm))./abs.(L)
	redm = abs.(d - value.(dm))./abs.(d)
#	rea = abs.(a0 - value.(a))./abs.(a0)
#	ref = abs.(f0 - value.(f))./abs.(f0)
#	rep = abs.(p0 - value.(phi))./abs.(p0)

	aeLm = abs.(L - value.(Lm))
	aedm = abs.(d - value.(dm))
#	aea = abs.(a0 - value.(a))
#	aef = abs.(f0 - value.(f))
#	aep = abs.(p0 - value.(phi))
	
	figure()
	subplot(211)
	PyPlot.plot(1:n^2, vec(reLm), "o", label="M^{-1}*L")
	PyPlot.plot(n^2 .+ (1:n), redm, "o", label="d/m")
#	PyPlot.plot(n^2+n .+ (1:n), rea, "o", label="a0")
#	PyPlot.plot(n^2+2*n .+ (1:n), ref, "o", label="f0")
#	PyPlot.plot(n^2+3*n .+ (1:n), rep, "o", label="phi0")
	title("Relative error")

	subplot(212)
	PyPlot.plot(1:n^2, vec(aeLm), "o", label="M^{-1}*L")
	PyPlot.plot(n^2 .+ (1:n), aedm, "o", label="d/m")
#	PyPlot.plot(n^2+n .+ (1:n), aea, "o", label="a0")
#	PyPlot.plot(n^2+2*n .+ (1:n), aef, "o", label="f0")
#	PyPlot.plot(n^2+3*n .+ (1:n), aep, "o", label="phi0")
	title("Absolute error")
	legend()

        return (Lm = value.(Lm), dm = value.(dm))
end


function test_wo_forcing_lin(L::Array{Float64,2}, d::Array{Float64,1}, Xs::Array{Float64,2}, dt::Float64, l::Float64=.01)
	nn,T = size(Xs)
	n = Int.(nn/2)

	system_id = Model(with_optimizer(Ipopt.Optimizer))
        
	@variable(system_id, Lm[i = 1:n, j = 1:n])
	@variable(system_id, dm[i = 1:n])
#	@variable(system_id, a[i = 1:n] >= 0)
#	@variable(system_id, f[i = 1:n] >= 0)
#	@variable(system_id, 0 <= phi[i = 1:n] <= 2pi)
        
	for i in 1:n-1
                for j in i+1:n
                        @constraint(system_id, Lm[i,j] <= 0)
                        @constraint(system_id, Lm[j,i] <= 0)
                end
        end
        
        for i in 1:n
                @constraint(system_id, sum(Lm[i,:]) == 0)
        end

#	@NLexpression(system_id, err[i = 1:n, t = 1:T-1], Xs[n+i,t+1] - Xs[n+i,t] - dt * (sum(-Lm[i,k]*Xs[k,t] for k = 1:n) - dm[i]*Xs[n+i,t] + a[i]*cos(2*pi*f[i]*t*dt + phi[i])))
#	@expression(system_id, err[i = 1:n, t = 1:T-1], Xs[n+i,t+1] - Xs[n+i,t] - dt * (sum(-Lm[i,k]*Xs[k,t] for k = 1:n) - dm[i]*Xs[n+i,t]))

#	@objective(system_id, Min, sum(err[i,t]^2 for i = 1:n for t = 1:T-1) - l * sum(Lm[i,j] + Lm[j,i] for i = 1:n-1 for j = i+1:n)) 
@objective(system_id, Min, sum((Xs[n+i,t+1] - Xs[n+i,t] - dt * (sum(-Lm[i,k]*Xs[k,t] for k = 1:n) - dm[i]*Xs[n+i,t]))*(Xs[n+i,t+1] - Xs[n+i,t] - dt * (sum(-Lm[i,k]*Xs[k,t] for k = 1:n) - dm[i]*Xs[n+i,t])) for i = 1:n for t = 1:T-1) - l * sum(Lm[i,j] + Lm[j,i] for i = 1:n-1 for j = i+1:n))

        optimize!(system_id)

	reLm = abs.(L - value.(Lm))./abs.(L)
	redm = abs.(d - value.(dm))./abs.(d)
#	rea = abs.(a0 - value.(a))./abs.(a0)
#	ref = abs.(f0 - value.(f))./abs.(f0)
#	rep = abs.(p0 - value.(phi))./abs.(p0)

	aeLm = abs.(L - value.(Lm))
	aedm = abs.(d - value.(dm))
#	aea = abs.(a0 - value.(a))
#	aef = abs.(f0 - value.(f))
#	aep = abs.(p0 - value.(phi))
	
	figure()
	subplot(211)
	PyPlot.plot(1:n^2, vec(reLm), "o", label="M^{-1}*L")
	PyPlot.plot(n^2 .+ (1:n), redm, "o", label="d/m")
#	PyPlot.plot(n^2+n .+ (1:n), rea, "o", label="a0")
#	PyPlot.plot(n^2+2*n .+ (1:n), ref, "o", label="f0")
#	PyPlot.plot(n^2+3*n .+ (1:n), rep, "o", label="phi0")
	title("Relative error")

	subplot(212)
	PyPlot.plot(1:n^2, vec(aeLm), "o", label="M^{-1}*L")
	PyPlot.plot(n^2 .+ (1:n), aedm, "o", label="d/m")
#	PyPlot.plot(n^2+n .+ (1:n), aea, "o", label="a0")
#	PyPlot.plot(n^2+2*n .+ (1:n), aef, "o", label="f0")
#	PyPlot.plot(n^2+3*n .+ (1:n), aep, "o", label="phi0")
	title("Absolute error")
	legend()

        return (Lm = value.(Lm), dm = value.(dm))
end

