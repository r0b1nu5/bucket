using DelimitedFiles, PyPlot, FFTW, Statistics, LinearAlgebra, JuMP, Ipopt




################## TODO: run_location_large_ntw ######################################




# Identifies dynamics and forcing for a small network, base only on measurements.
#
# INPUT:
# Xs: Time series of the phase angles (rows 1:n) and of the phase frequencies (n+1:2*n).
# dt: Time step.
# Df: Radius of data over which the median of neighboring Fourier modes in computed. 
# n_period: Number of period considered in each data chunk.
# mu: Initial value of the barrier parameter (in IPOPT).
# bp: Initial value of the bound_push parameter (in IPOPT).
# Zro: Tolerance for zero.
#
# OUTPUT:
# Lm: Laplacian matrix normalized by the inertias.
# dm: Damping over inertia ratios.
# a: Vector of amplitudes of the forcing.
# f: Vector of frequencies of the forcing.
# p: Vector of phases of the forcing.

function run_location_small_ntw(Xs::Array{Float64,2}, dt::Float64, Df::Int64=10, n_period::Float64=1., mu::Float64=1e-5, bp::Float64=1e-5, Zro::Float64=1e-5)
	fh = get_fh_fourier(Xs,dt,Df)
	
	Ah,Lh,dh = get_Ah_correl(Xs,dt,fh)

	ah = get_ah(Xs,dt,fh)

	Lm,dm,a,f,ps = optim_chunks(Xs,dt,Lh,dh,ah,fh,n_period,mu,bp,Zro)

	p = optim_phase(Xs,dt,Lm,dm,a,f)

	return Lm,dm,a,f,p
end


# Estimates the forcings frequency, based on the Fourier Transform of the times series of the phase frequencies. A peak is identified by dividing the value of each Fourier mode by the median value of its 2*Df neighbors.
#
# INPUT:
# Xs: Time series of the phase angles (rows 1:n) and of the phase frequencies (n+1:2*n).
# dt: Time step.
# Df: Radius of data over which the median in computed. 
#
# OUTPUT:
# fh: Estimate of the forcing's frequency.

function get_fh_fourier(Xs::Array{Float64,2}, dt::Float64, Df::Int64=10)
	nn,T = size(Xs)
	n = Int(nn/2)

	fX = zeros(Complex{Float64},nn,T)
	for i in 1:n
		fX[i,:] = fft(Xs[n+i,:]).*dt./pi
	end

	nfX = norm.(fX)
	
	mfX = zeros(n,T-2*Df)
	for i in 1:n
		for j in Df+1:T-Df
			mfX[i,j-Df] = nfX[i,j]/median([nfX[i,j-Df:j-1];nfX[i,j+1:j+Df]])
		end
	end

	fs = (0:T-1)./(dt*T)
	df = fs[2]-fs[1]

	freqs = Array{Float64,1}()
	maxs = Array{Float64,1}()

	for i in 1:n
		ma,id = findmax(mfX[i,1:Int(T/2)])

		ma2,id2 = findmax([mfX[i,1:id-1];0.;mfX[i,id+1:Int(T/2)]])

		if (ma2 > .5*ma) && (abs(id - id2) == 1)
			push!(freqs,mean([fs[Df+id],fs[Df+id2]]))
		elseif (ma2 > .8*ma)
			push!(freqs,NaN)
			@info "Fourier Transform: inconclusive!"
		else
			push!(freqs,fs[Df+id])
		end
	end

	fh = median(freqs)
	tf = sum((abs.(freqs .- fh) .> 2*df))

	if tf/n > .05
		@info "Fourier Transform: no clear result."

		return NaN
	else
		return fh
	end
end


# Estimates the forcing's frequency based on the autocorrelation of the times series of the phase angles.
#
# INPUT:
# Xs: Time series of the phase angles (rows 1:n) and of the phase frequencies (n+1:2*n).
# dt: Time step.
# p: Proportion of the measure time over which the autocorrelation is computed.
#
# OUTPUT: 
# fh: Estimate of the forcing's frequency.

function get_fh_autocorr(Xs::Array{Float64,2}, dt::Float64, p::Float64=.1)
	nn,T = size(Xs)
	n = Int(nn/2)

	ac = Array{Float64,2}(undef,n,0)
	leng = floor(Int,p*T)
	for i in 1:leng
		XX = Xs[1:n,1:end-i].*Xs[1:n,i+1:end]
		ac = [ac vec(sum(XX,dims=2))./(T-i)]
	end
	@info "Autocorrelation computed."

	fh = Array{Float64,1}()
	for i in 1:n
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

	return median(fh)
end


# Estimates the dynamics matrix bases on Lokhov18.
#
# INPUT:
# Xs: Time series of the phase angles (rows 1:n) and of the phase frequencies (n+1:2*n).
# dt: Time step.
# fh: Estimate of forcing's frequency.
#
# OUTPUT: 
# Ah: Estimate of the full dynamics matrix.
# Lh: Estimate of the Laplacian matrix normalized by the inertias (M^{-1}*L).
# dh: Estimate of the damping over inertia ratios.

function get_Ah_correl(Xs::Array{Float64,2}, dt::Float64, fh::Float64)
	nn,T = size(Xs)
	n = Int(nn/2)

	Xt = filter_signal(Xs, dt, fh)

	S0 = Xt*Xt' ./ T
	S1 = Xt[:,2:T]*Xt[:,1:T-1]' ./ (T-1)

	Ah = S1*inv(S0)
	Lh = -Ah[n+1:2*n,1:n]/dt
	dh = (ones(n) - diag(Ah[n+1:2*n,n+1:2*n]))./dt

	return Ah,Lh,dh
end


# Removes the Fourier modes of frequency fh from the signal Xs.
#
# INPUT: 
# Xs: Time series of the phase angles (rows 1:n) and of the phase frequencies (n+1:2*n).
# dt: Time step.
# fh: Frequency to remove.
#
# OUPUT:
# Xt: Filtered signal.

function filter_signal(Xs::Array{Float64,2}, dt::Float64, fh::Float64)
	nn,T = size(Xs)

	fs = (0:T-1)./(dt*T)

	ta = ceil(Int64,T*5e-5)

	mi,id = findmin(abs.(fs .- fh))

	fX = zeros(Complex{Float64},nn,T)
	fXt = zeros(Complex{Float64},nn,T)
	Xt = zeros(nn,T)

	for i in 1:nn
		fX[i,:] = fft(Xs[i,:])
		fXt[i,:] = [fX[i,1:id-ta-1];zeros(2*ta+1);fX[i,id+ta+1:T-id+1-ta];zeros(2*ta+1);fX[i,T-id+3+ta:T]]
		Xt[i,:] = real.(ifft(fXt[i,:]))
	end

	return Xt
end


# Estimates the amplitude of the forcing. This is a very rough estimate, but we hope it is better than nothing.
#
# INPUT
# Xs: Time series of the phase angles (rows 1:n) and of the phase frequencies (n+1:2*n).
# dt: Time step.
# fh: Estimate of the forcing's frequency.
#
# OUTPUT:
# ah: Estimate of the forcing's amplitude.

function get_ah(Xs::Array{Float64,2}, dt::Float64, fh::Float64)
	nn,T = size(Xs)
	n = Int(nn/2)

	ah = maximum(abs.(Xs[1:n,:])) * sqrt(n) * (2*pi*fh)

	return ah
end


# Optimizes the dynamics matrix and forcing's paramters on chunks of the time series. The optimization allows different phase lags on each chunk, but requires the same amplitude and frequency.
#
# INPUT:
# Xs: Time series of the phase angles (rows 1:n) and of the phase frequencies (n+1:2*n).
# dt: Time step.
# Lh: Estimate of the normalized Laplacian matrix.
# dh: Estimate of the damping over inertia ratios.
# ah: Estimate of the forcing's amplitude.
# fh: Estimate of the forcing's frequency.
# n_period: Number of period considered in each data chunk.
# mu: Initial value of the barrier parameter (in IPOPT).
# bp: Initial value of the bound_push parameter (in IPOPT).
# Zro: Tolerance for zero.
#
# OUTPUT:
# Lm: Optimized normalized Laplacian matrix.
# dm: Optimized damping of inertia ratios.
# a: Optimized forcing amplitudes.
# f: Optimized forcing frequencies.
# p: Optimized phase lags in each chunk.

function optim_chunks(Xs::Array{Float64,2}, dt::Float64, Lh::Array{Float64,2}, dh::Array{Float64,1}, ah::Float64, fh::Float64, n_period::Float64, mu::Float64=1e-1, bp::Float64=1e-1, Zro::Float64=1e-5)
	nn,T = size(Xs)
	n = Int(nn/2)

	Dt = ceil(Int64,1/(fh*dt))
	n_bin = floor(Int64,T/Dt)

	system_id = Model(with_optimizer(Ipopt.Optimizer, mu_init = mu, bound_push = bp))

	@variable(system_id, Lm[i = 1:n, j = 1:n])
	set_start_value.(Lm,Lh)
	@variable(system_id, dm[i = 1:n])
	set_start_value.(dm,dh)
	@variable(system_id, a[i = 1:n])	
	set_start_value.(a,ah*ones(n))
	@variable(system_id, f[i = 1:n])
	set_start_value.(f,fh*ones(n))
	@variable(system_id, p[i = 1:n, j = 1:n_bin])

	for i in 1:n-1
		for j in i+1:n
			@constraint(system_id, Lm[i,j] <= Zro)
			@constraint(system_id, Lm[j,i] <= Zro)
		end
	end
	for i in 1:n
		@constraint(system_id, sum(Lm[i,:]) == 0)
	end
	
	@NLexpression(system_id, err[i = 1:n, t = 1:Dt-1, j = 1:n_bin], Xs[n+i,(j-1)*Dt + t + 1] - Xs[n+i,(j-1)*Dt + t] - dt * (sum(-Lm[i,k]*Xs[k,(j-1)*Dt + t] for k = 1:n) - dm[i]*Xs[n+i,(j-1)*Dt + t] + a[i]*cos(2*pi*f[i]*((j-1)*Dt + t)*dt + p[i,j])))

	@NLobjective(system_id, Min, sum(err[i,t,j]^2 for i = 1:n for t = 1:Dt-1 for j = 1:n_bin)/T)

	optimize!(system_id)

	return value.(Lm), value.(dm), value.(a), value.(f), value.(p)
end


# Optimizes the dynamics matrix and forcing's paramters on the whole time series. 
#
# INPUT:
# Xs: Time series of the phase angles (rows 1:n) and of the phase frequencies (n+1:2*n).
# dt: Time step.
# Lh: Estimate of the normalized Laplacian matrix.
# dh: Estimate of the damping over inertia ratios.
# ah: Estimate of the forcing's amplitude.
# fh: Estimate of the forcing's frequency.
# mu: Initial value of the barrier parameter (in IPOPT).
# bp: Initial value of the bound_push parameter (in IPOPT).
# Zro: Tolerance for zero.
#
# OUTPUT:
# Lm: Optimized normalized Laplacian matrix.
# dm: Optimized damping of inertia ratios.
# a: Optimized forcing amplitudes.
# f: Optimized forcing frequencies.
# p: Optimized phase lags.
function optim_all(Xs::Array{Float64,2}, dt::Float64, Lh::Array{Float64,2}, dh::Array{Float64,1}, ah::Array{Float64,1}, fh::Array{Float64,1}, mu::Float64=1e-1, bp::Float64=1e-1, Zro::Float64=1e-5)
	nn,T = size(Xs)
	n = Int(nn/2)

	system_id = Model(with_optimizer(Ipopt.Optimizer, mu_init = mu, bound_push = bp))

	@variable(system_id, Lm[i = 1:n, j = 1:n])
	set_start_value.(Lm,Lh)
	@variable(system_id, dm[i = 1:n])
	set_start_value.(dm,dh)
	@variable(system_id, a[i = 1:n])
	set_start_value.(a,ah)
	@variable(system_id, f[i = 1:n])
	set_start_value.(f,fh)
	@variable(system_id, p[i = 1:n])

	for i in 1:n-1
		for j in i+1:n
			@constraint(system_id, Lm[i,j] <= Zro)
			@constraint(system_id, Lm[j,i] <= Zro)
		end
	end
	for i in 1:n
		@constraint(system_id, sum(Lm[i,:]) == 0)
	end

	@NLexpression(system_id, err[i = 1:n, t = 1:T-1], Xs[n+i,t+1] - Xs[n+i,t] - dt * (sum(-Lm[i,k]*Xs[k,t] for k = 1:n) - dm[i]*Xs[n+i,t] + a[i]*cos(2*pi*f[i]*t*dt + p[i])))

	@NLobjective(system_id, Min, sum(err[i,t]^2 for i = 1:n for t = 1:T-1)/T)

	optimize!(system_id)

	return value.(Lm), value.(dm), value.(a), value.(f), value.(p)
end


# Optimizes the phase lag of the forcing.
#
# INPUT:
# Xs: Time series of the phase angles (rows 1:n) and of the phase frequencies (n+1:2*n).
# dt: Time step.
# Lm: Estimate of the normalized Laplacian matrix.
# dm: Estimate of the damping over inertia ratios.
# a: Estimate of the forcing's amplitude.
# f: Estimate of the forcing's frequency.
# mu: Initial value of the barrier parameter (in IPOPT).
# bp: Initial value of the bound_push parameter (in IPOPT).
#
# OUTPUT: 
# p: Vector of phases of the forcing.

function optim_phase(Xs::Array{Float64,2}, dt::Float64, Lm::Array{Float64,2}, dm::Array{Float64,1}, a::Array{Float64,1}, f::Array{Float64,1}, mu::Float64=1e-1, bp::Float64=1e-1)
	n = length(dm)

	A = diagm(0 => ones(2*n)) + dt * [zeros(n,n) diagm(0 => ones(n));-Lm -diagm(0 => dm)]
	AX = A*Xs

	phase_id = Model(with_optimizer(Ipopt.Optimizer, mu_init = mu, bound_push = bp))

	@variable(phase_id, c[i = 1:n])
	@variable(phase_id, s[i = 1:n])

	@objective(phase_id, Min, sum((Xs[n+i,t+1] - AX[n+i,t] - dt * (a[i]*c[i]*cos(f[i]*2*pi*dt*t) - a[i]*s[i]*sin(f[i]*2*pi*dt*t))) * (Xs[n+i,t+1] - AX[n+i,t] - dt * (a[i]*c[i]*cos(f[i]*2*pi*dt*t) - a[i]*s[i]*sin(f[i]*2*pi*dt*t))) for i = 1:n for t = 1:T-1))

	optimize!(phase_id)

	return value.(p)
end


# Optimizes the forcing's amplitude and phase lag, assuming known dynamics matrix and frequency. These assumptions allow to write the problem as a Linear program, i.e., accelerates it a lot and makes it more scalable. 
#
# INPUT:
# Xs: Time series of the phase angles (rows 1:n) and of the phase frequencies (n+1:2*n).
# dt: Time step.
# Lh: Estimate of the normalized Laplacian matrix.
# dh: Estimate of the damping over inertia ratios.
# f: Estimate of the forcing's frequency.
# mu: Initial value of the barrier parameter (in IPOPT).
# bp: Initial value of the bound_push parameter (in IPOPT).
#
# OUTPUT:
# a: Optimized forcing amplitudes.
# p: Optimized phase lags.

function optim_all_lin(Xs::Array{Float64,2}, dt::Float64, Lh::Array{Float64,2}, dh::Array{Float64,1}, f::Float64, mu::Float64=1e-1, bp::Float64=1e-1)
	nn,T = size(Xs)
	n = Int(nn/2)

	Ah = diagm(0 => ones(nn)) + dt*[zeros(n,n) diagm(0 => ones(n));-Lh -diagm(0 => dh)]
	AX = Ah*Xs

	system_id = Model(with_optimizer(Ipopt.Optimizer, mu_init = mu, bound_push = bp))

	@variable(system_id, ac[i = 1:n])
	@variable(system_id, as[i = 1:n])

	@objective(system_id, Min, sum((Xs[n+i,t+1] - AX[n+i,t] - dt * (ac[i]*cos(f*2*pi*dt*t) - as[i]*sin(f*2*pi*dt*t))) * (Xs[n+i,t+1] - AX[n+i,t] - dt * (ac[i]*cos(f*2*pi*dt*t) - as[i]*sin(f*2*pi*dt*t))) for i = 1:n for t = 1:T-1))

	optimize!(system_id)
	
	a = sqrt.(value.(ac).^2 + value.(as).^2)
	p = angle.(value.(ac) + im.*value.(as))
	
	return a,p
end






