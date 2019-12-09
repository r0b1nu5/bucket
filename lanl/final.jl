using DelimitedFiles, PyPlot, FFTW, Statistics, LinearAlgebra, JuMP, Ipopt, Distributed



"""
    run_location_small_ntw(Xs::Array{Float64,2}, dt::Float64, Df::Int64=10, n_period::Float64=1., mu::Float64=1e-5, bp::Float64=1e-5, Zro::Float64=1e-5)

Identifies dynamics and forcing for a small network, base only on measurements.

_INPUT_:
`Xs`: Time series of the phase angles (rows 1:n) and of the phase frequencies (n+1:2*n).
`dt`: Time step.
`Df`: Radius of data over which the median of neighboring Fourier modes in computed. 
`n_period`: Number of period considered in each data chunk.
`mu`: Initial value of the barrier parameter (in IPOPT).
`bp`: Initial value of the bound_push parameter (in IPOPT).
`Zro`: Tolerance for zero.

_OUTPUT_:
`Lm`: Laplacian matrix normalized by the inertias.
`dm`: Damping over inertia ratios.
`a`: Vector of amplitudes of the forcing.
`f`: Vector of frequencies of the forcing.
`p`: Vector of phases of the forcing.
"""
function run_location_small_ntw(Xs::Array{Float64,2}, dt::Float64, Df::Int64=10, n_period::Float64=1., mu::Float64=1e-5, bp::Float64=1e-5, Zro::Float64=1e-5)
	nn,T = size(Xs)
	n = Int(nn/2)

	if T%2 == 0
		Xs = Xs[:,1:end-1]
		nn,T = size(Xs)
	end
	
	fh,amps = get_fh_fourier(Xs,dt,Df)
	
	Ah,Lh,dh = get_Ah_correl(Xs,dt,fh)

	ah = get_ah(Xs,dt,fh)

	Lm,dm,a,f,ps = optim_chunks(Xs,dt,Lh,dh,ah,fh,n_period,mu,bp,Zro)

	p = optim_phase(Xs,dt,Lm,dm,a,f)

	AA = sortslices([abs.(a) 1:n],dims=1,rev=true)
	PP = AA[:,1]./sum(AA[:,1])
	id1 = Int(AA[1,2])
	id2 = Int(AA[2,2])
	id3 = Int(AA[3,2])
	
	@info "========================================================================="
	@info "1. Forcing index: $id1, a=$(round(a[id1],digits=3)), f=$(round(f[id1],digits=3)), φ=$(round(p[id1]/pi,digits=2))*π,  (confidence: $(round(PP[1],digits=3)*100)%)"
	@info "2. Forcing index: $id2, a=$(round(a[id2],digits=3)), f=$(round(f[id2],digits=3)), φ=$(round(p[id2]/pi,digits=2))*π,  (confidence: $(round(PP[2],digits=3)*100)%)"
	@info "3. Forcing index: $id3, a=$(round(a[id3],digits=3)), f=$(round(f[id3],digits=3)), φ=$(round(p[id3]/pi,digits=2))*π,  (confidence: $(round(PP[3],digits=3)*100)%)"
	@info "========================================================================="

	return Lm,dm,a,f,p
end

"""
    run_location_large_ntw(Xs::Array{Float64,2}, dt::Float64, n_ref::Int64, Df::Int64=10, n_period::Float64=1., mu::Float64=1e-5, bp::Float64=1e-5, Zro::Float64=1e-5)

Identifies dynamics and forcing for a small network, base only on measurements.

_INPUT_:
`Xs`: Time series of the phase angles (rows 1:n) and of the phase frequencies (n+1:2*n).
`dt`: Time step.
`n_ref`: Number of nodes to use in the reduced network used in the nonlinear optimization.
`Df`: Radius of data over which the median of neighboring Fourier modes in computed. 
`n_period`: Number of period considered in each data chunk.
`mu`: Initial value of the barrier parameter (in IPOPT).
`bp`: Initial value of the bound_push parameter (in IPOPT).
`Zro`: Tolerance for zero.

_OUTPUT_:
`Lm`: Laplacian matrix normalized by the inertias.
`dm`: Damping over inertia ratios.
`a`: Vector of amplitudes of the forcing.
`f`: Vector of frequencies of the forcing.
`p`: Vector of phases of the forcing.
"""
function run_location_large_ntw(Xs::Array{Float64,2}, dt::Float64, n_ref::Int64=5, Df::Int64=10, n_period::Float64=1., mu::Float64=1e-5, bp::Float64=1e-5, Zro::Float64=1e-5)
	nn,T = size(Xs)
	n = Int(nn/2)
	
	if T%2 == 0
		Xs = Xs[:,1:end-1]
		nn,T = size(Xs)
	end
	
	fh,amps = get_fh_fourier(Xs,dt,Df)
	
	iids = Int.(sortslices([amps 1:n],dims=1,rev=true)[:,2])
	ids = [iids[1],]
	ri = 1
	while length(ids) < n_ref
		tryy = true
		while tryy
			tryy = false
			ri += 1
			for k in 1:length(ids)
				tryy = tryy || Xs[ids[k],:] == Xs[iids[ri],:] || Xs[ids[k]+n,:] == Xs[iids[ri]+n,:]
			end
		end
		push!(ids,iids[ri])
	end

	#ids = Int.(sortslices([amps 1:n],dims=1,rev=true)[1:min(n_ref,n),2])
	test = 

	Ah,Lh,dh = get_Ah_correl(Xs[[ids;ids.+n],:],dt,fh)

	ah = get_ah(Xs[[ids;ids.+n],:],dt,fh)

	Lm,dm,a,f,ps = optim_chunks(Xs[[ids;ids.+n],:],dt,Lh,dh,ah,fh,n_period,mu,bp,Zro)
	AA = sortslices([abs.(a) ids 1:n_ref],dims=1,rev=true)
	fh = f[Int(AA[1,3])]
	ff = fh*ones(n)
	ff[ids] = f
	
	Ah,Lh,dh = get_Ah_correl(Xs,dt,fh)
	
	a,p = optim_all_lin(Xs,dt,Lh,dh,fh,mu,bp)

	AA = sortslices([abs.(a) 1:n],dims=1,rev=true)
	PP = AA[:,1]./sum(AA[:,1])
	id1 = Int(AA[1,2])
	id2 = Int(AA[2,2])
	id3 = Int(AA[3,2])
	id4 = Int(AA[4,2])
	id5 = Int(AA[5,2])

	@info "========================================================================="
	@info "1. Forcing index: $id1, a=$(round(a[id1],digits=3)), f=$(round(ff[id1],digits=3)), φ=$(round(p[id1]/pi,digits=2))*π,  (confidence: $(round(PP[1],digits=3)*100)%)"
	@info "2. Forcing index: $id2, a=$(round(a[id2],digits=3)), f=$(round(ff[id2],digits=3)), φ=$(round(p[id2]/pi,digits=2))*π,  (confidence: $(round(PP[2],digits=3)*100)%)"
	@info "3. Forcing index: $id3, a=$(round(a[id3],digits=3)), f=$(round(ff[id3],digits=3)), φ=$(round(p[id3]/pi,digits=2))*π,  (confidence: $(round(PP[3],digits=3)*100)%)"
	@info "4. Forcing index: $id4, a=$(round(a[id4],digits=3)), f=$(round(ff[id4],digits=3)), φ=$(round(p[id4]/pi,digits=2))*π,  (confidence: $(round(PP[4],digits=3)*100)%)"
	@info "5. Forcing index: $id5, a=$(round(a[id5],digits=3)), f=$(round(ff[id5],digits=3)), φ=$(round(p[id5]/pi,digits=2))*π,  (confidence: $(round(PP[5],digits=3)*100)%)"
	@info "========================================================================="

	return Lh,dh,a,ff,p
end

"""
    get_fh_fourier(Xs::Array{Float64,2}, dt::Float64, Df::Int64=10)

Estimates the forcings frequency, based on the Fourier Transform of the times series of the phase frequencies. A peak is identified by dividing the value of each Fourier mode by the median value of its 2*Df neighbors.

_INPUT_:
`Xs`: Time series of the phase angles (rows 1:n) and of the phase frequencies (n+1:2*n).
`dt`: Time step.
`Df`: Radius of data over which the median in computed. 

_OUTPUT_:
`fh`: Estimate of the forcing's frequency.
"""
function get_fh_fourier(Xs::Array{Float64,2}, dt::Float64, Df::Int64=10)
	nn,T = size(Xs)
	n = Int(nn/2)

	fX = zeros(Complex{Float64},n,T)
	for i in 1:n
		fX[i,:] = fft(Xs[n+i,:]).*dt./pi
	end
	
	nfX = [zeros(n) norm.(fX[:,2:end])]
	
	mfX = zeros(n,T)
	for i in 1:n
		for j in 2:T
			mfX[i,j] = nfX[i,j]/max(1e-10,median([nfX[i,max(j-Df,1):j-2];nfX[i,j+2:min(j+Df,T)]]))
		end
	end

	fs = (0:T-1)./(dt*T)
	df = fs[2]-fs[1]

	freqs = Array{Float64,1}()
	maxs = Array{Float64,1}()
	for i in 1:n
		ma,id = findmax(mfX[i,1:ceil(Int,T/2)])

		ma2,id2 = findmax([mfX[i,1:id-1];0.;mfX[i,id+1:ceil(Int,T/2)]])

		if (ma2 > .5*ma) && (abs(id - id2) == 1)
			push!(freqs,mean([fs[id],fs[id2]]))
			push!(maxs,ma)
		else
			push!(freqs,fs[id])
			push!(maxs,ma)
		end
	end

	fh = median(freqs)
	tf = sum((abs.(freqs .- fh) .> 2*df))

	if tf/n > .05
		ids = sortslices([maxs 1:n],dims=1,rev=true)[:,2]
		fh2 = median(freqs[Int.(ids[1:10])])
		tf2 = sum((abs.(freqs[Int.(ids[1:10])] .- fh2) .> 2*df))

		if tf2/n > .05
			@info "WARNING: Fourier Transform: no clear result."
		else
			@info "WARNING: Fourier Transform: result not completely clear."
		end	
	end
	
	return fh,maxs
end


"""
    get_fh_autocorr(Xs::Array{Float64,2}, dt::Float64, p::Float64=.1)

Estimates the forcing's frequency based on the autocorrelation of the times series of the phase angles.

_INPUT_:
`Xs`: Time series of the phase angles (rows 1:n) and of the phase frequencies (n+1:2*n).
`dt`: Time step.
`p`: Proportion of the measure time over which the autocorrelation is computed.

_OUTPUT_: 
`fh`: Estimate of the forcing's frequency.
"""
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


"""
    get_Ah_correl(Xs::Array{Float64,2}, dt::Float64, fh::Float64)

Estimates the dynamics matrix bases on Lokhov18.

_INPUT_:
`Xs`: Time series of the phase angles (rows 1:n) and of the phase frequencies (n+1:2*n).
`dt`: Time step.
`fh`: Estimate of forcing's frequency.

_OUTPUT_: 
`Ah`: Estimate of the full dynamics matrix.
`Lh`: Estimate of the Laplacian matrix normalized by the inertias (M^{-1}*L).
`dh`: Estimate of the damping over inertia ratios.
"""
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


"""
    filter_signal(Xs::Array{Float64,2}, dt::Float64, fh::Float64)

Removes the Fourier modes of frequency fh from the signal Xs.

_INPUT_: 
`Xs`: Time series of the phase angles (rows 1:n) and of the phase frequencies (n+1:2*n).
`dt`: Time step.
`fh`: Frequency to remove.

_OUPUT_:
`Xt`: Filtered signal.
"""
function filter_signal(Xs::Array{Float64,2}, dt::Float64, fh::Float64)
	nn,T = size(Xs)
	nf = floor(Int,T/2)

	fs = (0:T-1)./(dt*T)

	ta = ceil(Int64,T*5e-5)

	mi,id = findmin(abs.(fs .- fh))

	fX = zeros(Complex{Float64},nn,T)
	fXt = zeros(Complex{Float64},nn,T)
	Xt = zeros(nn,T)
	
	for i in 1:nn
		fX[i,:] = fftshift(fft(Xs[i,:]))
		if nf-id-ta >= 0 && nf-id+ta+3 <= nf+id-ta-1
			fXt[i,:] = [fX[i,1:nf+2-id+ta-2*ta-1];zeros(2*ta+1);fX[i,nf+2-id+ta+1:nf+id-ta-1];zeros(2*ta+1);fX[i,nf+id-ta+2*ta+1:end]]
		elseif nf-id-ta >= 0
			fXt[i,:] = [fX[i,1:nf+2-id+ta-2*ta-1];zeros(2*id+2*ta-1);fX[i,nf+id-ta+2*ta+1:end]]
		else
			fXt[i,:] = [zeros(Int.((T-2*id+2*ta+3)/2));fX[i,nf+2-id+ta+1:nf+id-ta-1];zeros(Int.((T-2*id+2*ta+3)/2))]
		end

		Xt[i,:] = real.(ifft(ifftshift(fXt[i,:])))
	end

	return Xt
end


"""
    get_ah(Xs::Array{Float64,2}, dt::Float64, fh::Float64)

Estimates the amplitude of the forcing. This is a very rough estimate, but we hope it is better than nothing.

_INPUT_:
`Xs`: Time series of the phase angles (rows 1:n) and of the phase frequencies (n+1:2*n).
`dt`: Time step.
`fh`: Estimate of the forcing's frequency.

_OUTPUT_:
`ah`: Estimate of the forcing's amplitude.
"""
function get_ah(Xs::Array{Float64,2}, dt::Float64, fh::Float64)
	nn,T = size(Xs)
	n = Int(nn/2)

	ah = maximum(abs.(Xs[1:n,:])) * sqrt(n) * (2*pi*fh)

	return ah
end


"""
    optim_chunks(Xs::Array{Float64,2}, dt::Float64, Lh::Array{Float64,2}, dh::Array{Float64,1}, ah::Float64, fh::Float64, n_period::Float64, mu::Float64=1e-1, bp::Float64=1e-1, Zro::Float64=1e-5)

Optimizes the dynamics matrix and forcing's paramters on chunks of the time series. The optimization allows different phase lags on each chunk, but requires the same amplitude and frequency.

_INPUT_:
`Xs`: Time series of the phase angles (rows 1:n) and of the phase frequencies (n+1:2*n).
`dt`: Time step.
`Lh`: Estimate of the normalized Laplacian matrix.
`dh`: Estimate of the damping over inertia ratios.
`ah`: Estimate of the forcing's amplitude.
`fh`: Estimate of the forcing's frequency.
`n_period`: Number of period considered in each data chunk.
`mu`: Initial value of the barrier parameter (in IPOPT).
`bp`: Initial value of the bound_push parameter (in IPOPT).
`Zro`: Tolerance for zero.

_OUTPUT_:
`Lm`: Optimized normalized Laplacian matrix.
`dm`: Optimized damping of inertia ratios.
`a`: Optimized forcing amplitudes.
`f`: Optimized forcing frequencies.
`p`: Optimized phase lags in each chunk.
"""
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


"""
    optim_all(Xs::Array{Float64,2}, dt::Float64, Lh::Array{Float64,2}, dh::Array{Float64,1}, ah::Array{Float64,1}, fh::Array{Float64,1}, mu::Float64=1e-1, bp::Float64=1e-1, Zro::Float64=1e-5)

Optimizes the dynamics matrix and forcing's paramters on the whole time series. 

_INPUT_:
`Xs`: Time series of the phase angles (rows 1:n) and of the phase frequencies (n+1:2*n).
`dt`: Time step.
`Lh`: Estimate of the normalized Laplacian matrix.
`dh`: Estimate of the damping over inertia ratios.
`ah`: Estimate of the forcing's amplitude.
`fh`: Estimate of the forcing's frequency.
`mu`: Initial value of the barrier parameter (in IPOPT).
`bp`: Initial value of the bound_push parameter (in IPOPT).
`Zro`: Tolerance for zero.

_OUTPUT_:
`Lm`: Optimized normalized Laplacian matrix.
`dm`: Optimized damping of inertia ratios.
`a`: Optimized forcing amplitudes.
`f`: Optimized forcing frequencies.
`p`: Optimized phase lags.
"""
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


"""
    optim_chunks_l0(id::Int64, Xs::Array{Float64,2}, dt::Float64, Lh::Array{Float64,2}, dh::Array{Float64,1}, ah::Float64, fh::Float64, n_period::Float64, mu::Float64=1e-1, bp::Float64=1e-1, Zro::Float64=1e-5)

Optimizes the dynamics matrix and forcing's paramters on chunks of the time series, assuming the forcing is at node id. The optimization allows different phase lags on each chunk, but requires the same amplitude and frequency.

_INPUT_:
`id`: Index of the node where the forcing is assumed to be.
`Xs`: Time series of the phase angles (rows 1:n) and of the phase frequencies (n+1:2*n).
`dt`: Time step.
`Lh`: Estimate of the normalized Laplacian matrix.
`dh`: Estimate of the damping over inertia ratios.
`ah`: Estimate of the forcing's amplitude.
`fh`: Estimate of the forcing's frequency.
`n_period`: Number of period considered in each data chunk.
`mu`: Initial value of the barrier parameter (in IPOPT).
`bp`: Initial value of the bound_push parameter (in IPOPT).
`Zro`: Tolerance for zero.

_OUTPUT_:
`Lm`: Optimized normalized Laplacian matrix.
`dm`: Optimized damping of inertia ratios.
`a`: Optimized forcing amplitudes.
`f`: Optimized forcing frequencies.
`p`: Optimized phase lags in each chunk.
`obj`: Value of the objective function.
"""
function optim_chunks_l0(id::Int64, Xs::Array{Float64,2}, dt::Float64, Lh::Array{Float64,2}, dh::Array{Float64,1}, ah::Float64, fh::Float64, n_period::Float64, mu::Float64=1e-1, bp::Float64=1e-1, Zro::Float64=1e-5)
	nn,T = size(Xs)
	n = Int(nn/2)
	
	Dt = ceil(Int64,1/(fh*dt))
	n_bin = floor(Int64,T/Dt)

	system_id = Model(with_optimizer(Ipopt.Optimizer, mu_init = mu, bound_push = bp))

	@variable(system_id, Lm[i = 1:n, j = 1:n])
	set_start_value.(Lm,Lh)
	@variable(system_id, dm[i = 1:n])
	set_start_value.(dm,dh)
	@variable(system_id, a)
	set_start_value(a,ah)
	@variable(system_id, f)
	set_start_value(f,fh)
	@variable(system_id, p[j = 1:n_bin])

	for i in 1:n-1
		for j in i+1:n
			@constraint(system_id, Lm[i,j] >= Zro)
			@constraint(system_id, Lm[j,i] >= Zro)
		end
	end
	for i in 1:n
		@constraint(system_id, sum(Lm[i,:]) == 0)
	end

	@NLexpression(system_id, err[i = [1:id-1;id+1:n], t = 1:Dt-1, j = 1:n_bin], Xs[n+i,(j-1)*Dt + t + 1] - Xs[n+i,(j-1)*Dt + t] - dt * (sum(-Lm[i,k]*Xs[k,(j-1)*Dt + t] for k = 1:n) - dm[i]*Xs[n+i,(j-1)*Dt + t]))
	@NLexpression(system_id, erri[t = 1:Dt-1, j = 1:n_bin], Xs[n+id,(j-1)*Dt + t + 1] - Xs[n+id,(j-1)*Dt + t] - dt * (sum(-Lm[id,k]*Xs[k,(j-1)*Dt + t] for k = 1:n) - dm[id]*Xs[n+id,(j-1)*Dt + t] + a*cos(2*pi*f*dt*((j-1)*Dt + t) + p[j])))

	@NLobjective(system_id, Min, (sum(err[i,t,j]^2 for i = [1:id-1;id+1:n] for t = 1:Dt-1 for j = 1:n_bin) + sum(erri[t,j]^2 for t = 1:Dt-1 for j = 1:n_bin))/T)

	optimize!(system_id)

	return value.(Lm),value.(dm),value(a),value(f),value.(p),objective_value(system_id) 
end

"""
Intermediate function.
"""
function interim(x::Tuple{Int64, Array{Float64,2}, Float64, Array{Float64,2}, Array{Float64,1}, Float64, Float64, Float64, Float64, Float64, Float64})
	return optim_chunks_l0(x[1],x[2],x[3],x[4],x[5],x[6],x[7],x[8],x[9],x[10],x[11])
end


"""
    parallel_optim(n_thread::Int64, Xs::Array{Float64,2}, dt::Float64, Lh::Array{Float64,2}, dh::Array{Float64,1}, ah::Float64, fh::Float64, n_period::Float64, mu::Float64=1e-1, bp::Float64=1e-1, Zro::Float64=1e-5)

Runs the function `optim_chunks_l0` in parallel for each id = 1:n.
   
_INPUT_:
`n_thread`: Number of threads to parallelize on.
`Xs`: Time series of the phase angles (rows 1:n) and of the phase frequencies (n+1:2*n).
`dt`: Time step.
`Lh`: Estimate of the normalized Laplacian matrix.
`dh`: Estimate of the damping over inertia ratios.
`ah`: Estimate of the forcing's amplitude.
`fh`: Estimate of the forcing's frequency.
`n_period`: Number of period considered in each data chunk.
`mu`: Initial value of the barrier parameter (in IPOPT).
`bp`: Initial value of the bound_push parameter (in IPOPT).
`Zro`: Tolerance for zero.

_OUTPUT_:
`sols`: Array of the outputs of `optim_chunks_l0` for each id = 1:n.
"""
function parallel_optim(n_thread::Int64, Xs::Array{Float64,2}, dt::Float64, Lh::Array{Float64,2}, dh::Array{Float64,1}, ah::Float64, fh::Float64, n_period::Float64, mu::Float64=1e-1, bp::Float64=1e-1, Zro::Float64=1e-5)
	n = length(dh)
	nw = nworkers()

	if n_thread > nw
		addprocs(n_thread - nw)
	end

	@everywhere include("final.jl")

	args = Array{Tuple{Int64, Array{Float64,2}, Float64, Array{Float64,2}, Array{Float64,1}, Float64, Float64, Float64, Float64, Float64, Float64},1}()
	for i in 1:n
		push!(args,(i,Xs,dt,Lh,dh,ah,fh,n_period,mu,bp,Zro))
	end
	sols = pmap(interim,args)
#	sols = pmap(optim_chunks_l0,1:n,[Xs,],[dt,],[Lh,],[dh,],[ah,],[dh,],[n_period,],[mu,],[bp,],[Zro,])
	return sols
end



"""
    optim_phase(Xs::Array{Float64,2}, dt::Float64, Lm::Array{Float64,2}, dm::Array{Float64,1}, a::Array{Float64,1}, f::Array{Float64,1}, mu::Float64=1e-1, bp::Float64=1e-1)

Optimizes the phase lag of the forcing.

_INPUT_:
`Xs`: Time series of the phase angles (rows 1:n) and of the phase frequencies (n+1:2*n).
`dt`: Time step.
`Lm`: Estimate of the normalized Laplacian matrix.
`dm`: Estimate of the damping over inertia ratios.
`a`: Estimate of the forcing's amplitude.
`f`: Estimate of the forcing's frequency.
`mu`: Initial value of the barrier parameter (in IPOPT).
`bp`: Initial value of the bound_push parameter (in IPOPT).

_OUTPUT_: 
`p`: Vector of phases of the forcing.
"""
function optim_phase(Xs::Array{Float64,2}, dt::Float64, Lm::Array{Float64,2}, dm::Array{Float64,1}, a::Array{Float64,1}, f::Array{Float64,1}, mu::Float64=1e-1, bp::Float64=1e-1)
	nn,T = size(Xs)
	n = length(dm)

	A = diagm(0 => ones(2*n)) + dt * [zeros(n,n) diagm(0 => ones(n));-Lm -diagm(0 => dm)]
	AX = A*Xs

	phase_id = Model(with_optimizer(Ipopt.Optimizer, mu_init = mu, bound_push = bp))

	@variable(phase_id, c[i = 1:n])
	@variable(phase_id, s[i = 1:n])

	@objective(phase_id, Min, sum((Xs[n+i,t+1] - AX[n+i,t] - dt * (a[i]*c[i]*cos(f[i]*2*pi*dt*t) - a[i]*s[i]*sin(f[i]*2*pi*dt*t))) * (Xs[n+i,t+1] - AX[n+i,t] - dt * (a[i]*c[i]*cos(f[i]*2*pi*dt*t) - a[i]*s[i]*sin(f[i]*2*pi*dt*t))) for i = 1:n for t = 1:T-1))

	optimize!(phase_id)

	p = angle.(value.(c) + im*value.(s))

	return p
end


"""
    optim_all_lin(Xs::Array{Float64,2}, dt::Float64, Lh::Array{Float64,2}, dh::Array{Float64,1}, f::Float64, mu::Float64=1e-1, bp::Float64=1e-1)

Optimizes the forcing's amplitude and phase lag, assuming known dynamics matrix and frequency. These assumptions allow to write the problem as a Linear program, i.e., accelerates it a lot and makes it more scalable. 

_INPUT_:
`Xs`: Time series of the phase angles (rows 1:n) and of the phase frequencies (n+1:2*n).
`dt`: Time step.
`Lh`: Estimate of the normalized Laplacian matrix.
`dh`: Estimate of the damping over inertia ratios.
`f`: Estimate of the forcing's frequency.
`mu`: Initial value of the barrier parameter (in IPOPT).
`bp`: Initial value of the bound_push parameter (in IPOPT).

_OUTPUT_:
`a`: Optimized forcing amplitudes.
`p`: Optimized phase lags.
"""
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






