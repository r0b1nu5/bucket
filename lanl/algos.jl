using LinearAlgebra, JuMP, Ipopt, Dates, PyPlot, FFTW, Dates

#### STEP 1: SYSTEM IDENTIFICATION AND FORCING'S FREQUENCY #############

## INPUTS
# Xs: Time series of measurments (2*n x T). The first n rows are the angle time series and the n last rows are the frequency time series. T is the number of time steps
# dt: Step size [s].

## OUTPUTS
# Ah: Estimate of the matrix of dynamics.
# fh: Estimate of the frequency of the forcing.

function find_A_n_f(X::Array{Float64,2}, dt::Float64)
	nn,T = size(X)
	n = Int(nn/2)
	
	fX = zeros(Complex{Float64},nn,T)
	for i in 1:nn
		fX[i,:] = fft(X[i,:]).*dt./pi
	end
	
	nfX = norm.(fX)
	fXt = copy(fX)
	
	fs = (0:T-1)./(dt*T)
	
	freqs = Array{Float64,1}()
	maxs = Array{Float64,1}()
	
	for i in 1:n
		ma,id = findmax(nfX[n+i,2:Int(T/2)])
		id += 1
		
		mas = [ma,]
		ids = [id,]
		test = true
		
		fXt[:,id] = zeros(Complex{Float64},nn)
		fXt[:,T-id+2] = zeros(Complex{Float64},nn)
		
		while test == true
			test = false
			for k in [minimum(ids)-1,maximum(ids)+1]
				if 1 < k <= T
					if nfX[n+i,k] > .5*ma
						push!(mas,nfX[n+i,k])
						push!(ids,k)
						fXt[:,k] = zeros(Complex{Float64},nn)
						fXt[:,T-k+2] = zeros(Complex{Float64},nn)
						test = true
					end
				end
			end
		end
		
		push!(freqs, sum(mas.^2 .* fs[ids])./(sum(mas.^2)))
		push!(maxs,ma)
	end
	
	m,k = findmax(maxs)
	fh = freqs[k]

	Xt = zeros(nn,T)
	for i in 1:n
		Xt[i,:] = real.(ifft(fXt[i,:]))
		Xt[i+n,:] = real.(ifft(fXt[i+n,:]))
	end
	
	Ah, Adh = system_identification_correl(Xt, dt)
	
	return Ah, fh
end

#### STEP 2: DETERMINE FORCINGS AMPLITUDE AND PHASE LAG ###############3

#### WARNING !!! THE ACCURACY OF THIS SCRIPT IS VERY SENSITIVE TO ERRORS IN THE FREQUENCY !!! ############

## INPUTS
# Xs: Time series of measurments (2*n x T). The first n rows are the angle time series and the n last rows are the frequency time series. T is the number of time steps
# dt: Step size [s].
# f: Forcing's frequency [s^{-1}].
# A: Dynamics matrix.
# l: Penalty coefficient to the 1-norm of a.

## OUTPUTS
# ah: Estimate of the vector of forcing amplitude.
# ph: Estimate of the phase lag in the forcing.

function locate_f_n_lag(X::Array{Float64,2}, dt::Float64, f::Float64, A::Array{Float64,2}, l::Float64 = .1, verb::Bool = false)
	nn,T = size(X)
	n = Int(nn/2)
	
	AX = A*X
	
	system_id = Model(with_optimizer(Ipopt.Optimizer))
	
	@variable(system_id, 0 <= ac[i = 1:n])
	@variable(system_id, as[i = 1:n])
	
	if verb
		@info "$(now()) -- Variables"
	end
	
	@objective(system_id, Min, sum((X[n+i,t+1] - AX[n+i,t] - dt * (ac[i]*cos(f*2*pi*dt*t) - as[i]*sin(f*2*pi*dt*t))) * (X[n+i,t+1] - AX[n+i,t] - dt * (ac[i]*cos(f*2*pi*dt*t) - as[i]*sin(f*2*pi*dt*t))) for i = 1:n for t = 1:T-1) + l*sum(ac[i] for i = 1:n))
	
	if verb
		@info "$(now()) -- Objective"
	end
	
	optimize!(system_id)
	
	a = sqrt.(value.(ac).^2 + value.(as).^2)
	phi = angle.(value.(ac) + im.*value.(as))

	return a, phi
end


