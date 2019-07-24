using LinearAlgebra, JuMP, Ipopt, Dates, PyPlot, FFTW, Dates

# Based on time series of data, returns the estimated frequency of a forcing term.

## INPUTS
# X: 	Matrix of measurements (2*n x T). The n first rows contain the phase measurements and the n last rows contain the frequency measurements. Each column is a time step.
# dt: 	Step size [s].
# plot: If true, plots the DFT.

## OUTPUTS 
# f: Estimated frequency of the forcing term.

function forcing_frequency(X::Array{Float64,2}, dt::Float64, plot::Bool = false)
	nn,T = size(X)
	n = Int(nn/2)
	
	fX = zeros(Complex,nn,T)
	for i in 1:nn
		fX[i,:] = fft(X[i,:]).*dt./pi
	end
	
	nfX = norm.(fX)
	
	fs = (0:T-1)./(dt*T)
	
	if plot
		for i in 1:n
			PyPlot.plot(fs,nfX[n+i,:])
		end
	end
	
	freqs = Array{Float64,1}()
	maxs = Array{Float64,1}()
	
	for i in 1:n
		ma,id = findmax(nfX[n+i,2:Int(T/2)])
		id += 1
		
		mas = [ma,]
		ids = [id,]
		test = true
		
		while test == true
			test = false
			for k in [minimum(ids)-1,maximum(ids)+1]
				if 1 <= k <= T
					if nfX[n+i,k] > .5*ma
						push!(mas,nfX[n+i,k])
						push!(ids,k)
						test = true
					end
				end
			end
		end
		
		push!(freqs, sum(mas.^2 .* fs[ids])./(sum(mas.^2)))
		push!(maxs,ma)
	end
	
	return sum(maxs.^2 .* freqs)/sum(maxs.^2)
end


# Based on time series of data and on the forcing's frequency, return the matrix of dynamics, and the amplitude and phase lag of the forcing.

## INPUTS
# X: 	Matrix of measurements (2*n x T). The n first rows contain the phase measurements and the n last rows contain the frequency measurements. Each column is a time step.
# dt: 	Step size [s].
# f: 	Frequency of the forcing [s^{-1}].
# l: 	Weight of the 1-norm penalty.

## OUTPUTS
# A: 	Matrix of dynamics (2*n x 2*n).
# a: 	Vector of amplitudes of the forcing.
# phi:	Vector of phase lags. 

function system_id_n_locate_f_n_lag(X::Array{Float64,2}, dt::Float64, f::Float64, l::Float64 = .1)
	nn,T = size(X)
	n = Int(nn/2)
	
	system_id = Model(with_optimizer(Ipopt.Optimizer))
	
	@variable(system_id, MiL[i = 1:n, j = 1:n])
	@variable(system_id, mid[i = 1:n])
	@variable(system_id, 0 <= ac[i = 1:n])
	@variable(system_id, 0 <= as[i = 1:n])
@info "$(now()) -- Variables"
	
	for i in 1:n-1
		for j in i+1:n
			@constraint(system_id, MiL[i,j] <= 0)
			@constraint(system_id, MiL[i,j] <= 0)
			
		end
		@constraint(system_id, MiL[i,i] == -sum(MiL[i,j] for j = 1:i-1) - sum(MiL[i,j] for j = i+1:n))
	end
	@constraint(system_id, MiL[n,n] == -sum(MiL[n,j] for j = 1:n-1))
@info "$(now()) -- Constraints"
	

	@objective(system_id, Min, sum(
(X[n+i,t+1] - X[n+i,t] - dt * (sum(-MiL[i,k]*X[k,t] for k = 1:n) - mid[i]*X[n+i,t] + ac[i]*cos(f*2*pi*dt*t) - as[i]*sin(f*2*pi*dt*t))) * 
(X[n+i,t+1] - X[n+i,t] - dt * (sum(-MiL[i,k]*X[k,t] for k = 1:n) - mid[i]*X[n+i,t] + ac[i]*cos(f*2*pi*dt*t) - as[i]*sin(f*2*pi*dt*t))) 
for i = 1:n for t = 1:T-1) + l*sum(ac[i] for i = 1:n))
@info "$(now()) -- Objective"
	
	optimize!(system_id)
	
	I = diagm(0 => ones(n))
	
	A = [I dt*I; (-dt*value.(MiL)) (I - dt*diagm(0 => value.(mid)))]
	a = sqrt.(value.(ac).^2 + value.(as).^2)
	phi = angle.(value.(ac) + im.*value.(as))

	return A, a, phi
end



