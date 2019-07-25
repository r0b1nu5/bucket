using LinearAlgebra, JuMP, Ipopt, Dates, PyPlot, FFTW, Dates

#### STEP 1: FORCING FREQUENCY #######

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

#### STEP (1-)2: SYSTEM IDENTIFICATION (RETURN THE FORCING'S FREQUENCY AS WELL) #############3

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
#	fh = sum(maxs.^2 .* freqs)/sum(maxs.^2)

	Xt = zeros(nn,T)
	for i in 1:n
		Xt[i,:] = real.(ifft(fXt[i,:]))
		Xt[i+n,:] = real.(ifft(fXt[i+n,:]))
	end
	
	Ah, Adh = system_identification_correl(Xt, dt)
	
	return Ah, fh
end

#### STEP 3: DETERMINE FORCINGS AMPLITUDE AND PHASE LAG ###############3

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

##################### OLD CODES #############################

#### SYSTEM IDENTIFICATION USING CORRELATION MATRIX ########

function system_identification_correl(X::Array{Float64,2}, dt::Float64, l::Float64=.01)
	@info "$(now()) -- Start..."
	
	nn,T = size(X)
	n = Int(nn/2)

	S0 = X[:,1:T-1]*X[:,1:T-1]' ./ (T-1)
	S1 = X[:,2:T]*X[:,1:T-1]' ./ (T-1)
	
	Ah = S1*inv(S0)
	Adh = (Ah - diagm(0 => ones(2*n)))./dt
	
	@info "$(now()) -- Stop."
	
	return Ah,Adh
end


function find_forcing_freq(X::Array{Float64,2}, dt::Float64)
	n,T = size(X)
	
	for i in 1:n
	end
end




#### IF THE FORCING'S FREQUENCY IS MUCH SMALLER THAN d/λ, WE CAN DETERMINE THE PHASE LAG QUITE EFFICIENTLY

function forcing_frequency_n_lag(X::Array{Float64,2}, dt::Float64, plot::Bool = false)
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
	mids = Array{Int64,1}()
	fxs = Array{Complex{Float64},1}()
	
	for i in 1:n
		ma,id = findmax(nfX[n+i,2:Int(T/2)])
		id += 1
		
		mas = [ma,]
		ids = [id,]
		push!(fxs,fX[i,id])
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
	
	w0 = sum(maxs.^2 .* freqs)/sum(maxs.^2)
	
	
	Phis = pi/2 .- abs.(angle.(fxs))
	phi0 = mean(Phis)
	
	
	return w0, phi0
end

#### STEP 2: SYSTEM IDENTIFICATION #######

function system_identification(X::Array{Float64,2}, dt::Float64)
	nn,T = size(X)
	n = Int(nn/2)
	
	system_id = Model(with_optimizer(Ipopt.Optimizer))
	
	@variable(system_id, MiL[i = 1:n, j = 1:n])
	@variable(system_id, mid[i = 1:n])
@info "$(now()) -- Variables"
	
	for i in 1:n-1
		for j in i+1:n
			@constraint(system_id, MiL[i,j] <= 0)
			@constraint(system_id, MiL[j,i] <= 0)
		end
		@constraint(system_id, MiL[i,i] == -sum(MiL[i,j] for j = 1:i-1) - sum(MiL[i,j] for j = i+1:n))
	end
	@constraint(system_id, MiL[n,n] == -sum(MiL[n,j] for j = 1:n-1))
@info "$(now()) -- Constraints"
	
	@objective(system_id, Min, sum((X[n+i,t+1] - X[n+i,t] - dt * (sum(-MiL[i,k]*X[k,t] for k = 1:n) - mid[i]*X[n+i,t])) * (X[n+i,t+1] - X[n+i,t] - dt * (sum(-MiL[i,k]*X[k,t] for k = 1:n) - mid[i]*X[n+i,t])) for i = 1:n for t = 1:T-1))
@info "$(now()) -- Objective"
	
	optimize!(system_id)
@info "$(now()) -- Optimized"
	
	I = diagm(0 => ones(n))
	
	Ah = [I dt*I; (-dt*value.(MiL)) (I - dt*diagm(0 => value.(mid)))]
	
	return Ah
end


#### STEP 3: LOCATING FORCING ###########

function locate_forcing(X::Array{Float64,2}, A::Array{Float64,2}, dt::Float64, f::Float64, nf::Int64 = 1, l::Float64 = 1.)
	nn,T = size(X)
	n = Int(nn/2)
	
	AX = A*X
	
	locate_f = Model(with_optimizer(Ipopt.Optimizer))	

	@variable(locate_f, 0 <= a[i = 1:n] <= 10)
@info "$(now()) -- Variables"	
	
	@objective(locate_f, Min, sum((X[n+i,t+1] - AX[n+i,t] - dt*a[i]*cos(f*2*pi*t*dt))*(X[n+i,t+1] - AX[n+i,t] - dt*a[i]*cos(f*2*pi*t*dt)) for i in 1:n for t in 1:T-1) + l*sum(a[i] for i in 1:n))
@info "$(now()) -- Objective"
		
	optimize!(locate_f)
@info "$(now()) -- Optimized"

	return value.(a)
end



#### STEPS 2-3 ########################

function system_id_n_locate_f(X::Array{Float64,2}, dt::Float64, f::Float64, l::Float64 = .1)
	nn,T = size(X)
	n = Int(nn/2)
	
	system_id = Model(with_optimizer(Ipopt.Optimizer))
	
	@variable(system_id, MiL[i = 1:n, j = 1:n])
	@variable(system_id, mid[i = 1:n])
	@variable(system_id, 0 <= a[i = 1:n])
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
(X[n+i,t+1] - X[n+i,t] - dt * (sum(-MiL[i,k]*X[k,t] for k = 1:n) - mid[i]*X[n+i,t] + a[i]*cos(f*2*pi*dt*t))) * 
(X[n+i,t+1] - X[n+i,t] - dt * (sum(-MiL[i,k]*X[k,t] for k = 1:n) - mid[i]*X[n+i,t] + a[i]*cos(f*2*pi*dt*t))) for i = 1:n for t = 1:T-1) + l*sum(a[i] for i = 1:n))
@info "$(now()) -- Objective"
	
	optimize!(system_id)
	
	I = diagm(0 => ones(n))
	
	A = [I dt*I; (-dt*value.(MiL)) (I - dt*diagm(0 => value.(mid)))]
	a = value.(a)

	return A, a	
end



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



