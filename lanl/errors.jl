using LinearAlgebra, JuMP, Ipopt, Dates, PyPlot, FFTW, Dates

function errors_uk_inhomog(ids::Array{Int64,1}, a::Array{Float64,1}, f::Float64, phi::Array{Float64,1}, T::Int64, dt::Float64, eps::Float64 = 1e-5)
	n = 120	
	I = diagm(0 => ones(n))
	
	eA2 = Array{Float64,1}()
	eAi = Array{Float64,1}()
	ef = Array{Float64,1}()
	ea2 = Array{Float64,1}()
	eai = Array{Float64,1}()
	ep2 = Array{Float64,1}()
	epi = Array{Float64,1}()
	
	for id in ids
		L = readdlm("data/uk$(id)_L.csv",',')
		m = vec(readdlm("data/uk$(id)_m.csv",','))
		d = vec(readdlm("data/uk$(id)_d.csv",','))
		
		Mi = diagm(0 => 1 ./ m)
		D = diagm(0 => d)
		
		A = [I dt*I; (-dt*Mi*L) (I - dt*Mi*D)]
		
		Xs = readdlm("data/uk$(id)_forced_$(f)_$(T)_$(dt).csv",',')
		
		Ah, fh = find_A_n_f(Xs, dt)
		Ahr = zeros(2*n,2*n)
		for i in 1:2*n
			for j in 1:2*n
				if abs(Ah[i,j]) > eps
					Ahr[i,j] = Ah[i,j]
				end
			end
		end
		
		ah, ph = locate_f_n_lag(Xs, dt, fh, Ah)
		ahr = zeros(n)
		phr = zeros(n)
		for i in 1:n
			if abs(ah[i])/maximum(abs.(ah)) > .01
				ahr[i] = ah[i]
				phr[i] = ph[i]
			end
		end
		
		push!(eA2, sqrt(sum((Ahr - A).^2))/sqrt(sum(A.^2)))
		push!(eAi, maximum(abs.(Ahr - A))/maximum(abs.(A)))
		
		push!(ef, abs(fh - f)/abs(f))
		
		push!(ea2, norm(ahr - a)/norm(a))
		push!(eai, maximum(abs.(ahr- a))/maximum(abs.(a)))
		
		push!(ep2, norm(phr - phi)/norm(phi))
		push!(epi, maximum(abs.(phr - phi))/maximum(abs.(phi)))
		
	end
	
	return eA2, eAi, ef, ea2, eai, ep2, epi
end



function errors_uk_dt(dts::Array{Float64,1}, a::Array{Float64,1}, f::Float64, phi::Array{Float64,1}, T::Int64, eps::Float64 = 1e-5)
	n = 120	
	I = diagm(0 => ones(n))
	
	include("load_uk.jl")
	
	Mi = diagm(0 => 1 ./ m)
	D = diagm(0 => d)
		
	eA2 = Array{Float64,1}()
	eAi = Array{Float64,1}()
	ef = Array{Float64,1}()
	ea2 = Array{Float64,1}()
	eai = Array{Float64,1}()
	ep2 = Array{Float64,1}()
	epi = Array{Float64,1}()
	
	for dt in dts
		A = [I dt*I; (-dt*Mi*L) (I - dt*Mi*D)]
		
		Xs = readdlm("data/uk_forced_$(f)_$(T)_$(dt).csv",',')
		
		Ah, fh = find_A_n_f(Xs, dt)
		Ahr = zeros(2*n,2*n)
		for i in 1:2*n
			for j in 1:2*n
				if abs(Ah[i,j]) > eps
					Ahr[i,j] = Ah[i,j]
				end
			end
		end
		
		ah, ph = locate_f_n_lag(Xs, dt, fh, Ah)
		ahr = zeros(n)
		phr = zeros(n)
		for i in 1:n
			if abs(ah[i])/maximum(abs.(ah)) > .01
				ahr[i] = ah[i]
				phr[i] = ph[i]
			end
		end
		
		push!(eA2, sqrt(sum((Ahr - A).^2))/sqrt(sum(A.^2)))
		push!(eAi, maximum(abs.(Ahr - A))/maximum(abs.(A)))
		
		push!(ef, abs(fh - f)/abs(f))
		
		push!(ea2, norm(ahr - a)/norm(a))
		push!(eai, maximum(abs.(ahr- a))/maximum(abs.(a)))
		
		push!(ep2, norm(phr - phi)/norm(phi))
		push!(epi, maximum(abs.(phr - phi))/maximum(abs.(phi)))
		
	end
	
	return eA2, eAi, ef, ea2, eai, ep2, epi
end


function errors_uk_T(Ts::Array{Int64,1}, a::Array{Float64,1}, f::Float64, phi::Array{Float64,1}, dt::Float64, eps::Float64 = 1e-5)
	n = 120	
	I = diagm(0 => ones(n))
	
	include("load_uk.jl")
	
	Mi = diagm(0 => 1 ./ m)
	D = diagm(0 => d)
		
	eA2 = Array{Float64,1}()
	eAi = Array{Float64,1}()
	ef = Array{Float64,1}()
	ea2 = Array{Float64,1}()
	eai = Array{Float64,1}()
	ep2 = Array{Float64,1}()
	epi = Array{Float64,1}()
	
	for T in Ts
		A = [I dt*I; (-dt*Mi*L) (I - dt*Mi*D)]
		
		Xs = readdlm("data/uk_forced_$(f)_$(T)_$(dt).csv",',')
		
		Ah, fh = find_A_n_f(Xs, dt)
		Ahr = zeros(2*n,2*n)
		for i in 1:2*n
			for j in 1:2*n
				if abs(Ah[i,j]) > eps
					Ahr[i,j] = Ah[i,j]
				end
			end
		end
		
		ah, ph = locate_f_n_lag(Xs, dt, fh, Ah)
		ahr = zeros(n)
		phr = zeros(n)
		for i in 1:n
			if abs(ah[i])/maximum(abs.(ah)) > .01
				ahr[i] = ah[i]
				phr[i] = ph[i]
			end
		end
		
		push!(eA2, sqrt(sum((Ahr - A).^2))/sqrt(sum(A.^2)))
		push!(eAi, maximum(abs.(Ahr - A))/maximum(abs.(A)))
		
		push!(ef, abs(fh - f)/abs(f))
		
		push!(ea2, norm(ahr - a)/norm(a))
		push!(eai, maximum(abs.(ahr- a))/maximum(abs.(a)))
		
		push!(ep2, norm(phr - phi)/norm(phi))
		push!(epi, maximum(abs.(phr - phi))/maximum(abs.(phi)))
		
	end
	
	return eA2, eAi, ef, ea2, eai, ep2, epi
end


function errors_uk_dt_T(dts_Ts::Array{Tuple{Float64, Int64},1}, a::Array{Float64,1}, f::Float64, phi::Array{Float64,1}, eps::Float64 = 1e-5)
	n = 120	
	I = diagm(0 => ones(n))
	
	include("load_uk.jl")
	
	Mi = diagm(0 => 1 ./ m)
	D = diagm(0 => d)
		
	eA2 = Array{Float64,1}()
	eAi = Array{Float64,1}()
	ef = Array{Float64,1}()
	ea2 = Array{Float64,1}()
	eai = Array{Float64,1}()
	ep2 = Array{Float64,1}()
	epi = Array{Float64,1}()
	
	for (dt,T) in dts_Ts
		A = [I dt*I; (-dt*Mi*L) (I - dt*Mi*D)]
		
		Xs = readdlm("data/uk_forced_$(f)_$(T)_$(dt).csv",',')
		
		Ah, fh = find_A_n_f(Xs, dt)
		Ahr = zeros(2*n,2*n)
		for i in 1:2*n
			for j in 1:2*n
				if abs(Ah[i,j]) > eps
					Ahr[i,j] = Ah[i,j]
				end
			end
		end
		
		ah, ph = locate_f_n_lag(Xs, dt, fh, Ah)
		ahr = zeros(n)
		phr = zeros(n)
		for i in 1:n
			if abs(ah[i])/maximum(abs.(ah)) > .01
				ahr[i] = ah[i]
				phr[i] = ph[i]
			end
		end
		
		push!(eA2, sqrt(sum((Ahr - A).^2))/sqrt(sum(A.^2)))
		push!(eAi, maximum(abs.(Ahr - A))/maximum(abs.(A)))
		
		push!(ef, abs(fh - f)/abs(f))
		
		push!(ea2, norm(ahr - a)/norm(a))
		push!(eai, maximum(abs.(ahr- a))/maximum(abs.(a)))
		
		push!(ep2, norm(phr - phi)/norm(phi))
		push!(epi, maximum(abs.(phr - phi))/maximum(abs.(phi)))
		
	end
	
	return eA2, eAi, ef, ea2, eai, ep2, epi
end


