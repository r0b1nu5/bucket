using LinearAlgebra, JuMP, Ipopt, Dates, PyPlot, FFTW, Dates, DelimitedFiles, Statistics

include("algos.jl")

function errors_uk_inhomog(ids::Array{Int64,1}, a::Array{Float64,1}, f::Float64, phi::Array{Float64,1}, T::Int64, dt::Float64, eps::Float64 = 1e-5)
	n = 120	
	I = diagm(0 => ones(n))

	eAs = Array{Float64,1}()
	Ahs = Array{Array{Float64,2},1}()
	ef = Array{Float64,1}()
	fhs = Array{Float64,1}()
	eas = Array{Float64,1}()
	ahs = Array{Array{Float64,1},1}()
	eps = Array{Float64,1}()
	phs = Array{Array{Float64,1},1}()

	for id in ids
		L = readdlm("data/uk$(id)_L.csv",',')
		m = vec(readdlm("data/uk$(id)_m.csv",','))
		d = vec(readdlm("data/uk$(id)_d.csv",','))
		
		Mi = diagm(0 => 1 ./ m)
		D = diagm(0 => d)
		
		A = [I dt*I; (-dt*Mi*L) (I - dt*Mi*D)]
		
		Xs = readdlm("data2/uk$(id)_forced_$(f)_$(T)_$(dt).csv",',')
		
		Ah, fh = find_A_n_f(Xs, dt)
		
		ah, ph = locate_f_n_lag(Xs, dt, fh, Ah)
		
		eA = zeros(2*n,2*n)
		for i in 1:2*n
			for j in 1:2*n
				if A[i,j] > 0
					eA[i,j] = abs(A[i,j] - Ah[i,j])/abs(A[i,j])
				end
			end
		end
		push!(eAs, eA)
		push!(Ahs, Ahr)
		
		push!(ef, abs(fh - f)/abs(f))
		push!(fhs, fh)

		ea = zeros(n)
		for i in 1:n
			if abs(a[i]) > 0
				ea[i] = abs(a[i] - ah[i])/abs(a[i])
			end
		end
		push!(ahs, ahr)

		ep = zeros(n)
		for i in 1:n
			if abs(a[i]) > 0
				ep[i] = abs(p[i] - ph[i])/(abs(p[i]))
			end
		end
		push!(phs, phr)
	end
	
	return (Ahs, eA), (fhs, ef), (ahs, ea), (phs, ep)
end



function errors_uk_dt(dts::Array{Float64,1}, a::Array{Float64,1}, f::Float64, phi::Array{Float64,1}, T::Int64, plot::Bool = false, eps::Float64 = 1e-5)
	n = 120	
	I = diagm(0 => ones(n))
	
	include("load_uk.jl")
	
	Mi = diagm(0 => 1 ./ m)
	D = diagm(0 => d)
		
	eAs = Array{Array{Float64,1},1}()
	zAs = Array{Array{Float64,1},1}()
	Ahs = Array{Array{Float64,2},1}()
	efs = Array{Float64,1}()
	fhs = Array{Float64,1}()
	eas = Array{Array{Float64,1},1}()
	zas = Array{Array{Float64,1},1}()
	ahs = Array{Array{Float64,1},1}()
	eps = Array{Array{Float64,1},1}()
	zps = Array{Array{Float64,1},1}()
	phs = Array{Array{Float64,1},1}()
	
	for dt in dts
		A = [I dt*I; (-dt*Mi*L) (I - dt*Mi*D)]
		
		Xs = readdlm("data2/uk_forced_$(f)_$(T)_$(dt).csv",',')
		
		Ah, fh = find_A_n_f(Xs, dt)
		
		ah, ph = locate_f_n_lag(Xs, dt, fh, Ah)
		
		eA = Array{Float64,1}()
		zA = Array{Float64,1}()
		for i in 1:2*n
			for j in 1:2*n
				if abs(A[i,j]) > 0
					push!(eA, abs(A[i,j] - Ah[i,j])/abs(A[i,j]))
				else
					push!(zA, abs(Ah[i,j]))
				end
			end
		end
		push!(eAs, eA)
		push!(zAs, zA)

		push!(Ahs, Ah)

		push!(efs, abs(f - fh)/abs(f))
		push!(fhs, fh)

		ea = Array{Float64,1}()
		za = Array{Float64,1}()
		for i in 1:n
			if abs(a[i]) > 0
				push!(ea, abs(a[i] - ah[i])/abs(a[i]))
			else
				push!(za, abs(ah[i]))
			end
		end
		push!(eas, ea)
		push!(zas, za)

		push!(ahs, ah)

		ep = Array{Float64,1}()
		for i in 1:n
			if a[i] > 0
				push!(ep, abs(p[i] - ph[i]))
			end
		end
		push!(eps, ep)

		push!(phs, ph)
	end
	
	if plot
# Error in A
		e000 = Array{Float64,1}()
		e025 = Array{Float64,1}()
		e050 = Array{Float64,1}()
		e075 = Array{Float64,1}()
		e100 = Array{Float64,1}()
		z000 = Array{Float64,1}()
		z025 = Array{Float64,1}()
		z050 = Array{Float64,1}()
		z075 = Array{Float64,1}()
		z100 = Array{Float64,1}()
		
		for i in 1:length(dts)
			push!(e000,minimum(eAs[i]))
			push!(e025,quantile(eAs[i],.25))
			push!(e050,quantile(eAs[i],.5))
			push!(e075,quantile(eAs[i],.75))
			push!(e100,maximum(eAs[i]))
			push!(z000,minimum(zAs[i]))
			push!(z025,quantile(zAs[i],.25))
			push!(z050,quantile(zAs[i],.5))
			push!(z075,quantile(zAs[i],.75))
			push!(z100,maximum(zAs[i]))
		end
		figure()
		PyPlot.plot(dts,e000,"xr")
		PyPlot.plot(dts,e025,"--r")
		PyPlot.plot(dts,e050, "r")
		PyPlot.plot(dts,e075,"--r")
		PyPlot.plot(dts,e100,"xr")
		xlabel("δt")
		ylabel("relative error (non-zero terms)")
		title("UK 45, error in A, (a=$(a[45]), f=$(f), phi=$(p[45]))")
		figure()
		PyPlot.plot(dts,z000,"xr")
		PyPlot.plot(dts,z025,"--r")
		PyPlot.plot(dts,z050, "r")
		PyPlot.plot(dts,z075,"--r")
		PyPlot.plot(dts,z100,"xr")
		xlabel("δt")
		ylabel("error (zero terms)")
		title("UK 45, error in A, (a=$(a[45]), f=$(f), phi=$(p[45]))")
		
# Error in a
		e000 = Array{Float64,1}()
		e025 = Array{Float64,1}()
		e050 = Array{Float64,1}()
		e075 = Array{Float64,1}()
		e100 = Array{Float64,1}()
		z000 = Array{Float64,1}()
		z025 = Array{Float64,1}()
		z050 = Array{Float64,1}()
		z075 = Array{Float64,1}()
		z100 = Array{Float64,1}()
		
		for i in 1:length(dts)
			push!(e000,minimum(eas[i]))
			push!(e025,quantile(eas[i],.25))
			push!(e050,quantile(eas[i],.5))
			push!(e075,quantile(eas[i],.75))
			push!(e100,maximum(eas[i]))
			push!(z000,minimum(zas[i]))
			push!(z025,quantile(zas[i],.25))
			push!(z050,quantile(zas[i],.5))
			push!(z075,quantile(zas[i],.75))
			push!(z100,maximum(zas[i]))
		end
		figure()
		PyPlot.plot(dts,e000,"xb")
		PyPlot.plot(dts,e025,"--b")
		PyPlot.plot(dts,e050, "b")
		PyPlot.plot(dts,e075,"--b")
		PyPlot.plot(dts,e100,"xb")
		xlabel("δt")
		ylabel("relative error (non-zero terms)")
		title("UK 45, error in a, (a=$(a[45]), f=$(f), phi=$(p[45]))")
		figure()
		PyPlot.plot(dts,z000,"xb")
		PyPlot.plot(dts,z025,"--b")
		PyPlot.plot(dts,z050, "b")
		PyPlot.plot(dts,z075,"--b")
		PyPlot.plot(dts,z100,"xb")
		xlabel("δt")
		ylabel("error (zero terms)")
		title("UK 45, error in a, (a=$(a[45]), f=$(f), phi=$(p[45]))")

# Error in phi
		e = Array{Float64,1}()
		
		for i in 1:length(dts)
			push!(e, eps[2][1][1])
		end
		figure()
		PyPlot.plot(dts,e,"-og")
		xlabel("δt")
		ylabel("error")
		title("UK 45, error in phi, (a=$(a[45]), f=$(f), phi=$(p[45]))")
	end

	return (Ahs, eAs, zAs), (fhs, efs), (ahs, eas, zas), (phs, eps)
end


function errors_uk_T(Ts::Array{Int64,1}, a::Array{Float64,1}, f::Float64, phi::Array{Float64,1}, dt::Float64, plot::Bool = false, eps::Float64 = 1e-5)
	n = 120	
	I = diagm(0 => ones(n))
	
	include("load_uk.jl")
	
	Mi = diagm(0 => 1 ./ m)
	D = diagm(0 => d)
		
	eAs = Array{Array{Float64,1},1}()
	zAs = Array{Array{Float64,1},1}()
	Ahs = Array{Array{Float64,2},1}()
	efs = Array{Float64,1}()
	fhs = Array{Float64,1}()
	eas = Array{Array{Float64,1},1}()
	zas = Array{Array{Float64,1},1}()
	ahs = Array{Array{Float64,1},1}()
	eps = Array{Array{Float64,1},1}()
	zps = Array{Array{Float64,1},1}()
	phs = Array{Array{Float64,1},1}()
	
	for T in Ts
		A = [I dt*I; (-dt*Mi*L) (I - dt*Mi*D)]
		
		Xs = readdlm("data2/uk_forced_$(f)_$(T)_$(dt).csv",',')
		
		Ah, fh = find_A_n_f(Xs, dt)
		
		ah, ph = locate_f_n_lag(Xs, dt, fh, Ah)
		
		eA = Array{Float64,1}()
		zA = Array{Float64,1}()
		for i in 1:2*n
			for j in 1:2*n
				if abs(A[i,j]) > 0
					push!(eA, abs(A[i,j] - Ah[i,j])/abs(A[i,j]))
				else
					push!(zA, abs(Ah[i,j]))
				end
			end
		end
		push!(eAs, eA)
		push!(zAs, zA)

		push!(Ahs, Ah)

		push!(efs, abs(f - fh)/abs(f))
		push!(fhs, fh)

		ea = Array{Float64,1}()
		za = Array{Float64,1}()
		for i in 1:n
			if abs(a[i]) > 0
				push!(ea, abs(a[i] - ah[i])/abs(a[i]))
			else
				push!(za, abs(ah[i]))
			end
		end
		push!(eas, ea)
		push!(zas, za)

		push!(ahs, ah)

		ep = Array{Float64,1}()
		for i in 1:n
			if a[i] > 0
				push!(ep, abs(p[i] - ph[i]))
			end
		end
		push!(eps, ep)

		push!(phs, ph)
	end
	
	if plot
# Error in A
		e000 = Array{Float64,1}()
		e025 = Array{Float64,1}()
		e050 = Array{Float64,1}()
		e075 = Array{Float64,1}()
		e100 = Array{Float64,1}()
		z000 = Array{Float64,1}()
		z025 = Array{Float64,1}()
		z050 = Array{Float64,1}()
		z075 = Array{Float64,1}()
		z100 = Array{Float64,1}()
		
		for i in 1:length(Ts)
			push!(e000,minimum(eAs[i]))
			push!(e025,quantile(eAs[i],.25))
			push!(e050,quantile(eAs[i],.5))
			push!(e075,quantile(eAs[i],.75))
			push!(e100,maximum(eAs[i]))
			push!(z000,minimum(zAs[i]))
			push!(z025,quantile(zAs[i],.25))
			push!(z050,quantile(zAs[i],.5))
			push!(z075,quantile(zAs[i],.75))
			push!(z100,maximum(zAs[i]))
		end
		figure()
		PyPlot.plot(Ts,e000,"xr")
		PyPlot.plot(Ts,e025,"--r")
		PyPlot.plot(Ts,e050, "r")
		PyPlot.plot(Ts,e075,"--r")
		PyPlot.plot(Ts,e100,"xr")
		xlabel("T")
		ylabel("relative error (non-zero terms)")
		title("UK 45, error in A, (a=$(a[45]), f=$(f), phi=$(p[45]))")
		figure()
		PyPlot.plot(Ts,z000,"xr")
		PyPlot.plot(Ts,z025,"--r")
		PyPlot.plot(Ts,z050, "r")
		PyPlot.plot(Ts,z075,"--r")
		PyPlot.plot(Ts,z100,"xr")
		xlabel("T")
		ylabel("error (zero terms)")
		title("UK 45, error in A, (a=$(a[45]), f=$(f), phi=$(p[45]))")
		
# Error in a
		e000 = Array{Float64,1}()
		e025 = Array{Float64,1}()
		e050 = Array{Float64,1}()
		e075 = Array{Float64,1}()
		e100 = Array{Float64,1}()
		z000 = Array{Float64,1}()
		z025 = Array{Float64,1}()
		z050 = Array{Float64,1}()
		z075 = Array{Float64,1}()
		z100 = Array{Float64,1}()
		
		for i in 1:length(Ts)
			push!(e000,minimum(eas[i]))
			push!(e025,quantile(eas[i],.25))
			push!(e050,quantile(eas[i],.5))
			push!(e075,quantile(eas[i],.75))
			push!(e100,maximum(eas[i]))
			push!(z000,minimum(zas[i]))
			push!(z025,quantile(zas[i],.25))
			push!(z050,quantile(zas[i],.5))
			push!(z075,quantile(zas[i],.75))
			push!(z100,maximum(zas[i]))
		end
		figure()
		PyPlot.plot(Ts,e000,"xb")
		PyPlot.plot(Ts,e025,"--b")
		PyPlot.plot(Ts,e050, "b")
		PyPlot.plot(Ts,e075,"--b")
		PyPlot.plot(Ts,e100,"xb")
		xlabel("T")
		ylabel("relative error (non-zero terms)")
		title("UK 45, error in a, (a=$(a[45]), f=$(f), phi=$(p[45]))")
		figure()
		PyPlot.plot(Ts,z000,"xb")
		PyPlot.plot(Ts,z025,"--b")
		PyPlot.plot(Ts,z050, "b")
		PyPlot.plot(Ts,z075,"--b")
		PyPlot.plot(Ts,z100,"xb")
		xlabel("T")
		ylabel("error (zero terms)")
		title("UK 45, error in a, (a=$(a[45]), f=$(f), phi=$(p[45]))")

# Error in phi
		e = Array{Float64,1}()
		
		for i in 1:length(Ts)
			push!(e, eps[2][1][1])
		end
		figure()
		PyPlot.plot(Ts,e,"-og")
		xlabel("T")
		ylabel("error")
		title("UK 45, error in phi, (a=$(a[45]), f=$(f), phi=$(p[45]))")
	end

	return (Ahs, eAs, zAs), (fhs, efs), (ahs, eas, zas), (phs, eps)
end


function errors_uk_Tdt(Tdts::Array{Tuple{Int64,Float64},1}, a::Array{Float64,1}, f::Float64, phi::Array{Float64,1}, plot::Bool = false, eps::Float64 = 1e-5)
	n = 120	
	I = diagm(0 => ones(n))
	
	include("load_uk.jl")
	
	Mi = diagm(0 => 1 ./ m)
	D = diagm(0 => d)
		
	eAs = Array{Array{Float64,1},1}()
	zAs = Array{Array{Float64,1},1}()
	Ahs = Array{Array{Float64,2},1}()
	efs = Array{Float64,1}()
	fhs = Array{Float64,1}()
	eas = Array{Array{Float64,1},1}()
	zas = Array{Array{Float64,1},1}()
	ahs = Array{Array{Float64,1},1}()
	eps = Array{Array{Float64,1},1}()
	zps = Array{Array{Float64,1},1}()
	phs = Array{Array{Float64,1},1}()
	
	for (T,dt) in Tdts
		A = [I dt*I; (-dt*Mi*L) (I - dt*Mi*D)]
		
		Xs = readdlm("data2/uk_forced_$(f)_$(T)_$(dt).csv",',')
		
		Ah, fh = find_A_n_f(Xs, dt)
		
		ah, ph = locate_f_n_lag(Xs, dt, fh, Ah)
		
		eA = Array{Float64,1}()
		zA = Array{Float64,1}()
		for i in 1:2*n
			for j in 1:2*n
				if abs(A[i,j]) > 0
					push!(eA, abs(A[i,j] - Ah[i,j])/abs(A[i,j]))
				else
					push!(zA, abs(Ah[i,j]))
				end
			end
		end
		push!(eAs, eA)
		push!(zAs, zA)

		push!(Ahs, Ah)

		push!(efs, abs(f - fh)/abs(f))
		push!(fhs, fh)

		ea = Array{Float64,1}()
		za = Array{Float64,1}()
		for i in 1:n
			if abs(a[i]) > 0
				push!(ea, abs(a[i] - ah[i])/abs(a[i]))
			else
				push!(za, abs(ah[i]))
			end
		end
		push!(eas, ea)
		push!(zas, za)

		push!(ahs, ah)

		ep = Array{Float64,1}()
		for i in 1:n
			if a[i] > 0
				push!(ep, abs(p[i] - ph[i]))
			end
		end
		push!(eps, ep)

		push!(phs, ph)
	end
	
	if plot
		Ts = Array{Int64,1}()
		for (T,dt) in Tdts
			push!(Ts,T)
		end
# Error in A
		e000 = Array{Float64,1}()
		e025 = Array{Float64,1}()
		e050 = Array{Float64,1}()
		e075 = Array{Float64,1}()
		e100 = Array{Float64,1}()
		z000 = Array{Float64,1}()
		z025 = Array{Float64,1}()
		z050 = Array{Float64,1}()
		z075 = Array{Float64,1}()
		z100 = Array{Float64,1}()
		
		for i in 1:length(Ts)
			push!(e000,minimum(eAs[i]))
			push!(e025,quantile(eAs[i],.25))
			push!(e050,quantile(eAs[i],.5))
			push!(e075,quantile(eAs[i],.75))
			push!(e100,maximum(eAs[i]))
			push!(z000,minimum(zAs[i]))
			push!(z025,quantile(zAs[i],.25))
			push!(z050,quantile(zAs[i],.5))
			push!(z075,quantile(zAs[i],.75))
			push!(z100,maximum(zAs[i]))
		end
		figure()
		PyPlot.plot(Ts,e000,"xr")
		PyPlot.plot(Ts,e025,"--r")
		PyPlot.plot(Ts,e050, "r")
		PyPlot.plot(Ts,e075,"--r")
		PyPlot.plot(Ts,e100,"xr")
		xlabel("T")
		ylabel("relative error (non-zero terms)")
		title("UK 45, error in A (constant T*δt), (a=$(a[45]), f=$(f), phi=$(p[45]))")
		figure()
		PyPlot.plot(Ts,z000,"xr")
		PyPlot.plot(Ts,z025,"--r")
		PyPlot.plot(Ts,z050, "r")
		PyPlot.plot(Ts,z075,"--r")
		PyPlot.plot(Ts,z100,"xr")
		xlabel("T")
		ylabel("error (zero terms)")
		title("UK 45, error in A (constant T*δt), (a=$(a[45]), f=$(f), phi=$(p[45]))")
		
# Error in a
		e000 = Array{Float64,1}()
		e025 = Array{Float64,1}()
		e050 = Array{Float64,1}()
		e075 = Array{Float64,1}()
		e100 = Array{Float64,1}()
		z000 = Array{Float64,1}()
		z025 = Array{Float64,1}()
		z050 = Array{Float64,1}()
		z075 = Array{Float64,1}()
		z100 = Array{Float64,1}()
		
		for i in 1:length(Ts)
			push!(e000,minimum(eas[i]))
			push!(e025,quantile(eas[i],.25))
			push!(e050,quantile(eas[i],.5))
			push!(e075,quantile(eas[i],.75))
			push!(e100,maximum(eas[i]))
			push!(z000,minimum(zas[i]))
			push!(z025,quantile(zas[i],.25))
			push!(z050,quantile(zas[i],.5))
			push!(z075,quantile(zas[i],.75))
			push!(z100,maximum(zas[i]))
		end
		figure()
		PyPlot.plot(Ts,e000,"xb")
		PyPlot.plot(Ts,e025,"--b")
		PyPlot.plot(Ts,e050, "b")
		PyPlot.plot(Ts,e075,"--b")
		PyPlot.plot(Ts,e100,"xb")
		xlabel("T")
		ylabel("relative error (non-zero terms)")
		title("UK 45, error in a (constant T*δt), (a=$(a[45]), f=$(f), phi=$(p[45]))")
		figure()
		PyPlot.plot(Ts,z000,"xb")
		PyPlot.plot(Ts,z025,"--b")
		PyPlot.plot(Ts,z050, "b")
		PyPlot.plot(Ts,z075,"--b")
		PyPlot.plot(Ts,z100,"xb")
		xlabel("T")
		ylabel("error (zero terms)")
		title("UK 45, error in a (constant T*δt), (a=$(a[45]), f=$(f), phi=$(p[45]))")

# Error in phi
		e = Array{Float64,1}()
		
		for i in 1:length(Ts)
			push!(e, eps[2][1][1])
		end
		figure()
		PyPlot.plot(Ts,e,"-og")
		xlabel("T")
		ylabel("error")
		title("UK 45, error in phi (constant T*δt), (a=$(a[45]), f=$(f), phi=$(p[45]))")
	end

	return (Ahs, eAs, zAs), (fhs, efs), (ahs, eas, zas), (phs, eps)
end






