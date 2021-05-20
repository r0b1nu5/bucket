using PyPlot, DelimitedFiles

include("LoadTsatTxt.jl")
include("data_naspi/osc_init.jl")
include("data_naspi/naspi_files.jl")

function rm_nan(Xs::Array{Float64,2})
	n,T = size(Xs)

	Xt = Array{Float64,2}(undef,0,T)

	for i in 1:n
		x = Xs[i,:]
		xnan = isnan.(x)

		b_nan = Array{Int64,1}()
		e_nan = Array{Int64,1}()
		now_nan = false

		if isnan(x[1])
			now_nan = true
			push!(b_nan,1)
		end

		for j in 2:length(x)
			if now_nan && xnan[j-1:j] == [1,0]
				push!(e_nan,j-1)
				now_nan = false
			elseif !now_nan && xnan[j-1:j] == [0,1]
				push!(b_nan,j)
				now_nan = true
			end
		end

		if now_nan
			push!(e_nan,T)
		end

		for j in 1:length(b_nan)
			b = b_nan[j]
			e = e_nan[j]

			if b == 1
				x[b:e] = x[e+1]*ones(e-b+1)
			elseif e == T
				x[b:e] = x[b-1]*ones(e-b+1)
			else
				x[b:e] = (x[e+1]-x[b-1])/((e+1)-(b-1))*((b:e) .- (b-1)) .+ x[b-1]
			end
		end

		Xt = [Xt;x']
	end

	return Xt
end

# Estimates the derivatives of the time series x, with step size h, using the "five point stencil" method.

function stencil_5p(x::Array{Float64,1}, h::Float64)
	T = length(x)

	dx = Array{Float64,1}()
	for i in 3:T-2
		push!(dx,(x[i-2] - 8*x[i-1] + 8*x[i+1] - x[i+2])/(12*h))
	end

	return x[3:T-2],dx
end

function stencil_5p(X::Array{Float64,2}, h::Float64)
	n,T = size(X)

	dX = Array{Float64,2}(undef,0,T-4)
	for i in 1:n
		x,dx = stencil_5p(X[i,:],h)
		dX = [dX;dx']
	end

	return X[:,3:T-2],dX
end


function prepare_naspi(case::String, file::String, init::Int64)
	t,Xtt,tit = LoadTsatTxt(file*case*"/BusVolAng.txt")
	Xt = rm_nan(Xtt)

	τ = t[2]-t[1]
	
	X,dX = stencil_5p(Xt[:,init:end]*π/180,τ)

	Xs = [X;dX]

	writedlm("data_naspi/naspi_ids_"*case*".csv",tit,',')
	writedlm("data_naspi/Xs_"*case*".csv",Xs,',')

	return Xs
end
