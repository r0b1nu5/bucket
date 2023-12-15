using EDF, DelimitedFiles, FFTW

include("this.jl")

function read_eeg(file::String)
	xxx = EDF.read(file)

	s2signal = Dict{String,Vector{Float64}}()

	for sig in xxx.signals
		if typeof(sig) == EDF.Signal{Int16}
			s2signal[sig.header.label] = EDF.decode(sig)
		end
	end

	return s2signal
end

function average_over_zones(s2signal::Dict{String,Vector{Float64}}, s2z::Dict{String,Int64})
	l = maximum(values(s2z))
	T = length(s2signal[collect(keys(s2signal))[1]])

	c = zeros(l)
	asig = zeros(l,T)
	for k in keys(s2signal)
		z = s2z[k]
		c[z] += 1
		asig[z,:] = asig[z,:]*(c[z]-1)/c[z] + s2signal[k]/c[z]
	end

	return asig
end

function load_and_save_eeg_avg(subjects::Vector{String}, states::Vector{String}, nz::Int64)
	s = readdlm("eeg-data/sensors-$nz.csv",',',String)
	z = readdlm("eeg-data/zones-$nz.csv",',',Int64)
	s2z = Dict{String,Int64}(s[i] => z[i] for i in 1:length(s))
	for su in subjects
		for st in states
			@info "Working on S"*su*"R"*st
			s2signal = read_eeg("eeg-data/S"*su*"R"*st*".edf")
			asig = average_over_zones(s2signal,s2z)
			writedlm("eeg-data/S"*su*"R"*st*"-X.csv",asig,',')
		end
	end

	return nothing
end

function list_all_subjects(n::Int64)
	subjects = String[]
	for i in 1:min(9,n)
		push!(subjects,"00"*string(i))
	end
	for i in 10:min(99,n)
		push!(subjects,"0"*string(i))
	end
	for i in 100:min(999,n)
		push!(subjects,string(i))
	end

	return subjects
end

function plot_brain_3dge(ids::Vector{Int64}, p::Float64)
	plot_brain_3dge(ids,p,[(0.,0.,0.,1.) for i in 1:7])
end

function plot_brain_3dge(ids::Vector{Int64}, p::Float64, cols=Vector{NTuple{4,Float64}})
#	brain_bord = readdlm("eeg-data/brain-bord.csv",',')
	xy = readdlm("eeg-data/zones-7-xy.csv",',')

#	PyPlot.plot(brain_bord[:,1],brain_bord[:,2],"k")
	x = LinRange(-1,1,300)
	PyPlot.plot(x,-.8*sqrt.(1 .- x.^2),"k",x,.8*sqrt.(1 .- x.^2),"k")
	for i in 1:7
		PyPlot.plot(xy[i,1],xy[i,2],"o",color=cols[i],markersize=10)
	end
	PyPlot.fill(xy[ids,1],xy[ids,2],color="#999999")
	PyPlot.plot(xy[ids[1],1],xy[ids[1],2],"o",color=cols[ids[1]],markersize=20)
	PyPlot.plot(xy[ids[2],1],xy[ids[2],2],"o",color=cols[ids[2]],markersize=15)
	PyPlot.plot(xy[ids[3],1],xy[ids[3],2],"o",color=cols[ids[3]],markersize=10)
	title("p = $(round(p,digits=2)*100)%")
	xticks([])
	yticks([])
end

function vectorize_adj(A::Matrix{Float64})
	n = size(A)[1]
	a = Float64[]
	for i in 1:n-1
		for j in i+1:n
			push!(a,A[i,j])
			push!(a,A[j,i])
		end
	end
	return a
end

function vectorize_adj(A::Array{Float64,3})
	n = size(A)[1]
	a = Float64[]
	for i in 1:n-2
		for j in i+1:n-1
			for k in j+1:n
				push!(a,A[i,j,k])
				push!(a,A[i,k,j])
				push!(a,A[j,i,k])
				push!(a,A[j,k,i])
				push!(a,A[k,i,j])
				push!(a,A[k,j,i])
			end
		end
	end
	return a
end

# Use a sliding mean with Δt horizon
function smooth_time_series(x::Vector{Float64}, Δt::Int64=10)
	return [mean(x[t:t+Δt]) for t in 1:(length(x)-Δt)]
end

function smooth_time_series(x::Matrix{Float64}, Δt::Int64=10)
	n,T = size(x)
	y = zeros(T-Δt,0)
	for i in 1:n
		y = [y smooth_time_series(x[i,:],Δt)]
	end
	return Matrix(y')
end

function subsample_time_series(x::Vector{Float64}, Δt::Int64=10)
	return [mean(x[t:t+Δt]) for t in 1:Δt:length(x)-Δt]
end

function subsample_time_series(x::Matrix{Float64}, Δt::Int64=10)
	n,T = size(x)
	y = zeros(Int64(floor(T/Δt)),0)
	for i in 1:n
		y = [y subsample_time_series(x[i,:],Δt)]
	end
	return Matrix(y')
end

# Estimates the derivative from discrete data points with the 9-point stencil.
# Keeps data point every Δt time steps.
# Time step length is δt.
function stencil9(x::Vector{Float64}, δt::Float64, Δt::Int64=10)
	c = [1/280,-4/105,1/5,-4/5,0,4/5,-1/5,4/105,-1/280]./δt
	return x[5:Δt:length(x)-4], [dot(c,x[(-4:4) .+ t]) for t in 5:Δt:length(x)-4]
end

function stencil9(X::Matrix{Float64}, δt::Float64, Δt::Int64=10)
	n,T = size(X)
	x = zeros(length(5:Δt:T-4),0)
	dx = zeros(length(5:Δt:T-4),0)
	for i in 1:n
		a,b = stencil9(X[i,:],δt,Δt)
		x = [x a]
		dx = [dx b]
	end
	return Matrix(x'), Matrix(dx')
end

function denoise_fourier(x::Vector{Float64},nmodes::Int64)
	T = length(x)
	f = fft(x)
	g = zeros(Complex{Float64},T)
	g[1:nmodes+1] = f[1:nmodes+1]
	g[T-nmodes+1:T] = f[T-nmodes+1:T]
	return real.(ifft(g))
end

function denoise_fourier(X::Matrix{Float64},nmodes::Int64)
	n,T = size(X)
	Y = zeros(n,T)
	for i in 1:n
		Y[i,:] = denoise_fourier(X[i,:],nmodes)
	end
	return Y
end

function compare_eeg_trajectories(X::Matrix{Float64}, Y::Matrix{Float64}, Ξ::Matrix{Float64}, dmax::Int64)
	n,T = size(X)
	
	θ = get_θ(X,dmax)
	Yh = Ξ*θ 

	for i in 1:min(n,10)
		PyPlot.plot(Y[i,:],color="C$(i-1)")
		PyPlot.plot(Yh[i,:],"--",color="C$(i-1)")
	end

	return sum((Y-Yh).^2)
end

function cleaned_time_series(su::String, st::String, truncat::Int64)
	nz = 7
	s = readdlm("eeg-data/sensors-$nz.csv",',',String)
	z = readdlm("eeg-data/zones-$nz.csv",',',Int64)
	s2z = Dict{String,Int64}(s[i] => z[i] for i in 1:length(s))

	s2signal = read_eeg("eeg-data/S"*su*"R"*st*".edf")
	asig = average_over_zones(s2signal,s2z)

	dt = 1/160
	X0 = denoise_fourier(asig[:,1:truncat],100)
	Y0 = (X0[:,2:end] - X0[:,1:end-1])/dt
	X0 = X0[:,1:end-1]
	X = X0./mean(abs.(X0))
	Y = Y0./mean(abs.(Y0))

	return X,Y
end

function restrict_box_size(X::Matrix{Float64}, nstep::Int64)
	n,T = size(X)
	ns = min(T-1,nstep)

	mX = median(X,dims=2)
	dX = X - repeat(mX,1,T)
	nX = vec(sum(dX.^2,dims=1))
	sX = sort(nX)
	th = (sX[ns+1] + sX[ns])/2
	
	ids = (1:T)[nX .< th]

	return X[:,ids], ids
end




