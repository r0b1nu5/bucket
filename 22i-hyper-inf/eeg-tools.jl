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

function vectorize_adj(A::Array{Float64,4})
	n = size(A)[1]
	a = Float64[]
	for i in 1:n-3
		for j in i+1:n-2
			for k in j+1:n-1
				for l in k+1:n
					for p in permutations([i,j,k,l])
						push!(a,A[p[1],p[2],p[3],p[4]])
					end
				end
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


function contour_trigon_data(x::Vector{Float64}, y::Vector{Float64}, res::Int64=100, fig::String="trigon", scale::String="lin", cm::String="plasma", frame::Bool=true)
	N = length(x)

	M = [1 cos(π/3);0 sin(π/3)]
	Mi = inv(M)
	T = zeros(res,res)

 #=
	for i in 1:res
		for j in 1:res
			z = Mi*[i,j]
			if minimum(z) < 0 || sum(z) > res
				T[j,i] = NaN
			end
		end
	end		
# =#

	for i in 1:N
		#z = ceil.(Int64,res*M*[x[i],y[i]])
		z = floor.(Int64,1 .+ (res-1)*M*[x[i],y[i]])
		T[z[2],z[1]] += 1
	end

	figure(fig,(6,5))
	cmap = get_cmap(cm)
	if scale == "lin"
		PyPlot.contourf(LinRange(0,1,res),LinRange(0,1,res),T,cmap=cmap,res)
	elseif scale == "log"
		PyPlot.contourf(LinRange(0,1,res),LinRange(0,1,res),log.(10,T .+ .1),cmap=cmap,res)
	end
	
	if frame
		plot_trigon_frame(fig)
		if scale == "lin"
			colorbar(label="# occurences")
		elseif scale == "log"
			colorbar(label="log(# occurences)",ticks=(-1:round(Int64,log(10,maximum(T)))))
		end
	end

	return T
end

function contour_trigon_data(x::Matrix{Float64}, y::Matrix{Float64}, res::Int64=100, fig::String="trigon", scale::String="lin", cm::String="plasma", frame::Bool=true)
	contour_trigon_data(vec(x),vec(y),res,fig,scale,cm,frame)

	return nothing
end


function plot_trigon_data(x::Float64, y::Float64, fig::String="trigon", col::String="C0", al::Float64=1., frame::Bool=false)
	z = [1 cos(π/3);0 sin(π/3)]*[x,y]
	figure(fig)
	PyPlot.plot(z[1],z[2],".",color=col,alpha=al)
	if frame
		plot_trigon_frame(fig)
	end

	return nothing
end

function plot_trigon_data(x::Vector{Float64}, y::Vector{Float64}, fig::String="trigon", col::String="C0", al::Float64=1., frame::Bool=true)
	for i in 1:length(x)
		plot_trigon_data(x[i],y[i],fig,col,al,false)
	end
	if frame
		plot_trigon_frame(fig)
	end

	return nothing
end

function plot_trigon_data(x::Matrix{Float64}, y::Matrix{Float64}, fig::String="trigon", col::String="C0", al::Float64=1., frame::Bool=true)
	plot_trigon_data(vec(x),vec(y),fig,col,al,frame)

	return nothing
end

function plot_trigon_frame(fig::String="trigon")
	figure(fig)
	PyPlot.fill([-.01,cos(π/3),1.01,1.01,-.01,-.01],[-.01,sin(π/3),-.01,1.01,1.01,-.01],"w")
	PyPlot.plot([0,1,cos(π/3),0],[0,0,sin(π/3),0],"-k")
	axis([-.05,1.05,-.11,1-.01])
	xticks([])
	yticks([])

	return nothing
end

function plot_trigon_label(x::String="x", y::String="y", z::String="z")
	for pos in [.25,.5,.75]
		PyPlot.plot([pos,pos],
			    [0,-.02],
			    "-k",lw=1,)
		PyPlot.text(pos,-.05,"$pos",size="small",ha="center",va="center")

		PyPlot.plot([1+pos*cos(2π/3),1+pos*cos(2π/3)+.02*cos(π/6)],
			    [0+pos*sin(2π/3),0+pos*sin(2π/3)+.02*sin(π/6)],
			    "-k",lw=1)
		PyPlot.text(1+pos*cos(2π/3)+.05*cos(π/6),0+pos*sin(2π/3)+.05*sin(π/6),"$pos",size="small",rotation=-60.,rotation_mode="anchor",ha="center",va="center")

		PyPlot.plot([cos(π/3)+pos*cos(4π/3),cos(π/3)+pos*cos(4π/3)+.02*cos(5π/6)],
			    [sin(π/3)+pos*sin(4π/3),sin(π/3)+pos*sin(4π/3)+.02*sin(5π/6)],
			    "-k",lw=1)
		PyPlot.text(cos(π/3)+pos*cos(4π/3)+.05*cos(5π/6),sin(π/3)+pos*sin(4π/3)+.05*sin(5π/6),"$pos",size="small",rotation=60.,rotation_mode="anchor",ha="center",va="center")
	end

	PyPlot.text(.5,-.15,x,ha="center",va="center")
	PyPlot.text(1+.5*cos(2π/3)+.15*cos(π/6),0+.5*sin(2π/3)+.15*sin(π/6),y,rotation=-60.,rotation_mode="anchor",ha="center",va="center")
	PyPlot.text(cos(π/3)+.5*cos(4π/3)+.15*cos(5π/6),sin(π/3)+.5*sin(4π/3)+.15*sin(5π/6),z,rotation=60.,rotation_mode="anchor",ha="center",va="center")

	return nothing
end


function plot_relerr(col::Any="C0")
	E = zeros(218,0)

	for o in 2:7
		E = [E vec(readdlm("eeg-data/relative-error-$(o)$(o)x.csv",','))]
	end

	q00 = vec(minimum(E,dims=1))
	q10 = [quantile(E[:,o-1],.1) for o in 2:7]
	q25 = [quantile(E[:,o-1],.25) for o in 2:7]
	q50 = vec(median(E,dims=1))
	q75 = [quantile(E[:,o-1],.75) for o in 2:7]
	q90 = [quantile(E[:,o-1],.9) for o in 2:7]
	q100 = vec(maximum(E,dims=1))
	qm = vec(mean(E,dims=1))
	
	figure("EEG",(6,4))
	#PyPlot.plot(2:7,q00,"x",color=col)
	#PyPlot.plot(2:7,q100,"x",color=col)
	PyPlot.fill([2:7;7:-1:2],[q10;q90[6:-1:1]],color=col,alpha=.2,label="10%-90% quantiles")
	PyPlot.fill([2:7;7:-1:2],[q25;q75[6:-1:1]],color=col,alpha=.7,label="25%-75% quantiles")
	PyPlot.plot(2:7,q50,color=col,lw=2.,label="median")
	PyPlot.plot(2:7,qm,"--k",lw=2.,label="mean")

	xlabel("interaction order")
	ylabel("relative error")
	legend()
end


function get_d(n::Int64, dmax::Int64,i0::Int64=1)
	d = zeros(Int64,1,dmax)
	if dmax == 0
		return d 
	else
		for i in 1:n
			d0 = get_d(n-i+1,dmax-1,i)
			d = [d; d0 i*ones(Int64,size(d0)[1],1)]
		end
		d += (i0-1)*(d .> 0)
		
		return d
	end
end

function get_idx_mon(n::Int64, dmax::Int64)
	d = get_d(n,dmax)

	idx_mon = Dict{Int64,Vector{Int64}}()
	for i in 1:size(d)[1]
		mon = d[i,:][d[i,:] .!= 0]
		if length(mon) == length(union(mon))
			idx_mon[i] = sort(mon)
		end
	end

	return idx_mon
end




