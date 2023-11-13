using EDF, DelimitedFiles

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





