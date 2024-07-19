using PyPlot, DelimitedFiles, Statistics, LinearAlgebra

include("eeg-tools.jl")

n = 64
type = "full" 
data_exist = false
Tmax = 1000
thr2 = 1.
thr3 = 1.
thr4 = 1.

do_main_edges = false
do_violin = false
do_relative_error = false
do_hypergraph_density = false
do_trajectory = true

# Data to be loaded
subjects = list_all_subjects(109)
#subjects = ["001","002"]
states = ["01","02"]
#states = ["03","07","11"]

suffix = "222"

ρ3 = zeros(n,0)
relerr = zeros(length(states),0)
den2 = Float64[]
den3 = Float64[]

if do_trajectory
	su = "015"
	st = "01"
	i0 = 7

	s2signal = read_eeg("eeg-data/S"*su*"R"*st*".edf")
	sig = zeros(0,length(s2signal["Af3."]))
	for s in sort(collect(keys(s2signal)))
		global sig = [sig; s2signal[s]']
	end
		
	truncat = findmin(vec(maximum(abs.([sig zeros(n)]),dims=1)) .> 1e-6)[2] - 1
	
	dt = 1/160 

	X0 = denoise_fourier(sig[:,1:truncat],200)
	Y0 = (X0[:,2:end]-X0[:,1:end-1])./dt

	X0 = X0[:,1:end-1]
	X = X0./mean(abs.(X0))
	Y = Y0./mean(abs.(Y0))

	X,ids = restrict_box_size(X,1000)
	Y = Y[:,ids]

	y0 = readdlm("eeg-data/S"*su*"R"*st*"-"*type*"-coeff-"*suffix*".csv",',')[i0,1]
	A2 = readdlm("eeg-data/S"*su*"R"*st*"-"*type*"-A2-"*suffix*".csv",',')
	A3 = readdlm("eeg-data/S"*su*"R"*st*"-"*type*"-A3-"*suffix*".csv",',')

#	Ythis = zeros(1000)
	Ythis = y0*ones(1000)
	for l in 1:size(A2)[1]
		i,j = Int64.(A2[l,1:2])
		a = A2[l,3]
		if i == i0
			global Ythis += a*X[j,:]
		end
	end
	for l in 1:size(A3)[1]
		i,j,k = Int64.(A3[l,1:3])
		a = A3[l,4]
		if i == i0
			global Ythis += a*X[j,:].*X[k,:]
		end
	end

	re = round(readdlm("eeg-data/S"*su*"R"*st*"-"*type*"-re-"*suffix*".csv",',')[1],digits=2)

	figure("Fit",(13,4))
	PyPlot.plot(Y[i0,:])
	PyPlot.plot(Ythis,"--")
	xlabel("time steps (restricted to the box)")
	ylabel("Derivative of sensor $i0")
	title("Fitting the derivative of sensor $i0, for subject S"*su*"R"*st*" (relative error: $re)")
end


if do_hypergraph_density
	zer0 = 1e-6
	for su in subjects
		@info "S"*su*"Rxx"
		for st in states
			A2 = readdlm("eeg-data/S"*su*"R"*st*"-"*type*"-A2-"*suffix*".csv",',')
			A3 = readdlm("eeg-data/S"*su*"R"*st*"-"*type*"-A3-"*suffix*".csv",',')
			push!(den2,sum(A2[:,3] .> zer0)/(n*binomial(n-1,1)))
			push!(den3,sum(A3[:,4] .> zer0)/(n*binomial(n-1,2)))
		end
	end
	figure()
	ax = subplot(2,1,1)
	PyPlot.hist(den2,30)
	PyPlot.plot([0;sort(den2)],(0:218)*100/218)
	ylabel("occurences (out of 218)")
	secax = ax.secondary_yaxis("right",functions=(x->x/100,x->x*100))
	secax.set_ylabel("cumulative")
	title("Histogram: hypergraph density (top: pairwise, bottom: triadic)")
	ax = subplot(2,1,2)
	PyPlot.hist(den3,30)
	PyPlot.plot([0;sort(den3)],(0:218)*100/218)
	xlabel("density")
	ylabel("occurences (out of 218)")
	secax = ax.secondary_yaxis("right",functions=(x->x/100,x->x*100))
	secax.set_ylabel("cumulative")
end


if do_relative_error
	for su in subjects
		re = Float64[]
		for st in states
			push!(re,readdlm("eeg-data/S"*su*"R"*st*"-"*type*"-re-"*suffix*".csv",',')[1])
		end
		global relerr = [relerr re]
	end

	fig,ax = PyPlot.subplots()
	PyPlot.hist(vec(relerr),30)
	PyPlot.plot([0;sort(vec(relerr))],56*(0:218)./218,"-",lw=2)
	xlabel("Relative error")
	ylabel("Frequency")
	secax = ax.secondary_yaxis("right", functions=(x -> x/56, x -> 56*x))
	secax.set_ylabel("Cumulative")
	title("EEG 64 sensors: relative errors")
end


if do_main_edges
	n_edges = 10 # Number of most prominent edges to keep
	
	edge_score2 = Dict{Vector{Int64},Float64}()
	edge_score3 = Dict{Vector{Int64},Float64}()
	
	for su in subjects
		for st in states
			a2 = readdlm("eeg-data/S"*su*"R"*st*"-"*type*"-A2-"*suffix*".csv",',')
			t2 = norm(a2[:,3],1)
			a3 = readdlm("eeg-data/S"*su*"R"*st*"-"*type*"-A3-"*suffix*".csv",',')
			t3 = norm(a3[:,4],1)
		
			for l in 1:size(a2)[1]
				e = a2[l,1:2]
				a = a2[l,3]
				if haskey(edge_score2,e)
					edge_score2[e] += abs(a)/t2
				else
					edge_score2[e] = abs(a)/t2
				end
			end
	
			for l in 1:size(a3)[1]
				e = a3[l,1:3]
				a = a3[l,4]
				if haskey(edge_score3,e)
					edge_score3[e] += abs(a)/t3
				else
					edge_score3[e] = abs(a)/t3
				end
			end
		end
	end
	
	E2 = Vector{Int64}[]
	S2 = Float64[]
	for i in 1:length(edge_score2)
		sc,id = findmax(edge_score2)
		push!(E2,id)
		push!(S2,sc)
		pop!(edge_score2,id)
	end
	
	E3 = Vector{Int64}[]
	S3 = Float64[]
	for i in 1:length(edge_score3)
		sc,id = findmax(edge_score3)
		push!(E3,id)
		push!(S3,sc)
		pop!(edge_score3,id)
	end
	
	
	@info ""
	@info ""
	@info "The $n_edges most prominent k-edges are:"
	@info ""
	@info "| k |   2    |    3      |"
	@info "--------------------------"
	for i in 1:n_edges
		@info "|   | $(E2[i]) | $(E3[i]) |"
	end
	
	#@info ""; @info ""; @info "The $n_edges most prominent k-edges are:"; @info ""; @info "| k |   2    |    3      |     4        |"; @info "-----------------------------------------"; for i in 1:n_edges; @info "|   | $(E2[i]) | $(E3[i]) | $(E4[i]) |"; end
end



 
if do_violin
	c = 0
	for su in subjects
		for st in states
			global c += 1
			@info "Working on S"*su*"R"*st*": run "*suffix
			if !(data_exist)
				file = "eeg-data/S"*su*"R"*st*".edf"
				s2signal = read_eeg(file)
				sig = zeros(0,length(s2signal["Af3."]))
				for s in sort(collect(keys(s2signal)))
					sig = [sig; s2signal[s]']
				end
			
				# Finite differences
				dt = 1/160
					truncat = findmin(vec(maximum(abs.([sig zeros(n)]),dims=1)) .> 1e-6)[2] - 1
			
				X0 = denoise_fourier(sig[:,1:truncat],200)
				X0 = X0[:,1:end-1]
				X = X0./mean(abs.(X0))

				X,ids = restrict_box_size(X,1000)
				T = size(X)[2]
			end
	
			a2 = readdlm("eeg-data/S"*su*"R"*st*"-"*type*"-A2-"*suffix*".csv",',')
			t2 = norm(a2[:,3],1)
			a3 = readdlm("eeg-data/S"*su*"R"*st*"-"*type*"-A3-"*suffix*".csv",',')
			t3 = norm(a3[:,4],1)
	
			if !(data_exist)
				writedlm("eeg-data/T-S"*su*"R"*st*"-"*suffix*".csv",T,',')
				for t in 1:T
					z2 = zeros(n)
					for l in 1:size(a2)[1]
						i,j = Int64.(a2[l,1:2])
						e = [i,j]
						a = a2[l,3]
						z2[i] += abs(a*X[j,t])
						if haskey(edge_score2,e)
							edge_score2[e] += abs(a)/t2
						else
							edge_score2[e] = abs(a)/t2
						end
					end
					z3 = zeros(n)
					for l in 1:size(a3)[1]
						i,j,k = Int64.(a3[l,1:3])
						e = [i,j,k]
						a = a3[l,4]
						z3[i] += abs(a*X[j,t]*X[k,t])
						if haskey(edge_score3,e)
							edge_score3[e] += abs(a)/t3
						else
							edge_score3[e] = abs(a)/t3
						end
					end
					z = z2+z3
					ρ3 = [z[i] == 0. ? 0 : z3[i]/z[i] for i in 1:n]
					writedlm("eeg-data/violin-data3-"*type*"-$n-S"*su*"R"*st*"-"*suffix*"-$t.csv",ρ3,',')
				end
				ρ3 = zeros(n,0)
				for t in 1:T
					ρ3 = [ρ3 vec(readdlm("eeg-data/violin-data3-"*type*"-$n-S"*su*"R"*st*"-"*suffix*"-$t.csv",','))]
					rm("eeg-data/violin-data3-"*type*"-$n-S"*su*"R"*st*"-"*suffix*"-$t.csv")
					if t%1000 == 0
						writedlm("eeg-data/violin-data3-"*type*"-$n-S"*su*"R"*st*"-"*suffix*"-T$t.csv",ρ3,',')
						ρ3 = zeros(n,0)
					end
				end
			end
		end
	end

	ρ3 = zeros(n,0)
	for su in subjects
		for st in states
			global ρ3 = [ρ3 readdlm("eeg-data/violin-data3-"*type*"-$n-S"*su*"R"*st*"-"*suffix*"-T1000.csv",',')]
		end
	end
	
	contr_per_subject = zeros(n,0)
	q0 = zeros(n,0)
	q1 = zeros(n,0)
	q2 = zeros(n,0)
	q3 = zeros(n,0)
	q4 = zeros(n,0)
	nams = String[]
	for su in subjects
		for st in states
			@info "Loading S"*su*"R"*st
	
			T = Int64(1000*floor(readdlm("eeg-data/T-S"*su*"R"*st*"-"*suffix*".csv",',')[1]/1000))
			ρ = zeros(n,0)
			for t in 1000:1000:min(T,Tmax)
				ρ = [ρ readdlm("eeg-data/violin-data3-"*type*"-$n-S"*su*"R"*st*"-"*suffix*"-T$t.csv",',')]
			end
			global ρ3 = [ρ3 ρ]
			global contr_per_subject = [contr_per_subject mean(ρ,dims=2)]
			global q0 = [q0 [quantile(ρ[i,:],.00) for i in 1:n]]
			global q1 = [q1 [quantile(ρ[i,:],.25) for i in 1:n]]
			global q2 = [q2 [quantile(ρ[i,:],.50) for i in 1:n]]
			global q3 = [q3 [quantile(ρ[i,:],.75) for i in 1:n]]
			global q4 = [q4 [quantile(ρ[i,:],1.0) for i in 1:n]]
			push!(nams,"S"*su*"R"*st)
		end
	end
	
	figure("Violins-ter",(10,4))
	ζ = [vec(ρ3),]
	qs = quantile(vec(q2),[.05,.35,.65,.95])
	#qs = quantile(vec(q2),[.2,.4,.6,.8])
	m1,i1 = findmin(abs.(q2 .- qs[1]))
	m2,i2 = findmin(abs.(q2 .- qs[2]))
	m3,i3 = findmin(abs.(q2 .- qs[3]))
	m4,i4 = findmin(abs.(q2 .- qs[4]))
	for i in [i1,i2,i3,i4]
		push!(ζ,ρ3[i[1],(i[2]-1)*Tmax .+ (1:Tmax)])
	end
	fig = plt.violinplot(ζ,showextrema=false)
	cmap = get_cmap("magma")
	cols = [cmap(r) for r in [.1,.5,.6,.7,.8,.9]]
	c = 0
	for pc in fig["bodies"]
		global c += 1
		pc.set_facecolor(cols[c])
		pc.set_edgecolor("black")
		pc.set_alpha(.8)
	end
	cc = ["gray","black","black","black","black","black"]
	for i in 1:length(ζ)
		PyPlot.plot([i,i],[minimum(ζ[i]),maximum(ζ[i])],color=cc[i])
		PyPlot.plot([i,i],quantile(ζ[i],[.25,.75]),color=cc[i],linewidth=5)
		PyPlot.plot(i,median(ζ[i]),"o",color=cc[i],markersize=10)
	end
	#xticks([1,2,3,4,5,6],["All","1%","25%","50%","75%","99%"])
	#xticks([1,2,3,4,5,6],["All","0%","25%","50%","75%","100%"])
	#xticks([1,2,3,4,5],["All","20%","40%","60%","80%"])
	xticks([1,2,3,4,5],["All","5%","35%","65%","95%"])
	xlabel("Percentile")
	ylabel("Amount of the dynamics that is\nexplained by 3rd-order interactions")
end



