using PyPlot, DelimitedFiles, Statistics, LinearAlgebra

include("eeg-tools.jl")

n = 64
dmax = 3
type = "full" 
data_exist = true
Tmax = 1000
thr2 = 1.
thr3 = 1.
thr4 = 1.

do_main_edges = false
do_violin = false
do_violin_4 = false
do_triangle_full = true
do_triangle_part = false
do_relative_error = false
do_hypergraph_density = false
do_trajectory = false
do_list_bad_inference = false

# Data to be loaded
subjects = list_all_subjects(109)
#subjects = ["001","002"]
states = ["01","02"]
#states = ["03","07","11"]
#list_names = ["S"*su*"R"*st for su in subjects for st in states]
#list_names = vec(readdlm("eeg-data/nonzero-ts-666.csv",','))
list_names = vec(readdlm("eeg-data/names-relerr-10-80-666.csv",','))

suffix = "666"

ρ3 = zeros(n,0)
relerr = Float64[]
den2 = Float64[]
den3 = Float64[]

if do_list_bad_inference
	thresh = .9
	c = 0
	list = String[]
	rerr = Float64[]
#	for su in subjects
#		for st in states
#			name = "S"*su*"R"*st
	for name in list_names
			re = readdlm("eeg-data/"*name*"-"*type*"-re-"*suffix*".csv",',')[1]
			push!(rerr,re)
			if re > thresh
				@info "============================="
				@info name*": rel-err = $re"
				coeff = readdlm("eeg-data/"*name*"-"*type*"-coeff-"*suffix*".csv",',')
				@info "Max. coeff = $(maximum(abs.(coeff)))"
				figure(name)
				PyPlot.hist(vec(coeff),50)
				PyPlot.semilogy()
				global c += 1
				push!(list,name)
			end
#		end
	end
	@info "==========================="
	@info "Total: $c inferences with relative error higher than $thresh"
	cumrerr = [sum(rerr .> t) for t in LinRange(0,1,200)]
	PyPlot.plot(LinRange(0,1,200),cumrerr)
	xlabel("Relative error, r")
	ylabel("# of inferences with relative error > r")
end

if do_trajectory
	su = "015"
	st = "01"
	name = "S"*su*"R"*st
	i0 = 7

	s2signal = read_eeg("eeg-data/"*name*".edf")
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

	y0 = readdlm("eeg-data/"*name*"-"*type*"-coeff-"*suffix*".csv",',')[i0,1]
	A2 = readdlm("eeg-data/"*name*"-"*type*"-A2-"*suffix*".csv",',')
	A3 = readdlm("eeg-data/"*name*"-"*type*"-A3-"*suffix*".csv",',')

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

	re = round(readdlm("eeg-data/"*name*"-"*type*"-re-"*suffix*".csv",',')[1],digits=2)

	figure("Fit",(13,4))
	PyPlot.plot(Y[i0,:])
	PyPlot.plot(Ythis,"--")
	xlabel("time steps (restricted to the box)")
	ylabel("Derivative of sensor $i0")
	title("Fitting the derivative of sensor $i0, for subject "*name*" (relative error: $re)")
end


if do_hypergraph_density
	zer0 = 1e-6
#	for su in subjects
#		@info "S"*su*"Rxx"
#		for st in states
#			name = "S"*su*"R"*st
	for name in list_names
			A2 = readdlm("eeg-data/"*name*"-"*type*"-A2-"*suffix*".csv",',')
			A3 = readdlm("eeg-data/"*name*"-"*type*"-A3-"*suffix*".csv",',')
			push!(den2,sum(A2[:,3] .> zer0)/(n*binomial(n-1,1)))
			push!(den3,sum(A3[:,4] .> zer0)/(n*binomial(n-1,2)))
#		end
	end
	figure()
	ax = subplot(2,1,1)
	PyPlot.hist(den2,30)
	PyPlot.plot([0;sort(den2)],(0:length(den2))*100/length(den2))
	ylabel("occurences (out of $(length(den2)))")
	secax = ax.secondary_yaxis("right",functions=(x->x/100,x->x*100))
	secax.set_ylabel("cumulative")
	title("Histogram: hypergraph density (top: pairwise, bottom: triadic)")
	ax = subplot(2,1,2)
	PyPlot.hist(den3,30)
	PyPlot.plot([0;sort(den3)],(0:length(den3))*100/length(den3))
	xlabel("density")
	ylabel("occurences (out of $(length(list_names)))")
	secax = ax.secondary_yaxis("right",functions=(x->x/100,x->x*100))
	secax.set_ylabel("cumulative")
end


if do_relative_error
	remin = .1
	remax = .8
	keep_ts = String[]
#	for su in subjects
#		for st in states
#			name = "S"*su*"R"*st
	for name in list_names
			push!(relerr,readdlm("eeg-data/"*name*"-"*type*"-re-"*suffix*".csv",',')[1])
			if remin < relerr[end] < remax
				push!(keep_ts,name)
			end
#		end
	end

	fig,ax = PyPlot.subplots()
	PyPlot.hist(relerr,30)
	PyPlot.plot([0;sort(relerr)],56*(0:length(relerr))./length(relerr),"-",lw=2)
	PyPlot.plot([remin,remin],[0,100],"--k")
	PyPlot.plot([remax,remax],[0,100],"--k")
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
	
#	for su in subjects
#		for st in states
#			name = "S"*su*"R"*st
	for name in list_names
			a2 = readdlm("eeg-data/"*name*"-"*type*"-A2-"*suffix*".csv",',')
			t2 = norm(a2[:,3],1)
			a3 = readdlm("eeg-data/"*name*"-"*type*"-A3-"*suffix*".csv",',')
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
#		end
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

	# Sensor to zone pairing
	nz = 7
	s = readdlm("eeg-data/sensors-$nz.csv",',',String)
	z = readdlm("eeg-data/zones-$nz.csv",',',Int64)
	s2z = Dict{String,Int64}(s[i] => z[i] for i in 1:length(s))

	sensors = sort(collect(keys(s2z)))

	s2 = zeros(7)
	s3 = zeros(7)
	for k = 1:1000
		i2 = E2[k][1]
		z2 = s2z[sensors[i2]]
		s2[z2] += 1

		i3 = E3[k][1]
		z3 = s2z[sensors[i3]]
		s3[z3] += 1
	end

	@info "The 1000 most prominent 2-edges point towards area:"
	@info "1: $(s2[1])"
	@info "2: $(s2[2])"
	@info "3: $(s2[3])"
	@info "4: $(s2[4])"
	@info "5: $(s2[5])"
	@info "6: $(s2[6])"
	@info "7: $(s2[7])"
	
	@info "The 1000 most prominent 3-edges point towards area:"
	@info "1: $(s3[1])"
	@info "2: $(s3[2])"
	@info "3: $(s3[3])"
	@info "4: $(s3[4])"
	@info "5: $(s3[5])"
	@info "6: $(s3[6])"
	@info "7: $(s3[7])"

	
	#@info ""; @info ""; @info "The $n_edges most prominent k-edges are:"; @info ""; @info "| k |   2    |    3      |     4        |"; @info "-----------------------------------------"; for i in 1:n_edges; @info "|   | $(E2[i]) | $(E3[i]) | $(E4[i]) |"; end
end

 
if do_violin
	c = 0
	if !(data_exist)
#		for su in subjects
#			for st in states
#				name = "S"*su*"R"*st
		for name in list_names
				global c += 1
				@info "Working on "*name*": run "*suffix
				file = "eeg-data/"*name*".edf"
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
	
				z2 = zeros(n,T)
				z3 = zeros(n,T)
# #= NEW VERSION #####################

				idx_mon = get_idx_mon(n,dmax)
				coeff = readdlm("eeg-data/"*name*"-"*type*"-coeff-"*suffix*".csv",',')
				Ainf = Dict{Int64,Matrix{Float64}}(o => zeros(0,o+1) for o in 1:dmax+1)
				idx_coeff = Dict{Int64,Vector{Int64}}()
				nz_idx = Int64[]
				for i in keys(idx_mon)
					if length(idx_mon[i]) < 3
						aaa = setdiff((1:n)[abs.(coeff[:,i]) .> 1e-8],idx_mon[i])
						if !isempty(aaa)
							push!(nz_idx,i)
							idx_coeff[i] = aaa
						end
					end
				end

				for id in nz_idx
					ii = idx_coeff[id]
					jj = idx_mon[id]
					o = length(jj)+1
					a = coeff[ii,id]
					XX = prod(X[jj,:],dims=1)
					if o == 2
						z2[ii,:] += abs.(a*XX)
					elseif o == 3
						z3[ii,:] += abs.(a*XX)
					end
				end
# =#


 #= OLD VERSION #######################
				a2 = readdlm("eeg-data/S"*su*"R"*st*"-"*type*"-A2-"*suffix*".csv",',')
				t2 = norm(a2[:,3],1)
				a3 = readdlm("eeg-data/S"*su*"R"*st*"-"*type*"-A3-"*suffix*".csv",',')
				t3 = norm(a3[:,4],1)
	
				writedlm("eeg-data/T-S"*su*"R"*st*"-"*suffix*".csv",T,',')
				z2 = zeros(n,T)
				for l in 1:size(a2)[1]
					i,j = Int64.(a2[l,1:2])
					e = [i,j]
					a = a2[l,3]
					z2[i,:] += abs.(a*X[j,:])
				end
				z3 = zeros(n,T)
				for l in 1:size(a3)[1]
					i,j,k = Int64.(a3[l,1:3])
					e = [i,j,k]
					a = a3[l,4]
					z3[i,:] += abs.(a*X[j,:].*X[k,:])
				end
# =#

				z = z2+z3
				ρ3 = (z .> 0.).*(z3./z)
				writedlm("eeg-data/violin-data3-"*type*"-$n-"*name*"-"*suffix*"-rho3.csv",ρ3,',')
				ρ2 = (z .> 0.).*(z2./z)
				writedlm("eeg-data/violin-data2-"*type*"-$n-"*name*"-"*suffix*"-rho2.csv",ρ2,',')
#			end
		end
	end

	ρ3 = zeros(n,0)
	ρ2 = zeros(n,0)
	contr_per_subject = zeros(n,0)
	q0 = zeros(n,0)
	q1 = zeros(n,0)
	q2 = zeros(n,0)
	q3 = zeros(n,0)
	q4 = zeros(n,0)
	nams = String[]
#	for su in subjects
#		for st in states
#			name = "S"*su*"R"*st
	for name in list_names
			@info "Loading "*name
	
			#T = Int64(1000*floor(readdlm("eeg-data/T-"*name*"-"*suffix*".csv",',')[1]/1000))
			ρ = readdlm("eeg-data/violin-data3-"*type*"-$n-"*name*"-"*suffix*"-rho3.csv",',')
			global ρ3 = [ρ3 ρ]
			global ρ2 = [ρ2 readdlm("eeg-data/violin-data2-"*type*"-$n-"*name*"-"*suffix*"-rho2.csv",',')]
			global contr_per_subject = [contr_per_subject mean(ρ,dims=2)]
			global q0 = [q0 [quantile(ρ[i,:],.00) for i in 1:n]]
			global q1 = [q1 [quantile(ρ[i,:],.25) for i in 1:n]]
			global q2 = [q2 [quantile(ρ[i,:],.50) for i in 1:n]]
			global q3 = [q3 [quantile(ρ[i,:],.75) for i in 1:n]]
			global q4 = [q4 [quantile(ρ[i,:],1.0) for i in 1:n]]
			push!(nams,name)
#		end
	end
	
	figure("Violins-ter",(10,4))
	ζ = [vec(ρ3),]
#	ζ = [vec(ρ2),] # Does the same violins for the ratio of pairwise interactions
	qs = quantile(vec(q2),[.05,.35,.65,.95])
	#qs = quantile(vec(q2),[.2,.4,.6,.8])
	m1,i1 = findmin(abs.(q2 .- qs[1]))
	m2,i2 = findmin(abs.(q2 .- qs[2]))
	m3,i3 = findmin(abs.(q2 .- qs[3]))
	m4,i4 = findmin(abs.(q2 .- qs[4]))
	for i in [i1,i2,i3,i4]
		push!(ζ,ρ3[i[1],(i[2]-1)*Tmax .+ (1:Tmax)])
#		push!(ζ,ρ2[i[1],(i[2]-1)*Tmax .+ (1:Tmax)]) # Does the same violins for the ratio of pairwise interactions
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


if do_violin_4
	c = 0
	if !(data_exist)
#		for su in subjects
#			for st in states
#				name = "S"*su*"R"*st
		for name in list_names
				global c += 1
				@info "Working on "*name*": run "*suffix
				file = "eeg-data/"*name*".edf"
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
	
				z2 = zeros(n,T)
				z3 = zeros(n,T)
				z4 = zeros(n,T)

				idx_mon = get_idx_mon(n,dmax)
				coeff = readdlm("eeg-data/"*name*"-"*type*"-coeff-"*suffix*".csv",',')
				Ainf = Dict{Int64,Matrix{Float64}}(o => zeros(0,o+1) for o in 1:dmax+1)
				idx_coeff = Dict{Int64,Vector{Int64}}()
				nz_idx = Int64[]
				for i in keys(idx_mon)
					aaa = setdiff((1:n)[abs.(coeff[:,i]) .> 1e-8],idx_mon[i])
					if !isempty(aaa)
						push!(nz_idx,i)
						idx_coeff[i] = aaa
					end
				end

				for id in nz_idx
					ii = idx_coeff[id]
					jj = idx_mon[id]
					o = length(jj)+1
					a = coeff[ii,id]
					XX = prod(X[jj,:],dims=1)
					if o == 2
						z2[ii,:] += abs.(a*XX)
					elseif o == 3
						z3[ii,:] += abs.(a*XX)
					elseif o == 4
						z4[ii,:] += abs.(a*XX)
					end
				end
				
				z = z2+z3+z4
				ρ2 = (z .> 0.).*(z2./z)
				writedlm("eeg-data/trigon-data2-"*type*"-$n-"*name*"-"*suffix*"-rho2.csv",ρ2,',')
				ρ3 = (z .> 0.).*(z3./z)
				writedlm("eeg-data/trigon-data3-"*type*"-$n-"*name*"-"*suffix*"-rho3.csv",ρ3,',')
				ρ4 = (z .> 0.).*(z4./z)
				writedlm("eeg-data/trigon-data4-"*type*"-$n-"*name*"-"*suffix*"-rho4.csv",ρ4,',')
#			end
		end
	end


	ρ34 = zeros(n,0)
	ρ34vec = Float64[]
	ρ2 = zeros(n,0)
	contr_per_subject = zeros(n,0)
	q0 = zeros(n,0)
	q1 = zeros(n,0)
	q2 = zeros(n,0)
	q3 = zeros(n,0)
	q4 = zeros(n,0)
	nams = String[]
	nams_nz = String[]
#	for su in subjects
#		for st in states
#			name = "S"*su*"R"*st
	for name in list_names
			@info "Loading "*name
	
			#T = Int64(1000*floor(readdlm("eeg-data/T-"*name*"-"*suffix*".csv",',')[1]/1000))
			ρ = readdlm("eeg-data/trigon-data3-"*type*"-$n-"*name*"-"*suffix*"-rho3.csv",',') + readdlm("eeg-data/trigon-data4-"*type*"-$n-"*name*"-"*suffix*"-rho4.csv",',')
			ρi = readdlm("eeg-data/trigon-data2-"*type*"-$n-"*name*"-"*suffix*"-rho2.csv",',')
			ρ = (ρ+ρi .< .1).*1000. + (ρ+ρi .>= .1).*ρ
			#global ρ34 = hcat(ρ34,ρ)

			global ρ34 = hcat(ρ34,ρ)
			global ρ34vec = vcat(ρ34vec,ρ[(ρ+ρi) .< 2.])
			global ρ2 = hcat(ρ2,ρi)
			global contr_per_subject = [contr_per_subject mean(ρ,dims=2)]
			global q0 = [q0 [quantile(ρ[i,:],.00) for i in 1:n]]
			global q1 = [q1 [quantile(ρ[i,:],.25) for i in 1:n]]
			global q2 = [q2 [quantile(ρ[i,:],.50) for i in 1:n]]
			global q3 = [q3 [quantile(ρ[i,:],.75) for i in 1:n]]
			global q4 = [q4 [quantile(ρ[i,:],1.0) for i in 1:n]]
			push!(nams,name)
			if .1 < minimum([maximum(ρ[i,:]+ρi[i,:]) for i in 1:n]) < 1.1
				push!(nams_nz,name)
			end
#		end
	end
	
	figure("Violins-ter",(10,4))
	ζ = [ρ34vec,]
#	ζ = [vec(ρ2),] # Does the same violins for the ratio of pairwise interactions
#	qs = quantile(vec(q2),[.05,.35,.65,.95])
	qs = quantile(q2[q2 .< 2.],[.005,.3,.55,.8])
	m1,i1 = findmin(abs.(q2 .- qs[1]))
	m2,i2 = findmin(abs.(q2 .- qs[2]))
	m3,i3 = findmin(abs.(q2 .- qs[3]))
	m4,i4 = findmin(abs.(q2 .- qs[4]))
	for i in [i1,i2,i3,i4]
		push!(ζ,ρ34[i[1],(i[2]-1)*Tmax .+ (1:Tmax)])
#		push!(ζ,ρ2[i[1],(i[2]-1)*Tmax .+ (1:Tmax)]) # Does the same violins for the ratio of pairwise interactions
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
	#xticks([1,2,3,4,5],["All","5%","35%","65%","95%"])
	xticks([1,2,3,4,5],["All","0.5%","30%","55%","80%"])
	xlabel("Percentile")
	ylabel("Amount of the dynamics that is\nexplained by 3rd-order interactions")

	figure("Histogram: higher-order contribution")
	PyPlot.hist(ρ34vec,200)
	PyPlot.plot([.667,.667],[100,1e6],"--k")
	PyPlot.plot([.8547,.8547],[100,1e6],"--k")
	PyPlot.plot([.9109,.9109],[100,1e6],"--k")
	PyPlot.plot([.96,.96],[100,1e6],"--k")
	xlabel("Higher-order contribution")
	ylabel("Occurences")

	@info " 1% are <.667"
	@info " 5% are <.8547"
	@info "10% are <.9109"
	@info "20% are <.96"
	@info "36% are <1."
end

if do_triangle_full
	c = 0
	if !(data_exist)
#		for su in subjects
#			for st in states
#				name = "S"*su*"R"*st
		for name in list_names
				global c += 1
				@info "Working on "*name*": run "*suffix
				file = "eeg-data/"*name*".edf"
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
	
				z2 = zeros(n,T)
				z3 = zeros(n,T)
				z4 = zeros(n,T)

				idx_mon = get_idx_mon(n,dmax)
				coeff = readdlm("eeg-data/"*name*"-"*type*"-coeff-"*suffix*".csv",',')
				Ainf = Dict{Int64,Matrix{Float64}}(o => zeros(0,o+1) for o in 1:dmax+1)
				idx_coeff = Dict{Int64,Vector{Int64}}()
				nz_idx = Int64[]
				for i in keys(idx_mon)
					aaa = setdiff((1:n)[abs.(coeff[:,i]) .> 1e-8],idx_mon[i])
					if !isempty(aaa)
						push!(nz_idx,i)
						idx_coeff[i] = aaa
					end
				end

				for id in nz_idx
					ii = idx_coeff[id]
					jj = idx_mon[id]
					o = length(jj)+1
					a = coeff[ii,id]
					XX = prod(X[jj,:],dims=1)
					if o == 2
						z2[ii,:] += abs.(a*XX)
					elseif o == 3
						z3[ii,:] += abs.(a*XX)
					elseif o == 4
						z4[ii,:] += abs.(a*XX)
					end
				end
				
				z = z2+z3+z4
				ρ2 = (z .> 0.).*(z2./z)
				writedlm("eeg-data/trigon-data2-"*type*"-$n-"*name*"-"*suffix*"-rho2.csv",ρ2,',')
				ρ3 = (z .> 0.).*(z3./z)
				writedlm("eeg-data/trigon-data3-"*type*"-$n-"*name*"-"*suffix*"-rho3.csv",ρ3,',')
				ρ4 = (z .> 0.).*(z4./z)
				writedlm("eeg-data/trigon-data4-"*type*"-$n-"*name*"-"*suffix*"-rho4.csv",ρ4,',')
#			end
		end
	end

	ρ2 = zeros(n,0)
	ρ3 = zeros(n,0)
	ρ4 = zeros(n,0)
#	for su in subjects
#		for st in states
#			name = "S"*su*"R"*st
	for name in list_names
			@info "Loading "*name
			global ρ2 = [ρ2 readdlm("eeg-data/trigon-data2-"*type*"-$n-"*name*"-"*suffix*"-rho2.csv",',')]
			global ρ3 = [ρ3 readdlm("eeg-data/trigon-data3-"*type*"-$n-"*name*"-"*suffix*"-rho3.csv",',')]
			global ρ4 = [ρ4 readdlm("eeg-data/trigon-data4-"*type*"-$n-"*name*"-"*suffix*"-rho4.csv",',')]
#		end	
	end

	i_nz = (1:length(ρ2))[(vec(ρ2)+vec(ρ3)+vec(ρ4)) .> 1e-8]
	contour_trigon_data(vec(ρ3)[i_nz],vec(ρ4)[i_nz],100,"trigon","log","RdPu")
	plot_trigon_label("ρ3","ρ4","ρ2")
end



if do_triangle_part
	#nam = String[]; for su in list_all_subjects(109); for st in ["01","02"]; push!(nam,"S"*su*"R"*st); end; end
	#nam = nam[201:218]
	list_names = ["S001R01","S001R02","S002R01","S003R01","S005R02","S009R02","S021R01","S023R02","S024R01","S032R02","S034R01","S069R01","S069R02","S090R02"]
	c = 0
	if !(data_exist)
		for name in list_names
			global c += 1
			@info "Working on "*name*": run "*suffix
			file = "eeg-data/"*name*".edf"
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

			z2 = zeros(n,T)
			z3 = zeros(n,T)
			z4 = zeros(n,T)

			idx_mon = get_idx_mon(n,dmax)
			coeff = readdlm("eeg-data/"*name*"-"*type*"-coeff-"*suffix*".csv",',')
			Ainf = Dict{Int64,Matrix{Float64}}(o => zeros(0,o+1) for o in 1:dmax+1)
			idx_coeff = Dict{Int64,Vector{Int64}}()
			nz_idx = Int64[]
			for i in keys(idx_mon)
				aaa = setdiff((1:n)[abs.(coeff[:,i]) .> 1e-8],idx_mon[i])
				if !isempty(aaa)
					push!(nz_idx,i)
					idx_coeff[i] = aaa
				end
			end

			for id in nz_idx
				ii = idx_coeff[id]
				jj = idx_mon[id]
				o = length(jj)+1
				a = coeff[ii,id]
				XX = prod(X[jj,:],dims=1)
				if o == 2
					z2[ii,:] += abs.(a*XX)
				elseif o == 3
					z3[ii,:] += abs.(a*XX)
				elseif o == 4
					z4[ii,:] += abs.(a*XX)
				end
			end
			
			z = z2+z3+z4
			ρ2 = (z .> 0.).*(z2./z)
			writedlm("eeg-data/trigon-data2-"*type*"-$n-"*name*"-"*suffix*"-rho2.csv",ρ2,',')
			ρ3 = (z .> 0.).*(z3./z)
			writedlm("eeg-data/trigon-data3-"*type*"-$n-"*name*"-"*suffix*"-rho3.csv",ρ3,',')
			ρ4 = (z .> 0.).*(z4./z)
			writedlm("eeg-data/trigon-data4-"*type*"-$n-"*name*"-"*suffix*"-rho4.csv",ρ4,',')	
		end
	end

	for name in list_names
		@info "Loading "*name
		ρ2 = readdlm("eeg-data/trigon-data2-"*type*"-$n-"*name*"-"*suffix*"-rho2.csv",',')
		ρ3 = readdlm("eeg-data/trigon-data3-"*type*"-$n-"*name*"-"*suffix*"-rho3.csv",',')
		ρ4 = readdlm("eeg-data/trigon-data4-"*type*"-$n-"*name*"-"*suffix*"-rho4.csv",',')
		
		i_nz = (1:length(ρ2))[(vec(ρ2)+vec(ρ3)+vec(ρ4)) .> 1e-8]
		contour_trigon_data(vec(ρ3)[i_nz],vec(ρ4)[i_nz],100,"trigon: "*name,"log")
		plot_trigon_label("ρ3","ρ4","ρ2")
	end


end


