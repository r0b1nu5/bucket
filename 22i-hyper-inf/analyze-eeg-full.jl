using PyPlot, DelimitedFiles, Statistics, LinearAlgebra

include("eeg-tools.jl")

n = 64
type = "full" 
data_exist = false
Tmax = 1000
thr2 = 1.
thr3 = 1.
thr4 = 1.

n_edges = 10 # Number of most prominent edges to keep

# Data to be loaded
subjects = list_all_subjects(109)
#subjects = ["001","002"]
states = ["01","02"]
#states = ["03","07","11"]

suffix = "111"

#if data_exist
#	ρ3 = readdlm("eeg-data/violin-data-"*type*"-$n.csv",',')
#else
	ρ3 = zeros(n,0)
#end

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



 #= 
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

# =#


