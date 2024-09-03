# Same as "analyze-eeg-4th.jl" but assuming that adjacency tensors are stored as lists

using PyPlot, DelimitedFiles, Statistics, LinearAlgebra

include("eeg-tools.jl")

n = 7
type = "avg" # "avg" or "sub"
data_exist = true
Tmax = 1000

n_edges = 10 # Number of most prominent edges to keep

# Data to be loaded
subjects = list_all_subjects(109)
#subjects = ["001","002"]
states = ["01","02"]
#states = ["03","07","11"]

suffix = "44x"
λ = .1

# Sensor to zone pairing
s = readdlm("eeg-data/sensors-$n.csv",',',String)
z = readdlm("eeg-data/zones-$n.csv",',',Int64)
s2z = Dict{String,Int64}(s[i] => z[i] for i in 1:length(s))

edge_score2 = Dict{Vector{Int64},Float64}()
edge_score3 = Dict{Vector{Int64},Float64}()
edge_score4 = Dict{Vector{Int64},Float64}()

AA2 = zeros(n,n)
AAA2 = zeros(n,n,length(subjects)*length(states))
AA3 = zeros(n,n,n)
AAA3 = zeros(n,n,n,length(subjects)*length(states))
AA4 = zeros(n,n,n,n)
AAA4 = zeros(n,n,n,n,length(subjects)*length(states))
vA2 = Float64[]
vA3 = Float64[]
vA4 = Float64[]
contr_per_subject = zeros(n,0)

c = 0
constant_interaction_ratio = Float64[]

nthresh = 21
ths = [Inf;LinRange(.5,0,nthresh-1)]
for su in subjects
	for st in states
		global c += 1
		@info "Working on S"*su*"R"*st*": run "*suffix
		if !(data_exist)
			s2signal = read_eeg("eeg-data/S"*su*"R"*st*".edf")
			asig = average_over_zones(s2signal,s2z)
			dt = 1/160
			truncat = findmin(vec(maximum(abs.([asig zeros(n)]),dims=1)) .> 1e-6)[2] - 1
#			truncat = 3001
			X0 = denoise_fourier(asig[:,1:truncat],300)
			X0 = X0[:,1:end-1]
			X = X0./mean(abs.(X0))
		
			X,ids = restrict_box_size(X,1000)
			T = size(X)[2]
		end

		a2 = readdlm("eeg-data/S"*su*"R"*st*"-"*type*"-A2-"*suffix*".csv",',')
		t2 = norm(a2[:,3],1)
		a3 = readdlm("eeg-data/S"*su*"R"*st*"-"*type*"-A3-"*suffix*".csv",',')
		t3 = norm(a3[:,4],1)
		a4 = readdlm("eeg-data/S"*su*"R"*st*"-"*type*"-A4-"*suffix*".csv",',')
		t4 = norm(a4[:,5],1)

		if !(data_exist)
			z1 = abs.(readdlm("eeg-data/S"*su*"R"*st*"-avg-coeff-"*suffix*".csv",',')[:,1])
			writedlm("eeg-data/T-S"*su*"R"*st*"-"*suffix*".csv",T,',')
			for t in 1:T
				z2 = zeros(n)
				z3 = zeros(n)
				z4 = zeros(n)
			for p in 2:nthresh
				for l in (1:size(a2)[1])[ths[p] .< abs.(a2[:,3]) .< ths[p-1]]
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
				for l in (1:size(a3)[1])[ths[p] .< abs.(a3[:,4]) .< ths[p-1]]
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
				for l in (1:size(a4)[1])[ths[p] .< abs.(a4[:,5]) .< ths[p-1]]
					i,j,k,kk = Int64.(a4[l,1:4])
					e = [i,j,k,kk]
					a = a4[l,5]
					z4[i] += abs(a*X[j,t]*X[k,t]*X[kk,t])
					if haskey(edge_score4,e)
						edge_score4[e] += abs(a)/t4
					else
						edge_score4[e] = abs(a)/t4
					end
				end
				z = z2+z3+z4
				append!(constant_interaction_ratio,z1./z)

				ρ2 = [z[i] == 0. ? 0 : z2[i]/z[i] for i in 1:n]
				ρ3 = [z[i] == 0. ? 0 : z3[i]/z[i] for i in 1:n]
				ρ4 = [z[i] == 0. ? 0 : z4[i]/z[i] for i in 1:n]
				writedlm("eeg-data/trigon-data2-"*type*"-$n-S"*su*"R"*st*"-"*suffix*"-$t-th$p.csv",ρ2,',')
				writedlm("eeg-data/trigon-data3-"*type*"-$n-S"*su*"R"*st*"-"*suffix*"-$t-th$p.csv",ρ3,',')
				writedlm("eeg-data/trigon-data4-"*type*"-$n-S"*su*"R"*st*"-"*suffix*"-$t-th$p.csv",ρ4,',')
			end
			end
			ρ2 = zeros(n,0)
			ρ3 = zeros(n,0)
			ρ4 = zeros(n,0)
			for p in 2:nthresh
			for t in 1:T
				ρ2 = hcat(ρ2,vec(readdlm("eeg-data/trigon-data2-"*type*"-$n-S"*su*"R"*st*"-"*suffix*"-$t-th$p.csv",',')))
				rm("eeg-data/trigon-data2-"*type*"-$n-S"*su*"R"*st*"-"*suffix*"-$t-th$p.csv")
				ρ3 = hcat(ρ3,vec(readdlm("eeg-data/trigon-data3-"*type*"-$n-S"*su*"R"*st*"-"*suffix*"-$t-th$p.csv",',')))
				rm("eeg-data/trigon-data3-"*type*"-$n-S"*su*"R"*st*"-"*suffix*"-$t-th$p.csv")
				ρ4 = hcat(ρ4,vec(readdlm("eeg-data/trigon-data4-"*type*"-$n-S"*su*"R"*st*"-"*suffix*"-$t-th$p.csv",',')))
				rm("eeg-data/trigon-data4-"*type*"-$n-S"*su*"R"*st*"-"*suffix*"-$t-th$p.csv")
				if t%1000 == 0
					writedlm("eeg-data/trigon-data2-"*type*"-$n-S"*su*"R"*st*"-"*suffix*"-T$t-th$p.csv",ρ2,',')
					writedlm("eeg-data/trigon-data3-"*type*"-$n-S"*su*"R"*st*"-"*suffix*"-T$t-th$p.csv",ρ3,',')
					writedlm("eeg-data/trigon-data4-"*type*"-$n-S"*su*"R"*st*"-"*suffix*"-T$t-th$p.csv",ρ4,',')
					ρ2 = zeros(n,0)
					ρ3 = zeros(n,0)
					ρ4 = zeros(n,0)
				end
			end
			end
		end

	end
end

 #=
PyPlot.hist(constant_interaction_ratio,30)
xlabel("z1/(z2+z3+z4)")
ylabel("# occurences (out of 1,526,000)")
title("Histogram: relative contribution of the constant term")

re = readdlm("eeg-data/relative-error-"*suffix*".csv",',')
# =#

q20 = Float64[]
q21 = Float64[]
q22 = Float64[]
q23 = Float64[]
q24 = Float64[]
q30 = Float64[]
q31 = Float64[]
q32 = Float64[]
q33 = Float64[]
q34 = Float64[]
q40 = Float64[]
q41 = Float64[]
q42 = Float64[]
q43 = Float64[]
q44 = Float64[]
discard = Int64[]

for p in 3:nthresh
global ρ2 = zeros(n,0)
global ρ3 = zeros(n,0)
global ρ4 = zeros(n,0)
for su in subjects
	for st in states
		@info "Load: S"*su*"R"*st
		global ρ2 = hcat(ρ2,readdlm("eeg-data/trigon-data2-"*type*"-$n-S"*su*"R"*st*"-"*suffix*"-T1000-th$p.csv",','))
		global ρ3 = hcat(ρ3,readdlm("eeg-data/trigon-data3-"*type*"-$n-S"*su*"R"*st*"-"*suffix*"-T1000-th$p.csv",','))
		global ρ4 = hcat(ρ4,readdlm("eeg-data/trigon-data4-"*type*"-$n-S"*su*"R"*st*"-"*suffix*"-T1000-th$p.csv",','))
	end
end
v2 = vec(ρ2)
v3 = vec(ρ3)
v4 = vec(ρ4)
t = (v2+v3+v4 .> 1e-4)
r2 = v2[t]
r3 = v3[t]
r4 = v4[t]

if length(r2) > 100
push!(q20,quantile(r2,.10))
push!(q21,quantile(r2,.25))
push!(q22,quantile(r2,.50))
push!(q23,quantile(r2,.75))
push!(q24,quantile(r2,.90))
push!(q30,quantile(r3,.10))
push!(q31,quantile(r3,.25))
push!(q32,quantile(r3,.50))
push!(q33,quantile(r3,.75))
push!(q34,quantile(r3,.90))
push!(q40,quantile(r4,.10))
push!(q41,quantile(r4,.25))
push!(q42,quantile(r4,.50))
push!(q43,quantile(r4,.75))
push!(q44,quantile(r4,.90))
else
push!(discard,p)
end
end

contour_trigon_data(vec(ρ3),vec(ρ4),100,"trigon")
plot_trigon_label("ρ3","ρ4","ρ2")

card = setdiff(3:nthresh,discard)

figure("Contribution vs threshold",(5,8))
subplot(3,1,1)
PyPlot.fill([ths[card];ths[card[end:-1:1]]],[q20;q24[end:-1:1]],color="C0",alpha=.3)
PyPlot.fill([ths[card];ths[card[end:-1:1]]],[q21;q23[end:-1:1]],color="C0",alpha=.5)
PyPlot.plot(ths[card],q22,color="C0")
PyPlot.plot([λ,λ],[0,1],"--k")
ylabel("ρ2")
subplot(3,1,2)
PyPlot.fill([ths[card];ths[card[end:-1:1]]],[q30;q34[end:-1:1]],color="C1",alpha=.3)
PyPlot.fill([ths[card];ths[card[end:-1:1]]],[q31;q33[end:-1:1]],color="C1",alpha=.5)
PyPlot.plot(ths[card],q32,color="C1")
PyPlot.plot([λ,λ],[0,1],"--k")
ylabel("ρ3")
subplot(3,1,3)
PyPlot.fill([ths[card];ths[card[end:-1:1]]],[q40;q44[end:-1:1]],color="C2",alpha=.3)
PyPlot.fill([ths[card];ths[card[end:-1:1]]],[q41;q43[end:-1:1]],color="C2",alpha=.5)
PyPlot.plot(ths[card],q42,color="C2")
PyPlot.plot([λ,λ],[0,1],"--k")
xlabel("Threshold")
ylabel("ρ4")

 #=
ρ34 = zeros(n,0)
ρ34vec = Float64[]
contr_per_subject = zeros(n,0)
q0 = zeros(n,0)
q1 = zeros(n,0)
q2 = zeros(n,0)
q3 = zeros(n,0)
q4 = zeros(n,0)
nams = String[]
nams_nz = String[]
for su in subjects
	for st in states
		name = "S"*su*"R"*st
		@info "Loading "*name

		ρ = readdlm("eeg-data/trigon-data3-"*type*"-$n-"*name*"-"*suffix*"-T1000.csv",',') + readdlm("eeg-data/trigon-data4-"*type*"-$n-"*name*"-"*suffix*"-T1000.csv",',')
		ρi = readdlm("eeg-data/trigon-data2-"*type*"-$n-"*name*"-"*suffix*"-T1000.csv",',')
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
	end
end

figure("Violins-ter",(10,4))
ζ = [ρ34vec,]
#	ζ = [vec(ρ2),] # Does the same violins for the ratio of pairwise interactions
qs = quantile(vec(q2),[.05,.35,.65,.95])
m1,i1 = findmin(abs.(q2 .- qs[1]))
m2,i2 = findmin(abs.(q2 .- qs[2]))
m3,i3 = findmin(abs.(q2 .- qs[3]))
m4,i4 = findmin(abs.(q2 .- qs[4]))
for i in [i1,i2,i3,i4]
	push!(ζ,ρ34[i[1],(i[2]-1)*Tmax .+ (1:Tmax)])
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
ylabel("ρ(s,t)")




E2 = Vector{Int64}[]
S2 = Float64[]
for i in 1:min(n_edges,length(edge_score2))
	sc,id = findmax(edge_score2)
	push!(E2,id)
	push!(S2,sc)
	pop!(edge_score2,id)
end

E3 = Vector{Int64}[]
S3 = Float64[]
for i in 1:min(n_edges,length(edge_score3))
	sc,id = findmax(edge_score3)
	push!(E3,id)
	push!(S3,sc)
	pop!(edge_score3,id)
end

E4 = Vector{Int64}[]
S4 = Float64[]
for i in 1:min(n_edges,length(edge_score4))
	sc,id = findmax(edge_score4)
	push!(E4,id)
	push!(S4,sc)
	pop!(edge_score4,id)
end

@info ""
@info ""
@info "The $n_edges most prominent k-edges are:"
@info ""
@info "| k |   2    |    3      |     4        |"
@info "-----------------------------------------"
for i in 1:n_edges
	@info "|   | $(E2[i]) | $(E3[i]) | $(E4[i]) |"
end

#@info ""; @info ""; @info "The $n_edges most prominent k-edges are:"; @info ""; @info "| k |   2    |    3      |     4        |"; @info "-----------------------------------------"; for i in 1:n_edges; @info "|   | $(E2[i]) | $(E3[i]) | $(E4[i]) |"; end



# =#





