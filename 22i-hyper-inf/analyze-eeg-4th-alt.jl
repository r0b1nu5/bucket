# Same as "analyze-eeg-4th.jl" but assuming that adjacency tensors are stored as lists

using PyPlot, DelimitedFiles, Statistics, LinearAlgebra

include("eeg-tools.jl")

n = 7
type = "avg" # "avg" or "sub"
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

suffix = "44x"

# Sensor to zone pairing
s = readdlm("eeg-data/sensors-$n.csv",',',String)
z = readdlm("eeg-data/zones-$n.csv",',',Int64)
s2z = Dict{String,Int64}(s[i] => z[i] for i in 1:length(s))

#if data_exist
#	ρ3 = readdlm("eeg-data/violin-data-"*type*"-$n.csv",',')
#else
	ρ3 = zeros(n,0)
#end

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
				z4 = zeros(n)
				for l in 1:size(a4)[1]
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

				ρ3 = [z[i] == 0. ? 0 : z3[i]/z[i] for i in 1:n]
				ρ4 = [z[i] == 0. ? 0 : z4[i]/z[i] for i in 1:n]
				writedlm("eeg-data/trigon-data3-"*type*"-$n-S"*su*"R"*st*"-"*suffix*"-$t.csv",ρ3,',')
				writedlm("eeg-data/trigon-data4-"*type*"-$n-S"*su*"R"*st*"-"*suffix*"-$t.csv",ρ4,',')
			end
			ρ3 = zeros(n,0)
			ρ4 = zeros(n,0)
			for t in 1:T
				ρ3 = [ρ3 vec(readdlm("eeg-data/trigon-data3-"*type*"-$n-S"*su*"R"*st*"-"*suffix*"-$t.csv",','))]
				rm("eeg-data/trigon-data3-"*type*"-$n-S"*su*"R"*st*"-"*suffix*"-$t.csv")
				ρ4 = [ρ4 vec(readdlm("eeg-data/trigon-data4-"*type*"-$n-S"*su*"R"*st*"-"*suffix*"-$t.csv",','))]
				rm("eeg-data/trigon-data4-"*type*"-$n-S"*su*"R"*st*"-"*suffix*"-$t.csv")
				if t%1000 == 0
					writedlm("eeg-data/trigon-data3-"*type*"-$n-S"*su*"R"*st*"-"*suffix*"-T$t.csv",ρ3,',')
					writedlm("eeg-data/trigon-data4-"*type*"-$n-S"*su*"R"*st*"-"*suffix*"-T$t.csv",ρ4,',')
					ρ3 = zeros(n,0)
					ρ4 = zeros(n,0)
				end
			end

		end

	end
end

PyPlot.hist(constant_interaction_ratio,30)
xlabel("z1/(z2+z3+z4)")
ylabel("# occurences (out of 1,526,000)")
title("Histogram: relative contribution of the constant term")

re = readdlm("eeg-data/relative-error-"*suffix*".csv",',')

ρ3 = zeros(n,0)
ρ4 = zeros(n,0)
for su in subjects
	for st in states
		global ρ3 = [ρ3 readdlm("eeg-data/trigon-data3-"*type*"-$n-S"*su*"R"*st*"-"*suffix*"-T1000.csv",',')]
		global ρ4 = [ρ4 readdlm("eeg-data/trigon-data4-"*type*"-$n-S"*su*"R"*st*"-"*suffix*"-T1000.csv",',')]
	end
end
contour_trigon_data(vec(ρ3),vec(ρ4),100,"trigon")
plot_trigon_label("ρ3","ρ4","ρ2")



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
