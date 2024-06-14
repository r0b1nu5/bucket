using PyPlot, DelimitedFiles, Statistics

include("eeg-tools.jl")

n = 7
type = "avg" # "avg" or "sub"
data_exist = true
Tmax = 1000
thr2 = 1.
thr3 = 1.
thr4 = 1.

# Data to be loaded
#subjects = list_all_subjects(109)
subjects = ["001","002"]
states = ["01","02"]
#states = ["03","07","11"]

suffix = "x11"

# Sensor to zone pairing
s = readdlm("eeg-data/sensors-$n.csv",',',String)
z = readdlm("eeg-data/zones-$n.csv",',',Int64)
s2z = Dict{String,Int64}(s[i] => z[i] for i in 1:length(s))

#if data_exist
#	ρ3 = readdlm("eeg-data/violin-data-"*type*"-$n.csv",',')
#else
	ρ3 = zeros(n,0)
#end

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

		A2 = readdlm("eeg-data/S"*su*"R"*st*"-"*type*"-A2-"*suffix*".csv",',')
		a3 = readdlm("eeg-data/S"*su*"R"*st*"-"*type*"-A3-"*suffix*".csv",',')
		A3 = zeros(n,n,n)
		for k in 1:n
			A3[:,:,k] = a3[:,(k-1)*n .+ (1:n)]
		end
		a4 = readdlm("eeg-data/S"*su*"R"*st*"-"*type*"-A4-"*suffix*".csv",',')
		A4 = zeros(n,n,n,n)
		for i in 1:n
			for j in 1:n
				A4[:,:,i,j] = a4[:,(i-1)*n^2 + (j-1)*n .+ (1:n)]
			end
		end
		global vA2 = [vA2;vectorize_adj(A2)]
		global vA3 = [vA3;vectorize_adj(A3)]
		global vA4 = [vA4;vectorize_adj(A4)]

		if !(data_exist)
			writedlm("eeg-data/T-S"*su*"R"*st*"-"*suffix*".csv",T,',')
			for t in 1:T
				ρ3 = Float64[]
				ρ4 = Float64[]
				for i in 1:n
					z2 = 0.
					z3 = 0.
					z4 = 0.
					for j in 1:n
						for k in 1:n
							for l in 1:n
								z4 += abs(A4[i,j,k,l]*X[j,t]*X[k,t]*X[k,l])
							end
							z3 += abs(A3[i,j,k]*X[j,t]*X[k,t])
						end
						z2 += abs(A2[i,j]*X[j,t])
					end
					if z2+z3+z4 == 0.
						push!(ρ3,0.)
						push!(ρ4,0.)
					else
						push!(ρ3,z3/(z2+z3+z4))
						push!(ρ4,z4/(z2+z3+z4))
					end
				end
#				global ρ3 = [ρ3 ρ]
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

		global AA2 += (abs.(A2) .> thr2)
		global AAA2[:,:,c] = abs.(A2)
		global AA3 += (abs.(A3) .> thr3)
		global AAA3[:,:,:,c] = abs.(A3)
		global AA4 += (abs.(A4) .> thr4)
		global AAA4[:,:,:,:,c] = abs.(A4)
	end
end

MA2 = zeros(n,n)
MA3 = zeros(n,n,n)
MA4 = zeros(n,n,n,n)
for i in 1:n
	for j in 1:n
		MA2[i,j] = median(AAA2[i,j,:])
		for k in 1:n
			MA3[i,j,k] = median(AAA3[i,j,k,:])
			for l in 1:n
				MA4[i,j,k,l] = median(AAA4[i,j,k,l,:])
			end
		end
	end
end

re = readdlm("eeg-data/relative-error-"*suffix*".csv",',')








