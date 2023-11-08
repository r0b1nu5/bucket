using PyPlot, DelimitedFiles

include("eeg-tools.jl")

n = 7
type = "avg" # "avg" or "sub"

# Sensor to zone pairing
s = readdlm("eeg-data/sensors-$n.csv",',',String)
z = readdlm("eeg-data/zones-$n.csv",',',Int64)
s2z = Dict{String,Int64}(s[i] => z[i] for i in 1:length(s))

# Data to be loaded
subjects = String[]
for i in 1:9
	push!(subjects,"00"*string(i))
end
for i in 10:99
	push!(subjects,"0"*string(i))
end
for i in 100:109
	push!(subjects,string(i))
end
subjects = ["001","002"]
states = ["01","02"]

ρ3 = zeros(n,0)

AA2 = zeros(n,n)
AA3 = zeros(n,n,n)
for su in subjects
	for st in states
		if type == "avg"
			s2signal = read_eeg("eeg-data/S"*su*"R"*st*".edf")
			asig = average_over_zones(s2signal,s2z)
			dt = 1/160
			truncat = findmin(vec(maximum(abs.([asig zeros(n)]),dims=1)) .> 1e-6)[2] - 1
			X0 = asig[:,1:truncat]
			X = X0[:,1:end-1]
		elseif type == "sub"
			s2signal = read_eeg("eeg-data/S"*su*"R"*st*".edf")
			ssig = zeros(0,length(s2signal[keys(s2signal)[1]]))
			for s in slist
				ssig = [ssig;s2signal[s]']
			end
			dt = 1/160
			truncat = findmin(vec(maximum(abs.([ssig zeros(n)]),dims=1)) .> 1e-6)[2] - 1
			X0 = ssig[:,1:truncat]
			X = X0[:,1:end-1]
		end
		
		T = size(X)[2]

		A2 = readdlm("eeg-data/S"*su*"R"*st*"-"*type*"-A2-bis.csv",',')
		a3 = readdlm("eeg-data/S"*su*"R"*st*"-"*type*"-A3-bis.csv",',')
		A3 = zeros(n,n,n)
		for k in 1:n
			A3[:,:,k] = a3[:,(k-1)*n .+ (1:n)]
		end

		for t in 1:T
			ρ = Float64[]
			for i in 1:n
				z2 = 0.
				z3 = 0.
				for j in 1:n
					for k in 1:n
						z3 += abs(A3[i,j,k]*X[j,t]*X[k,t])
					end
					z2 += abs(A2[i,j]*X[j,t])
				end
				push!(ρ,z3/(z2+z3))
			end
			global ρ3 = [ρ3 ρ]
		end

		global AA2 += (abs.(A2) .> 1e-6)
		global AA3 += (abs.(A3) .> 1e-6)
	end
end

# Fig 1
figure("Violins")
#for i in 1:n
#	PyPlot.plot(i*ones(length(ρ3[i,:])),ρ3[i,:],"xk",alpha=.002)
#end
plt.violinplot(ρ3',showextrema=false,showmedians=true)
xlabel("Area")
ylabel("Amount of the dynamics that is explained by 3rd-order interactions")

# Fig 2
N = length(subjects)*length(states)
M2 = maximum(AA2)
M3 = maximum(AA3)

figure("Histograms - average - $n")
subplot(2,1,1)
PyPlot.hist(100*vec(AA2)./N,bins=100*((0:.5:M2+.5) .- .25)./N,density=true)
xticks(100*(0:ceil(M2/10):M2)./N)
ylabel("# of 2-edges")
subplot(2,1,2)
PyPlot.hist(100*vec(AA3)./N,bins=100*((0:.5:M3+.5) .- .25)./N)
xticks(100*(0:ceil(M3/10):M3)./N)
xlabel("Appears in x% of the inferred hypergraphs")
ylabel("# of 3-edges")








