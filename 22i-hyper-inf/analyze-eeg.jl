using PyPlot, DelimitedFiles

include("eeg-tools.jl")

n = 7
type = "avg" # "avg" or "sub"
data_exist = true

# Sensor to zone pairing
s = readdlm("eeg-data/sensors-$n.csv",',',String)
z = readdlm("eeg-data/zones-$n.csv",',',Int64)
s2z = Dict{String,Int64}(s[i] => z[i] for i in 1:length(s))

# Data to be loaded
subjects = list_all_subjects(109)
#subjects = ["001","002"]
states = ["01","02"]

if data_exist
	ρ3 = readdlm("eeg-data/violin-data-"*type*"-$n.csv",',')
else
	ρ3 = zeros(n,0)
end

AA2 = zeros(n,n)
AA3 = zeros(n,n,n)
for su in subjects
	for st in states
		@info "Working on S"*su*"R"*st
		if !(data_exist)
			if type == "avg"
				s2signal = read_eeg("eeg-data/S"*su*"R"*st*".edf")
				asig = average_over_zones(s2signal,s2z)
				dt = 1/160
				truncat = findmin(vec(maximum(abs.([asig zeros(n)]),dims=1)) .> 1e-6)[2] - 1
				X0 = asig[:,1:truncat]
				X = X0[:,1:end-1]
			elseif type == "sub"
				s2signal = read_eeg("eeg-data/S"*su*"R"*st*".edf")
				ssig = zeros(0,length(s2signal["Pz.."]))
				for s in slist
					ssig = [ssig;s2signal[s]']
				end
				dt = 1/160
				truncat = findmin(vec(maximum(abs.([ssig zeros(n)]),dims=1)) .> 1e-6)[2] - 1
				X0 = ssig[:,1:truncat]
				X = X0[:,1:end-1]
			end
			
			T = size(X)[2]
		end

		A2 = readdlm("eeg-data/S"*su*"R"*st*"-"*type*"-A2-bis.csv",',')
		a3 = readdlm("eeg-data/S"*su*"R"*st*"-"*type*"-A3-bis.csv",',')
		A3 = zeros(n,n,n)
		for k in 1:n
			A3[:,:,k] = a3[:,(k-1)*n .+ (1:n)]
		end

		if !(data_exist)
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
#				global ρ3 = [ρ3 ρ]
				writedlm("eeg-data/violin-data-"*type*"-$n-S"*su*"R"*st*".csv",ρ,',')
			end
		end

		global AA2 += (abs.(A2) .> 1e-6)
		global AA3 += (abs.(A3) .> 1e-6)
	end
end

for su in subjects
	for st in states
		global ρ3 = [ρ3 readdlm("eeg-data/violin-data-"*type*"-$n-S"*su*"R"*st*".csv",',')]
	end
end

if !(data_exist)
	writedlm("eeg-data/violin-data-"*type*"-$n.csv",ρ3,',')
end


# Fig 1
figure("Violins",(15,6))
#for i in 1:n
#	PyPlot.plot(i*ones(length(ρ3[i,:])),ρ3[i,:],"xk",alpha=.002)
#end
fig = plt.violinplot(ρ3',showextrema=false,showmedians=false)
for pc in fig["bodies"]
	pc.set_facecolor("#d41367")
	pc.set_edgecolor("black")
	pc.set_alpha(.3)
end
for i in 1:n
	PyPlot.plot([i,i],[minimum(ρ3[i,:]),maximum(ρ3[i,:])],"k")
	PyPlot.plot([i,i],quantile(ρ3[i,:],[.25,.75]),"k",linewidth=5)
	PyPlot.plot(i,median(ρ3[i,:]),"ok",markersize=10)
end

xlabel("Area")
ylabel("Amount of the dynamics that is explained by 3rd-order interactions")

# Fig 2
N = length(subjects)*length(states)
M2 = maximum(AA2)
M3 = maximum(AA3)

figure("Histograms - average - $n",(15,6))
subplot(2,1,1)
#PyPlot.hist(100*vec(AA2)./N,bins=100*((0:.5:M2+.5) .- .25)./N,color="#d41367")
PyPlot.hist(100*vec(AA2)./N,bins=100*((0:1:M2+.5) .- .25)./N,color="#d41367")
PyPlot.plot([-25/N,(100*M2+25)/N],[n^2,n^2],"--k")
dM = 10*(M2/N > .5) + 5*(.05 < M2/N < .5) + 1*(M2/N < .05)
xticks(0:dM:(100*M2/N + dM - 1e-4))
#xticks(100*(0:ceil(M2/10):M2)./N)
ylabel("# of 2-edges")
subplot(2,1,2)
#PyPlot.hist(100*vec(AA3)./N,bins=100*((0:.5:M3+.5) .- .25)./N,color="#d41367")
PyPlot.hist(100*vec(AA3)./N,bins=100*((0:1:M3+.5) .- .25)./N,color="#d41367")
PyPlot.plot([-25/N,(100*M3+25)/N],[n^3,n^3],"--k")
dM = 10*(M3/N > .5) + 5*(.05 < M3/N < .5) + 1*(M3/N < .05)
xticks(0:dM:(100*M2/N + dM - 1e-4))
#xticks(100*(0:ceil(M3/10):M3)./N)
xlabel("Appears in x% of the inferred hypergraphs")
ylabel("# of 3-edges")

# Fig 3
figure("Brains",(15,8))
nr = 3
nc = 4
c = 0
nrc = nr*nc

a3temp = copy(AA3)
for i in 1:nr
	for j in 1:nc
		global c += 1
		v,idx = findmax(a3temp)
		ids = [idx[k] for k in 1:3]
		a3temp[idx] = -1000.
		subplot(nr,nc,c)
		plot_brain_3dge(ids,v/N)
	end
end







