include("hyper_inf.jl")
include("eeg-tools.jl")

n = 7

# Sensor selection
slist = readdlm("eeg-data/slist-$n.csv",',',String)

# Data to be loaded
# #=
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
# =#
#subjects = ["001",]
states = ["01","02"]

AA2 = zeros(Int64,n,n)
AA3 = zeros(Int64,n,n,n)

for subject in subjects
	for state in states
		@info "----------------------------"
		@info "Running S"*subject*"R"*state

		# Loading data
		file = "eeg-data/S"*subject*"R"*state*".edf"
		s2signal = read_eeg(file)
		ssig = zeros(0,length(s2signal[keys(s2signal)[1]]))
		for s in slist
			ssig = [ssig;s2signal[s]']
		end

		# Finite differences
		dt = 1/160
		truncat = findmin(vec(maximum(abs.([ssig zeros(n)]),dims=1)) .> 1e-6)[2] - 1
		
		X0 = ssig[:,1:truncat]
		X = X0[:,1:end-1]
		Y = (X0[:,2:end]-X0[:,1:end-1])/dt
		
		for sen in slist
			X0 = [X0;s2signal[sen][1:end-truncate]']
		end
		X = X0[:,1:end-1]
		Y = (X0[:,2:end]-X0[:,1:end-1])/dt

		# Inference
		ooi = [2,3]
		dmax = 4
		xxx = hyper_inf(X,Y,ooi,dmax)

		# Retrieve adjacency tensors
		A2 = inferred_adj_2nd(xxx[1][2],n)[2]
		A3 = inferred_adj_3rd(xxx[1][3],n)[2]
		B2 = (abs.(A2) .> 1e-8)
		B3 = (abs.(A3) .> 1e-8)

		# Collect all boolean adjacency tensors
		global AA2 += B2
		global AA3 += B3

		# Plots (not too many)
		if length(subjects)*length(states) < 10
			figure("Histograms - subsampling-7 - S"*subject*"R"*state)
			subplot(2,1,1)
			PyPlot.hist(vec(abs.(A2)),20)
			subplot(2,1,2)
			PyPlot.hist(vec(abs.(A3)),20)
		end
		
		writedlm("eeg-data/S"*subject*"R"*state*"-sub-A2-bis.csv",A2,',')
		x = zeros(n,0)
		for i in 1:n
			x = [x A3[:,:,i]]
		end
		writedlm("eeg-data/S"*subject*"R"*state*"-sub-A3-bis.csv",x,',')
	end
end

N = length(subjects)*length(states)
M2 = maximum(AA2)
M3 = maximum(AA3)

figure("Histograms - subsampling-$n")
subplot(2,1,1)
PyPlot.hist(100*vec(AA2)./N,bins=100*((0:.5:M2+.5) .- .25)./N)
xticks(100*(0:ceil(M2/10):M2)./N)
ylabel("# of 2-edges")
subplot(2,1,2)
PyPlot.hist(100*vec(AA3)./N,bins=100*((0:.5:M3+.5) .- .25)./N)
xticks(100*(0:ceil(M3/10):M3)./N)
xlabel("Appears in x% of the inferred hypergraphs")
ylabel("# of 3-edges")




