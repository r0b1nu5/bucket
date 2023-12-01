include("hyper_inf.jl")
include("eeg-tools.jl")

n = 8

# Sensor selection
llist = readdlm("eeg-data/llist-$n.csv",',',String)

# Data to be loaded
#subjects = list_all_subjects(109)
subjects = ["001",]
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
		lsig = zeros(0,length(s2signal["Pz.."]))
		for s in llist
			lsig = [lsig;s2signal[s]']
		end

		# Finite differences
		dt = 1/160
		truncat = findmin(vec(maximum(abs.([lsig zeros(n)]),dims=1)) .> 1e-6)[2] - 1
		
		X0 = smooth_time_series(lsig[:,1:truncat])
		X = X0[:,1:end-1]
		Y = (X0[:,2:end]-X0[:,1:end-1])/dt
		
		# Inference
		ooi = [2,3]
		dmax = 4
		xxx = hyper_inf(X,Y,ooi,dmax,.01)

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
			figure("Histograms - local sampling - 8 - S"*subject*"R"*state)
			subplot(2,1,1)
			PyPlot.hist(vec(abs.(A2)),20)
			subplot(2,1,2)
			PyPlot.hist(vec(abs.(A3)),20)
		end
		
		writedlm("eeg-data/S"*subject*"R"*state*"-loc-A2-bis.csv",A2,',')
		x = zeros(n,0)
		for i in 1:n
			x = [x A3[:,:,i]]
		end
		writedlm("eeg-data/S"*subject*"R"*state*"-loc-A3-bis.csv",x,',')
	end
end

N = length(subjects)*length(states)
M2 = maximum(AA2)
M3 = maximum(AA3)

figure("Histograms - local sampling - $n")
subplot(2,1,1)
PyPlot.hist(100*vec(AA2)./N,bins=100*((0:.5:M2+.5) .- .25)./N)
xticks(100*(0:ceil(M2/10):M2)./N)
ylabel("# of 2-edges")
subplot(2,1,2)
PyPlot.hist(100*vec(AA3)./N,bins=100*((0:.5:M3+.5) .- .25)./N)
xticks(100*(0:ceil(M3/10):M3)./N)
xlabel("Appears in x% of the inferred hypergraphs")
ylabel("# of 3-edges")




