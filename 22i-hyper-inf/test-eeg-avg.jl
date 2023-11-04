include("hyper_inf.jl")
include("eeg-tools.jl")

nz = 7

# Sensor to zone pairing
s = readdlm("eeg-data/sensors-$nz.csv",',',String)
z = readdlm("eeg-data/zones-$nz.csv",',',Int64)
s2z = Dict{String,Int64}(s[i] => z[i] for i in 1:length(s))

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
subjects = ["001",]
states = ["01","02"]

AA2 = zeros(Int64,nz,nz)
AA3 = zeros(Int64,nz,nz,nz)

for subject in subjects
	for state in states
		@info "Running S"*subject*"R"*state

		# Loading data
		file = "eeg-data/S"*subject*"R"*state*".edf"
		s2signal = read_eeg(file)
		
		# Average signal over zones
		asig = average_over_zones(s2signal,s2z)
		
		# Finite differences
		dt = 1/160
		truncate = 128
		
		X0 = asig[:,1:end-truncate]
		X = X0[:,1:end-1]
		Y = (X0[:,2:end]-X0[:,1:end-1])/dt

		# Inference
		ooi = [2,3]
		dmax = 4
		xxx = hyper_inf(X,Y,ooi,dmax)
		
		# Retrieve adjacency tensors	
		A2 = inferred_adj_2nd(xxx[1][2],nz)[2]
		A3 = inferred_adj_3rd(xxx[1][3],nz)[2]
		B2 = (abs.(A2) .> 1e-8)
		B3 = (abs.(A3) .> 1e-8)

		# Collect all boolean adjacency tensors
		global AA2 += B2
		global AA3 += B3

		# Plots (not too many)
		if length(subjects)*length(states) < 10
			figure("Histograms - average-$nz - S"*subject*"R"*state)
			subplot(2,1,1)
			PyPlot.hist(vec(abs.(A2)),20)
			subplot(2,1,2)
			PyPlot.hist(vec(abs.(A3)),20)
		end
		
		writedlm("eeg-data/S"*subject*"R"*state*"-avg-A2.csv",A2,',')
		x = zeros(nz,0)
		for i in 1:nz
			x = [x A3[:,:,i]]
		end
		writedlm("eeg-data/S"*subject*"R"*state*"-avg-A3.csv",x,',')
	end
end

figure("Histograms - average-$nz")
subplot(2,1,1)
PyPlot.hist(vec(AA2),maximum(AA2)+1)
subplot(2,1,2)
PyPlot.hist(vec(AA3),maximum(AA3)+1)








