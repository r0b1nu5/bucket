include("hyper_inf.jl")
include("eeg-tools.jl")

nz = 7
λ = .1
ρ = .1

suffix = "x11"

# Data to be loaded
#subjects = list_all_subjects(109)
subjects = ["001","002"]
#subjects = ["021",]
states = ["01","02"]
#states = ["03","07","11"]
#states = ["01","02","03","07","11"]

# Sensor to zone pairing
s = readdlm("eeg-data/sensors-$nz.csv",',',String)
z = readdlm("eeg-data/zones-$nz.csv",',',Int64)
s2z = Dict{String,Int64}(s[i] => z[i] for i in 1:length(s))

AA2 = zeros(Int64,nz,nz)
AA3 = zeros(Int64,nz,nz,nz)
AA4 = zeros(Int64,nz,nz,nz,nz)

relerr = zeros(length(states),0)

for subject in subjects
	re = Float64[]
	for state in states
		@info "Running S"*subject*"R"*state

		# Loading data
		file = "eeg-data/S"*subject*"R"*state*".edf"
		s2signal = read_eeg(file)
		
		# Average signal over zones
		asig = average_over_zones(s2signal,s2z)
		
		# Finite differences
		dt = 1/160
		truncat = findmin(vec(maximum(abs.([asig zeros(nz)]),dims=1)) .> 1e-6)[2] - 1
#		truncat = 3001

#		X0 = asig[:,1:truncat]
#		X0 = denoise_fourier(asig[:,1:truncat],100)
		X0 = denoise_fourier(asig[:,1:truncat],200)
		Y0 = (X0[:,2:end]-X0[:,1:end-1])./dt
		X0 = X0[:,1:end-1]
		X = X0./mean(abs.(X0))
		Y = Y0./mean(abs.(Y0))

		X,ids = restrict_box_size(X,1000)
		Y = Y[:,ids]

		# Inference
		ooi = [2,3,4]
		dmax = 5
		xxx = hyper_inf(X,Y,ooi,dmax,λ,ρ)
		push!(re,xxx[4])
		
		# Retrieve adjacency tensors	
		A2 = inferred_adj_2nd(xxx[1][2],nz)[2]
		A3 = inferred_adj_3rd(xxx[1][3],nz)[2]
		A4 = inferred_adj_4th(xxx[1][4],nz)[2] ############################
		B2 = (abs.(A2) .> 1e-8)
		B3 = (abs.(A3) .> 1e-8)
		B4 = (abs.(A4) .> 1e-8) ############################

		# Collect all boolean adjacency tensors
		global AA2 += B2
		global AA3 += B3
		global AA4 += B4 ############################

		# Plots (not too many)
		if length(subjects)*length(states) < 1
			figure("Histograms - average-$nz - S"*subject*"R"*state)
			subplot(2,1,1)
			PyPlot.hist(vec(abs.(A2)),20)
			subplot(2,1,2)
			PyPlot.hist(vec(abs.(A3)),20)
		end
		
		writedlm("eeg-data/S"*subject*"R"*state*"-avg-A2-"*suffix*".csv",A2,',')
		x = zeros(nz,0)
		for i in 1:nz
			x = [x A3[:,:,i]]
		end
		writedlm("eeg-data/S"*subject*"R"*state*"-avg-A3-"*suffix*".csv",x,',')
		y = zeros(nz,0) ############################
		for i in 1:nz ############################
			for j in 1:nz ############################
				y = [y A4[:,:,i,j]] ############################
			end ############################
		end ############################
		writedlm("eeg-data/S"*subject*"R"*state*"-avg-A4-"*suffix*".csv",y,',') ############################

		writedlm("eeg-data/S"*subject*"R"*state*"-avg-coeff-"*suffix*".csv",xxx[2],',')	end
	global relerr = [relerr re]
end

writedlm("eeg-data/relative-error-"*suffix*".csv",relerr,',')

N = length(subjects)*length(states)
M2 = maximum(AA2)
M3 = maximum(AA3)
M4 = maximum(AA4) ############################

figure("Histograms - average-$nz")
subplot(3,1,1)
PyPlot.hist(100*vec(AA2)./N,bins=100*((0:.5:M2+.5) .- .25)./N)
xticks(100*(0:ceil(M2/10):M2)./N)
ylabel("# of 2-edges")
subplot(3,1,2)
PyPlot.hist(100*vec(AA3)./N,bins=100*((0:.5:M3+.5) .- .25)./N)
xticks(100*(0:ceil(M3/10):M3)./N)
xlabel("Appears in x% of the inferred hypergraphs")
ylabel("# of 3-edges")
subplot(3,1,3) ############################
PyPlot.hist(100*vec(AA4)./N,bins=100*((0:.5:M4+.5) .- .25)./N) ############################
xlabel("Appears in x% of the inferred hypergraphs") ############################
ylabel("# of 4-edges") ############################

figure("Relative error")
PyPlot.plot(vec(mean(relerr,dims=1)),"x")





