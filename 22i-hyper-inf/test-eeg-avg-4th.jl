include("hyper_inf.jl")
include("eeg-tools.jl")

nz = 7
λ = .1
ρ = .1

suffix = "77x"

# Data to be loaded
subjects = list_all_subjects(109)
#subjects = ["001","002"]
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
		ooi = [2,3,4,5,6,7]
		dmax = 7
		xxx = hyper_inf(X,Y,ooi,dmax,λ,ρ)
		push!(re,xxx[4])
		
		# Retrieve adjacency tensors	
		A2 = xxx[1][2]
		A3 = maximum(ooi) > 2 ? xxx[1][3] : zeros(0,4)
		A4 = maximum(ooi) > 3 ? xxx[1][4] : zeros(0,5) ############################

		writedlm("eeg-data/S"*subject*"R"*state*"-avg-A2-"*suffix*".csv",A2,',')
		writedlm("eeg-data/S"*subject*"R"*state*"-avg-A3-"*suffix*".csv",A3,',')
		writedlm("eeg-data/S"*subject*"R"*state*"-avg-A4-"*suffix*".csv",A4,',') ############################

		writedlm("eeg-data/S"*subject*"R"*state*"-avg-coeff-"*suffix*".csv",xxx[2],',')	
	end
	global relerr = [relerr re]
end

writedlm("eeg-data/relative-error-"*suffix*".csv",relerr,',')

figure("Relative error")
PyPlot.plot(vec(mean(relerr,dims=1)),"x")





