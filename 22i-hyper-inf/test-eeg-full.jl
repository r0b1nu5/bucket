include("hyper_inf.jl")
include("eeg-tools.jl")

nz = 64
λ = .1
ρ = .1
niter = 10

suffix = "666"

# Data to be loaded
subjects = list_all_subjects(109)
#subjects = ["001","002"]
#subjects = ["021",]
states = ["01","02"]
#states = ["02",]
#states = ["03","07","11"]
#states = ["01","02","03","07","11"]

AA2 = zeros(Int64,nz,nz)
AA3 = zeros(Int64,nz,nz,nz)

relerr = zeros(length(states),0)

for subject in subjects
	re = Float64[]
	for state in states
		@info "Running S"*subject*"R"*state

		# Loading data
		file = "eeg-data/S"*subject*"R"*state*".edf"
		s2signal = read_eeg(file)
		sig = zeros(0,length(s2signal["Af3."]))
		for s in sort(collect(keys(s2signal)))
			sig = [sig; s2signal[s]']
		end
		
		# Finite differences
		dt = 1/160
		truncat = findmin(vec(maximum(abs.([sig zeros(nz)]),dims=1)) .> 1e-6)[2] - 1
		
		X0 = denoise_fourier(sig[:,1:truncat],200)
		Y0 = (X0[:,2:end]-X0[:,1:end-1])./dt
		X0 = X0[:,1:end-1]
		X = X0./mean(abs.(X0))
		Y = Y0./mean(abs.(Y0))

		X,ids = restrict_box_size(X,1000)
		Y = Y[:,ids]

		# Inference
		ooi = [2,3]
		dmax = 3
		xxx = hyper_inf(X,Y,ooi,dmax,λ,ρ,niter)
		push!(re,xxx[4])
		
		# Retrieve adjacency tensors	
		A2 = xxx[1][2]
		A3 = xxx[1][3]

		writedlm("eeg-data/S"*subject*"R"*state*"-full-A2-"*suffix*".csv",A2,',')
		writedlm("eeg-data/S"*subject*"R"*state*"-full-A3-"*suffix*".csv",A3,',')
		writedlm("eeg-data/S"*subject*"R"*state*"-full-re-"*suffix*".csv",re,',')
		writedlm("eeg-data/S"*subject*"R"*state*"-full-coeff-"*suffix*".csv",xxx[2],',')
	end
#	global relerr = [relerr re]
end

#writedlm("eeg-data/relative-error-"*suffix*".csv",relerr,',')

#figure("Relative error")
#PyPlot.plot(vec(mean(relerr,dims=1)),"x")





