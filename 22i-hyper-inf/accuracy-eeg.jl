using PyPlot, DelimitedFiles, Statistics

include("hyper_inf.jl")
include("eeg-tools.jl")

n = 7
type = "avg" # "avg" or "sub"
data_exist = true
Tmax = 1000

# Sensor to zone pairing
s = readdlm("eeg-data/sensors-$n.csv",',',String)
z = readdlm("eeg-data/zones-$n.csv",',',Int64)
s2z = Dict{String,Int64}(s[i] => z[i] for i in 1:length(s))

# Data to be loaded
subjects = list_all_subjects(109)
#subjects = ["001","002"]
states = ["01","02"]

#if data_exist
#	ρ3 = readdlm("eeg-data/violin-data-"*type*"-$n.csv",',')
#else
	ρ3 = zeros(n,0)
#end

dmax = 4
# Defining the basis of functions to use, i.e., the monomials up to order 'dmax'.
@variables x[1:n]
prebasis = polynomial_basis([x[i] for i in 1:n],dmax)
basis = Basis(prebasis,[x[i] for i in 1:n])

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
				Y = (X0[:,2:end]-X0[:,1:end-1])/dt
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
				Y = (X0[:,2:end]-X0[:,1:end-1])/dt
			end
			
			T = size(X)[2]
		end
		
		coeff = readdlm("eeg-data/S"*su*"R"*st*"-"*type*"-coeff.csv",',')

		if !(data_exist)
			αs = zeros(n,0)
			for t in 1:T
				if t%1000 == 0
					writedlm("eeg-data/accuracy-data-"*type*"-$n-S"*su*"R"*st*"-T$t.csv",αs,',')
					αs = zeros(n,0)
				end
				b = basis(X[:,t])
				y = coeff*b
				αs = [αs abs.((y - Y[:,t])./Y[:,t])]
			end
		end
	end
end

αs = zeros(n,0)
for su in subjects
	for st in states
		for t in 1000:1000:Tmax
			global αs = [αs readdlm("eeg-data/accuracy-data-"*type*"-$n-S"*su*"R"*st*"-T$t.csv",',')]
		end
	end
end
PyPlot.hist(vec(αs),100)









