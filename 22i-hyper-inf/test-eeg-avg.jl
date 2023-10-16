include("hyper_inf.jl")
include("eeg-tools.jl")

nz = 7

# Sensor to zone pairing
s = readdlm("eeg-data/sensors-$nz.csv",',',String)
z = readdlm("eeg-data/zones-$nz.csv",',',Int64)
s2z = Dict{String,Int64}(s[i] => z[i] for i in 1:length(s))

# Loading data
subject = "003"
state = "02"
file = "eeg-data/S"*subject*"R"*state*".edf"

s2signal = read_eeg(file)
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

# Plots
figure("Histograms - average-7 - S"*subject*"R"*state)
A2 = inferred_adj_2nd(xxx[1][2],nz)[2]
subplot(2,1,1)
PyPlot.hist(vec(abs.(A2)),20)
A3 = inferred_adj_3rd(xxx[1][3],nz)[2]
subplot(2,1,2)
PyPlot.hist(vec(abs.(A3)),20)






