include("hyper_inf.jl")
include("eeg-tools.jl")

n = 7

# Sensor selection
slist = readdlm("eeg-data/slist-$n.csv",',',String)

# Loading data
subject = "001"
state = "03"
file = "eeg-data/S"*subject*"R"*state*".edf"

s2signal = read_eeg(file)

# Finite differences
dt = 1/160
truncate = 128

X0 = zeros(0,length(s2signal["Cp3."])-truncate)
for sen in slist
	global X0 = [X0;s2signal[sen][1:end-truncate]']
end
X = X0[:,1:end-1]
Y = (X0[:,2:end]-X0[:,1:end-1])/dt

# Inference
ooi = [2,3]
dmax = 4
xxx = hyper_inf(X,Y,ooi,dmax)

# Plots
figure("Histograms - subsampling-7 - S"*subject*"R"*state)
A2 = inferred_adj_2nd(xxx[1][2],nz)[2]
subplot(2,1,1)
PyPlot.hist(vec(abs.(A2)),20)
A3 = inferred_adj_3rd(xxx[1][3],nz)[2]
subplot(2,1,2)
PyPlot.hist(vec(abs.(A3)),20)






