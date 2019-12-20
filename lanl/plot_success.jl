using PyPlot, DelimitedFiles, LinearAlgebra

n_run = 100

ntw = "uk10"
a0 = .2
f0 = .001
p0 = pi/10

T = 100000
dt = .1

n = 120

found = zeros(Bool,n_run,n)
conf1 = zeros(n_run,n)
conf2 = zeros(n_run,n)
conf3 = zeros(n_run,n)

for run in 1:n_run
	for node in 1:n
		ah = readdlm("data/uk10_$(node)_$(run)_ah.csv",',')

		AA = sortslices([abs.(ah) 1:n],dims=1,rev=true)
		PP = AA[:,1]./sum(AA[:,1])

		if Int(AA[1,2]) == node
			found[run,node] = true
		end

		conf1[run,node] = PP[1]
		conf2[run,node] = PP[2]
		conf3[run,node] = PP[3]
	end
end

# TODO How to plot all this???







