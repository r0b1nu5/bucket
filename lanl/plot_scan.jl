using PyPlot, DelimitedFiles, LinearAlgebra

param = "amp" # ["amp","fre","pha","tim","ste","noi"]
pmin = .01
pmax = .31
np = 21

n_run = 100
node = 76

ntw = "uk10"
a0 = .2
f0 = .001
p0 = pi/10

T = 100000
dt = .1

n = 120

found = zeros(Bool,n_run,np)
conf1 = zeros(n_run,np)
conf2 = zeros(n_run,np)
conf3 = zeros(n_run,np)

for run in 1:n_run
	for val in 1:np
		ah = readdlm("data/uk10_"*param*"($pmin,$pmax,$np,$val)_run$(run)_ah.csv",',')

		AA = sortslices([abs.(ah) 1:n],dims=1,rev=true)
		PP = AA[:,1]./sum(AA[:,1])

		if Int(AA[1,2]) == node
			found[run,val] = true
		end

		conf1[run,val] = PP[1]
		conf2[run,val] = PP[2]
		conf3[run,val] = PP[3]
	end
end

mf = Array{Float64,1}()
m1 = Array{Float64,1}()
m2 = Array{Float64,1}()
m3 = Array{Float64,1}()
m12 = Array{Float64,1}()

for val in 1:np
	push!(mf,mean(found[:,val]))
	push!(m1,mean(conf1[:,val]))
	push!(m2,mean(conf2[:,val]))
	push!(m3,mean(conf3[:,val]))
	push!(m12,mean(conf1[:,val]-conf2[:,val]))
end

PyPlot.plot(LinRange(pmin,pmax,np),mf)
PyPlot.plot(LinRange(pmin,pmax,np),m1)
PyPlot.plot(LinRange(pmin,pmax,np),m2)
PyPlot.plot(LinRange(pmin,pmax,np),m3)
PyPlot.plot(LinRange(pmin,pmax,np),m12)


# TODO How to plot all this???






