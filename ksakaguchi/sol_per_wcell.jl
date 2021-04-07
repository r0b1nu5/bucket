using DelimitedFiles

include("ksakaguchi.jl")
include("tools.jl")

L = readdlm("ntw_data/ntw10_L.csv",',')
include("ntw_data/ntw10_cycles.jl")

n = size(L)[1]

ω = rand(n)
ω .-= mean(ω)
ω *= .1

T = 50.
α = .5
#α = 0.
tol = 1e-2

max_N = 30
max_iter = 10

init = Array{Float64,2}(undef,n,0)
sols = Array{Float64,2}(undef,n,0)
us = Array{Int64,2}(undef,3,0)
iter = 0

while size(sols)[2] < max_N && iter < max_iter
	global us,sols,init
	global iter += 1

	local θ0 = 2π*rand(n)

	local ts,θ = ksakaguchi_ND(L,ω,θ0,α,(0.,T),RK4())
#	local θ,dθ,err,it = ksakaguchi(L,ω,θ0,α,true)

	local θf = mod.(θ[:,end] .- θ[1,end],2π)

	d = 1000.
	dc = 0

	while d > tol && dc < size(sols)[2]
		dc += 1
		d = norm(θf - sols[:,dc],Inf)
	end

	if d > tol
		init = [init θ0]
		sols = [sols θf]
		us = [us winding(θf,C)]
	end

	if iter%100 == 0
		writedlm("temp_data/iter.csv",iter,',')
	end
end

writedlm("temp_data/omega.csv",ω,',')
writedlm("temp_data/sols.csv",sols,',')
writedlm("temp_data/init.csv",init,',')
writedlm("temp_data/us.csv",us,',')

