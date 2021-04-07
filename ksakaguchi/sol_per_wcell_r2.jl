using DelimitedFiles

include("ksakaguchi.jl")
include("tools.jl")

L = readdlm("ntw_data/ntw10_L.csv",',')
include("ntw_data/ntw10_cycles.jl")

ω = readdlm("temp_data/omega.csv",',')

T1 = 50.
T2 = 200.

α = .5
tol = 1e-2

init = readdlm("temp_data/init.csv",',')
sols = readdlm("temp_data/sols.csv",',')
us = Int.(readdlm("temp_data/us.csv",','))

n,N = size(sols)

sol2 = Array{Float64,2}(undef,n,0)
us2 = Array{Int64,2}(undef,3,0)

for i in 1:N
	global us,sols,init,sol2,us2

	θ1 = sols[:,i]

	local ts,θ = ksakaguchi_ND(L,ω,θ1,α,(T1,T2),RK4())

	local θf = mod.(θ[:,end] .- θ[1,end],2π)

	d = 1000
	dc = 0

	while d > tol && dc < size(sol2)[2]
		dc += 1
		d = norm(θf - sol2[:,dc],Inf)
	end

	if d > tol
		sol2 = [sol2 θf]
		us2 = [us winding(θf,C)]
	end
end

writedlm("temp_data/sol2.csv",sol2,',')
writedlm("temp_data/us2.csv",us2,',')

