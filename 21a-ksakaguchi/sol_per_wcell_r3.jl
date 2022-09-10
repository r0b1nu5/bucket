using DelimitedFiles

include("ksakaguchi.jl")
include("tools.jl")

L = readdlm("ntw_data/ntw10_L.csv",',')
include("ntw_data/ntw10_cycles.jl")

ω = vec(readdlm("temp_data/omega.csv",','))

T1 = 200.
T2 = 350.
dt = .001
dtmax = .0012
dtmin = 1e-8

α = .5
tol = 1e-2

init2 = readdlm("temp_data/init2.csv",',')
sol2 = readdlm("temp_data/sol2.csv",',')
us2 = Int.(readdlm("temp_data/us2.csv",','))

n,N = size(sol2)

init3 = Array{Float64,2}(undef,n,0)
sol3 = Array{Float64,2}(undef,n,0)
us3 = Array{Int64,2}(undef,3,0)

for i in 1:N
	@info "iter = $i/$N"

	global us2,sol2,init2,sol3,us3,init3,dt,dtmax,dtmin

	θ2 = sol2[:,i]

#	local ts,θ = ksakaguchi_ND(L,ω,θ2,α,(T1,T2),RK4())
	local ts,θ = ksakaguchi_ND(L,ω,θ2,α,(T1,T2),(dt,dtmax,dtmin),RK4())
#	local θ,dθ,err,it = ksakaguchi(L,ω,θ2,α,true)

	local θf = mod.(θ[:,end] .- θ[1,end],2π)

	d = 1000
	dc = 0

	while d > tol && dc < size(sol3)[2]
		dc += 1
		d = norm(θf - sol3[:,dc],Inf)
	end

	if d > tol
		init3 = [init3 θ2]
		sol3 = [sol3 θf]
		us3 = [us3 winding(θf,C)]
	end
end

writedlm("temp_data/init3.csv",init3,',')
writedlm("temp_data/sol3.csv",sol3,',')
writedlm("temp_data/us3.csv",us3,',')

