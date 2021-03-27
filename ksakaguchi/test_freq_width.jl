using DelimitedFiles

include("tools.jl")
include("ksakaguchi.jl")

n = 10
run = rand(1000:9999)

L = cyqle(n)

θ0 = zeros(n)
θ1 = -2π*Array(1:n)/n

 #=
ω = rand(n)
ω .-= mean(ω)
# =#

β0 = Array{Float64,2}(undef,2,0)
β1 = Array{Float64,2}(undef,2,0)

αs = [1.,]

for α in αs
	xxx = freq_width(L,ω,θ0,α,Array(1:n),true)
	global β0 = [β0 [xxx[1],xxx[2]]]
	xxx = freq_width(L,ω,θ1,α,Array(1:n),true)
	global β1 = [β1 [xxx[1],xxx[2]]]
end

 #=
writedlm("temp_data/αs_$(run).csv",αs,',')
writedlm("temp_data/β0min_$(run).csv",β0[1,:],',')
writedlm("temp_data/β0max_$(run).csv",β0[2,:],',')
writedlm("temp_data/β1min_$(run).csv",β1[1,:],',')
writedlm("temp_data/β1max_$(run).csv",β1[2,:],',')
# =#

PyPlot.plot(αs,β0[1,:],"o",color="C0")
PyPlot.plot(αs,β0[2,:],"o",color="C0")
PyPlot.plot(αs,β1[1,:],"o",color="C1")
PyPlot.plot(αs,β1[2,:],"o",color="C1")

