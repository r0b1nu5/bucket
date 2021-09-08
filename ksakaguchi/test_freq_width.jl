using DelimitedFiles

include("tools.jl")
include("ksakaguchi.jl")

n = 18
run = rand(1000:9999)

L = cyqle(n)

θ0 = zeros(n)
θ1 = -2π*Array(1:n)/n
θ2 = 2π*Array(1:n)/n

 #=
ω = rand(n)
ω .-= mean(ω)
# =#

ω = zeros(n)
ω[[1,7]] = [1,-1]/sqrt(2)

β0 = Array{Float64,1}()
β1 = Array{Float64,1}()
β2 = Array{Float64,1}()

αs = Array(LinRange(0,π/2,10))

for α in αs
	global β0,β1,β2
	xxx = freq_width(L,ω,θ0,α,[Array(1:n),],true)
	push!(β0,xxx[1])
	xxx = freq_width(L,ω,θ1,α,[Array(1:n),],true)
	push!(β1,xxx[1])
	xxx = freq_width(L,ω,θ2,α,[Array(1:n),],true)
	push!(β2,xxx[1])
end

 #=
writedlm("temp_data/αs_$(run).csv",αs,',')
writedlm("temp_data/β0min_$(run).csv",β0[1,:],',')
writedlm("temp_data/β0max_$(run).csv",β0[2,:],',')
writedlm("temp_data/β1min_$(run).csv",β1[1,:],',')
writedlm("temp_data/β1max_$(run).csv",β1[2,:],',')
# =#

figure()
PyPlot.fill([αs[1];αs;αs[end]],[0.;β0;0.],color="C0",alpha=.5,label="q=0")
PyPlot.plot(αs,β0,color="C0","-o")
PyPlot.fill([αs[1];αs;αs[end]],[0.;β1;0.],color="C1",alpha=.5,label="q=1")
PyPlot.plot(αs,β1,color="C1","-o")
PyPlot.fill([αs[1];αs;αs[end]],[0.;β2;0.],color="C2",alpha=.5,label="q=-1")
PyPlot.plot(αs,β2,color="C2","-o")
legend()

