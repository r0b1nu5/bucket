using DelimitedFiles

include("tools.jl")
include("ksakaguchi.jl")

n = 12
n2 = Int(n/2)
run = rand(1000:9999)

L = cyqle(n)
L[[1,n2+1],[1,n2+1]] += [1. -1.;-1. 1.]
C = [Array(1:n2+1),[Array(n2+1:n);1]]

θ0 = zeros(n)
θ1 = [-2π*Array(0:n2)/(n2+1);-2π*Array(n2-1:-1:1)/(n2+1)]
θ2 = [2π*Array(0:n2)/(n2+1);2π*Array(n2-1:-1:1)/(n2+1)]

 #=
ω = rand(n)
ω .-= mean(ω)
# =#

ω = zeros(n)
ω[[1,n2+1]] = [1,-1]/sqrt(2)

β0 = Array{Float64,2}(undef,2,0)
β1 = Array{Float64,2}(undef,2,0)
β2 = Array{Float64,2}(undef,2,0)

st = 10
αs = Array(LinRange(0,π/2,st))

for α in αs
	@info "================= α = $α ================"
	global β0,β1,β2
	xxx = freq_width(L,ω,θ0,α,C,true)
	β0 = [β0 [xxx[1],xxx[2]]]
	xxx = freq_width(L,ω,θ1,α,C,true)
	β1 = [β1 [xxx[1],xxx[2]]]
	xxx = freq_width(L,ω,θ2,α,C,true)
	β2 = [β2 [xxx[1],xxx[2]]]
end

 #=
writedlm("temp_data/αs_$(run).csv",αs,',')
writedlm("temp_data/β0min_$(run).csv",β0[1,:],',')
writedlm("temp_data/β0max_$(run).csv",β0[2,:],',')
writedlm("temp_data/β1min_$(run).csv",β1[1,:],',')
writedlm("temp_data/β1max_$(run).csv",β1[2,:],',')
# =#

figure()
PyPlot.fill([αs;αs[st:-1:1]],[β0[1,:];β0[2,(st:-1:1)]],color="C0",alpha=.5,label="q=[0,0]")
PyPlot.plot(αs,β0[1,:],color="C0","-o")
PyPlot.fill([αs;αs[st:-1:1]],[β1[1,:];β1[2,(st:-1:1)]],color="C1",alpha=.5,label="q=[1,-1]")
PyPlot.plot(αs,β1[1,:],color="C1","-o")
PyPlot.fill([αs;αs[st:-1:1]],[β2[1,:];β2[2,(st:-1:1)]],color="C2",alpha=.5,label="q=[-1,1]")
PyPlot.plot(αs,β2[1,:],color="C2","-o")
legend()

