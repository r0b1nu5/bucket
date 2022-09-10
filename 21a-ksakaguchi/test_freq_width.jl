using DelimitedFiles

include("freq_width.jl")

n = 18
run = rand(1000:9999)

L = cyqle(n)

C = [Array(1:n),]

θ0 = zeros(n)
θ1 = -2π*Array(1:n)/n
θ2 = 2π*Array(1:n)/n

 #=
ω = rand(n)
ω .-= mean(ω)
# =#

 #=
ω = zeros(n)
ω[[1,7]] = [1,-1]/sqrt(2)
# =#

# #=
ω = zeros(n)
ω[1] = 1.
ω .-= mean(ω)
ω /= norm(ω)
# =#

β0 = Array{Float64,2}(undef,3,0)
β1 = Array{Float64,2}(undef,3,0)
β2 = Array{Float64,2}(undef,3,0)

st = 10
αs = Array(LinRange(0,π/2,st))

for α in αs
	@info "================= α = $α ================"
	global β0,β1,β2
	xxx = freq_width(L,ω,θ0,α,C,true)
	β0 = [β0 [xxx[1],xxx[2],xxx[3]]]
	xxx = freq_width(L,ω,θ1,α,C,true)
	β1 = [β1 [xxx[1],xxx[2],xxx[3]]]
	xxx = freq_width(L,ω,θ2,α,C,true)
	β2 = [β2 [xxx[1],xxx[2],xxx[3]]]
end

 #=
writedlm("temp_data/αs_$(run).csv",αs,',')
writedlm("temp_data/β0_$(run).csv",β0,',')
writedlm("temp_data/β1_$(run).csv",β1,',')
writedlm("temp_data/β2_$(run).csv",β2,',')
writedlm("temp_data/L_$(run).csv",L,',')
writedlm("temp_data/w_$(run).csv",ω,',')
writedlm("temp_data/th0_$(run).csv",θ0,',')
writedlm("temp_data/th1_$(run).csv",θ1,',')
writedlm("temp_data/th2_$(run).csv",θ2,',')
# =#

figure()
PyPlot.fill([αs;αs[st:-1:1]],[β0[1,:];β0[3,(st:-1:1)]],color="C0",alpha=.5,label="q=[0,0]")
PyPlot.plot(αs,β0[1,:],color="C0","-o")
PyPlot.plot(αs,β0[2,:],color="C0","-")
PyPlot.fill([αs;αs[st:-1:1]],[β1[1,:];β1[3,(st:-1:1)]],color="C1",alpha=.5,label="q=[1,-1]")
PyPlot.plot(αs,β1[1,:],color="C1","-o")
PyPlot.plot(αs,β1[2,:],color="C1","-")
PyPlot.fill([αs;αs[st:-1:1]],[β2[1,:];β2[3,(st:-1:1)]],color="C2",alpha=.5,label="q=[-1,1]")
PyPlot.plot(αs,β2[1,:],color="C2","-o")
PyPlot.plot(αs,β2[2,:],color="C2","-")
legend()



