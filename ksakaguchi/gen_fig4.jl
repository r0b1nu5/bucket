using PyPlot,DelimitedFiles

include("ksakaguchi.jl")
include("L2B.jl")

name_cmap = "viridis"
cmap = get_cmap(name_cmap)
col1 = cmap(1/6)
col2 = cmap(3/6)
col3 = cmap(5/6)

 #=
run = 6179

n = 18

L = readdlm("temp_data/L_$run.csv",',')
b,w = L2B(L)
B = [b -b]
Bout = Float64.(B .> 1e-5)

θ0 = vec(readdlm("temp_data/th0_$run.csv",','))
θ1 = vec(readdlm("temp_data/th1_$run.csv",','))
θ2 = vec(readdlm("temp_data/th2_$run.csv",','))

α = .4

ps = [LinRange(0.,.74,40);LinRange(.74,1.15,50);LinRange(1.15,1.62,30);LinRange(1.62,1.72,50)]

ω = [1.;zeros(5);-1.;zeros(11)]

#=
λ0 = Vector{Float64}()
λ1 = Vector{Float64}()
λ2 = Vector{Float64}()
=#

λ0 = Matrix{Float64}(undef,n,0)
λ1 = Matrix{Float64}(undef,n,0)
λ2 = Matrix{Float64}(undef,n,0)

ω0 = Vector{Float64}()
ω1 = Vector{Float64}()
ω2 = Vector{Float64}()

for p in ps
	global L,B,Bout,θ0,θ1,θ2,λ0,λ1,λ2,ω0,ω1,ω2,α,ω

	@info "$(round(100*(p-minimum(ps))/(maximum(ps)-minimum(ps))))%"

	xxx = ksakaguchi(L,p*ω,θ0,α,false,false,.01,1e-10)
	θ0 = vec(xxx[1])
	λ0 = [λ0 sort(real.(eigvals(Bout*diagm(0 => -cos.(vec(B'*θ0 .- α)))*B')))]
#	push!(λ0,sort(real.(eigvals(Bout*diagm(0 => -cos.(vec(B'*θ0 .- α)))*B')))[end-1])
#	push!(λ0,eigvals(jacobian(L,θ0,α))[2])
	push!(ω0,mean(Bout*(sin.(B'*θ0 .- α) .+ sin(α))))

	xxx = ksakaguchi(L,p*ω,θ1,α,false,false,.01,1e-10)
	θ1 = vec(xxx[1])
	λ1 = [λ1 sort(real.(eigvals(Bout*diagm(0 => -cos.(vec(B'*θ1 .- α)))*B')))]
#	push!(λ1,sort(real.(eigvals(Bout*diagm(0 => -cos.(vec(B'*θ1 .- α)))*B')))[end-1])
#	push!(λ1,eigvals(jacobian(L,θ1,α))[2])
	push!(ω1,mean(Bout*(sin.(B'*θ1 .- α) .+ sin(α))))

	xxx = ksakaguchi(L,p*ω,θ2,α,false,false,.01,1e-10)
	θ2 = vec(xxx[1])
	λ2 = [λ2 sort(real.(eigvals(Bout*diagm(0 => -cos.(vec(B'*θ2 .- α)))*B')))]
#	push!(λ2,sort(real.(eigvals(Bout*diagm(0 => -cos.(vec(B'*θ2 .- α)))*B')))[end-1])
	push!(ω2,mean(Bout*(sin.(B'*θ2 .- α) .+ sin(α))))
end

writedlm("temp_data/fig4_ps_$(run).csv",ps,',')
writedlm("temp_data/fig4_l0_$(run).csv",λ0,',')
writedlm("temp_data/fig4_l1_$(run).csv",λ1,',')
writedlm("temp_data/fig4_l2_$(run).csv",λ2,',')
writedlm("temp_data/fig4_w0_$(run).csv",ω0,',')
writedlm("temp_data/fig4_w1_$(run).csv",ω1,',')
writedlm("temp_data/fig4_w2_$(run).csv",ω2,',')

# =#

# #=
run = 6786

n = 18

L = readdlm("temp_data/L_$run.csv",',')
b,w = L2B(L)
B = [b -b]
Bout = Float64.(B .> 1e-5)

θ0 = vec(readdlm("temp_data/th0_$run.csv",','))
θ1 = vec(readdlm("temp_data/th1_$run.csv",','))
θ2 = vec(readdlm("temp_data/th2_$run.csv",','))

α = .4

ps = [LinRange(.2,.3,40);LinRange(.3,1.09,30);LinRange(1.09,1.2,40);LinRange(1.2,1.95,40);LinRange(1.95,2.2,40)]

ω = [1.;zeros(8);-1.;zeros(8)]

#=
λ0 = Vector{Float64}()
λ1 = Vector{Float64}()
λ2 = Vector{Float64}()
=#

λ0 = Matrix{Float64}(undef,n,0)
λ1 = Matrix{Float64}(undef,n,0)
λ2 = Matrix{Float64}(undef,n,0)

ω0 = Vector{Float64}()
ω1 = Vector{Float64}()
ω2 = Vector{Float64}()

for p in ps
	global L,B,Bout,θ0,θ1,θ2,λ0,λ1,λ2,ω0,ω1,ω2,α,ω
1.683265306122449
	@info "$(round(100*(p-minimum(ps))/(maximum(ps)-minimum(ps))))%"

	xxx = ksakaguchi(L,p*ω,θ0,α,false,false,.01,1e-10)
	θ0 = vec(xxx[1])
	λ0 = [λ0 sort(real.(eigvals(Bout*diagm(0 => -cos.(vec(B'*θ0 .- α)))*B')))]
#	push!(λ0,sort(real.(eigvals(Bout*diagm(0 => -cos.(vec(B'*θ0 .- α)))*B')))[end-1])
#	push!(λ0,eigvals(jacobian(L,θ0,α))[2])
	push!(ω0,mean(Bout*(sin.(B'*θ0 .- α) .+ sin(α))))

	xxx = ksakaguchi(L,p*ω,θ1,α,false,false,.01,1e-10)
	θ1 = vec(xxx[1])
	λ1 = [λ1 sort(real.(eigvals(Bout*diagm(0 => -cos.(vec(B'*θ1 .- α)))*B')))]
#	push!(λ1,sort(real.(eigvals(Bout*diagm(0 => -cos.(vec(B'*θ1 .- α)))*B')))[end-1])
#	push!(λ1,eigvals(jacobian(L,θ1,α))[2])
	push!(ω1,mean(Bout*(sin.(B'*θ1 .- α) .+ sin(α))))

	xxx = ksakaguchi(L,p*ω,θ2,α,false,false,.01,1e-10)
	θ2 = vec(xxx[1])
	λ2 = [λ2 sort(real.(eigvals(Bout*diagm(0 => -cos.(vec(B'*θ2 .- α)))*B')))]
#	push!(λ2,sort(real.(eigvals(Bout*diagm(0 => -cos.(vec(B'*θ2 .- α)))*B')))[end-1])
	push!(ω2,mean(Bout*(sin.(B'*θ2 .- α) .+ sin(α))))
end

writedlm("temp_data/fig4_ps_$(run).csv",ps,',')
writedlm("temp_data/fig4_l0_$(run).csv",λ0,',')
writedlm("temp_data/fig4_l1_$(run).csv",λ1,',')
writedlm("temp_data/fig4_l2_$(run).csv",λ2,',')
writedlm("temp_data/fig4_w0_$(run).csv",ω0,',')
writedlm("temp_data/fig4_w1_$(run).csv",ω1,',')
writedlm("temp_data/fig4_w2_$(run).csv",ω2,',')

# =#

