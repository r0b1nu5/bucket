using PyPlot

include("iterations.jl")
include("ntw_data/ntw10_cycles.jl")

L = readdlm("ntw_data/ntw10_L.csv",',')
b,w = L2B(L)
B = [b -b]
Bout = B.*(B .> 0)
n,m = size(B)
m2 = Int(m/2)
Id = diagm(0 => ones(m))

ω = vec(readdlm("ntw_data/ntw10_om1.csv",','))

P = dir_cycle_proj(B,ones(m))

α = .1
γ = π/2 - α

hγm = h(-γ)
hγp = h(γ)
hγ = (hγm,hγp)

#u = [0,-1,0]
u = [0,0,0]

n_ini = 1000

F1 = (hγp - hγm)*rand(m,n_ini) .+ hγm

T11 = Array{Float64,2}(undef,m,0)
T12 = Array{Float64,2}(undef,m,0)
T13 = Array{Float64,2}(undef,m,0)
T14 = Array{Float64,2}(undef,m,0)

for i in 1:n_ini
	global T11 = [T11 Tu_dir(F1[:,i],u,P,C,.1*Id)]
	global T12 = [T12 Tu_dir(F1[:,i],u,P,C,.2*Id)]
	global T13 = [T13 Tu_dir(F1[:,i],u,P,C,.5*Id)]
	global T14 = [T14 Tu_dir(F1[:,i],u,P,C,Id)]
end

k = 1

figure()
PyPlot.plot(F1[k,:],F1[k+m2,:],".")
PyPlot.plot(T11[k,:],T11[k+m2,:],".")
PyPlot.plot(T12[k,:],T12[k+m2,:],".")
PyPlot.plot(T13[k,:],T13[k+m2,:],".")
PyPlot.plot(T14[k,:],T14[k+m2,:],".")


F2 = Array{Float64,2}(undef,m,0)
T21 = Array{Float64,2}(undef,m,0)
T22 = Array{Float64,2}(undef,m,0)
T23 = Array{Float64,2}(undef,m,0)
T24 = Array{Float64,2}(undef,m,0)
for i in 1:n_ini
	global F2 = [F2 rand_init(ω,Bout,hγ)]
	global T21 = [T21 Tu_dir(F2[:,i],u,P,C,.1*Id)]
	global T22 = [T22 Tu_dir(F2[:,i],u,P,C,.2*Id)]
	global T23 = [T23 Tu_dir(F2[:,i],u,P,C,.5*Id)]
	global T24 = [T24 Tu_dir(F2[:,i],u,P,C,Id)]
end

figure()
PyPlot.plot(F2[k,:],F2[k+m2,:],".")
PyPlot.plot(T21[k,:],T21[k+m2,:],".")
PyPlot.plot(T22[k,:],T22[k+m2,:],".")
PyPlot.plot(T23[k,:],T23[k+m2,:],".")
PyPlot.plot(T24[k,:],T24[k+m2,:],".")

F3 = Array{Float64,2}(undef,m,0)
T31 = Array{Float64,2}(undef,m,0)
T32 = Array{Float64,2}(undef,m,0)
T33 = Array{Float64,2}(undef,m,0)
T34 = Array{Float64,2}(undef,m,0)

nomi = sqrt(minimum(sum((F2[[k,k+m2],:]).^2,dims=1)))
for i in 1:n_ini
	no = norm([F2[k,i],F2[k+m2,i]])
	global F3 = [F3 F2[:,i]*nomi/no]
	global T31 = [T31 T21[:,i]*nomi/no]
	global T32 = [T32 T22[:,i]*nomi/no]
	global T33 = [T33 T23[:,i]*nomi/no]
	global T34 = [T34 T24[:,i]*nomi/no]
end

figure()
PyPlot.plot(F3[k,:],F3[k+m2,:],".")
PyPlot.plot(T31[k,:],T31[k+m2,:],".")
PyPlot.plot(T32[k,:],T32[k+m2,:],".")
PyPlot.plot(T33[k,:],T33[k+m2,:],".")
PyPlot.plot(T34[k,:],T34[k+m2,:],".")


F4 = copy(F3)
T41 = Array{Float64,2}(undef,m,0)
T42 = Array{Float64,2}(undef,m,0)
T43 = Array{Float64,2}(undef,m,0)
T44 = Array{Float64,2}(undef,m,0)
for i in 1:n_ini
	global T41 = [T41 Tu_dir(F4[:,i],u,P,C,.1*Id)]
	global T42 = [T42 Tu_dir(T41[:,i],u,P,C,.1*Id)]
	global T43 = [T43 Tu_dir(T42[:,i],u,P,C,.1*Id)]
	global T44 = [T44 Tu_dir(T43[:,i],u,P,C,.1*Id)]
end

figure()
PyPlot.plot(F4[k,:],F4[k+m2,:],".")
PyPlot.plot(T41[k,:],T41[k+m2,:],".")
PyPlot.plot(T42[k,:],T42[k+m2,:],".")
PyPlot.plot(T43[k,:],T43[k+m2,:],".")
PyPlot.plot(T44[k,:],T44[k+m2,:],".")



