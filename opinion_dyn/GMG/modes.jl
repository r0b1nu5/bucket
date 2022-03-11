using PyPlot, LinearAlgebra

include("distrib.jl")

n = 201

epss = LinRange(0,2,100)

xs = sort(bi_gauss(n,-.5,.2,.5,.2))

M1 = Array{Float64,2}(undef,n,0)
M2 = Array{Float64,2}(undef,n,0)
M3 = Array{Float64,2}(undef,n,0)
M4 = Array{Float64,2}(undef,n,0)
M5 = Array{Float64,2}(undef,n,0)

l2 = Array{Float64,1}()
l3 = Array{Float64,1}()
l4 = Array{Float64,1}()
l5 = Array{Float64,1}()
l6 = Array{Float64,1}()

for eps in epss
	global M1,M2,M3,M4,M5,l2,l3,l4,l5,l6,xs,n

	@info "$(round(100*eps/maximum(epss)))%"

	A = Float64.((0 .< abs.(repeat(xs,1,n) - repeat(xs',n,1)) .< eps))
	L = diagm(0 => vec(sum(A,dims=1))) - A

	ei = eigen(L)

	id = 1

#	while ei.values[id] < 1e-8 && id < n-3
#		id += 1
#	end
	id += 1

	M1 = [M1 abs.(ei.vectors[:,id])]
	M2 = [M2 abs.(ei.vectors[:,id+1])]
	M3 = [M3 abs.(ei.vectors[:,id+2])]
	M4 = [M4 abs.(ei.vectors[:,id+3])]
	M5 = [M5 abs.(ei.vectors[:,id+4])]

	push!(l2,ei.values[id])
	push!(l3,ei.values[id+1])
	push!(l4,ei.values[id+2])
	push!(l5,ei.values[id+3])
	push!(l6,ei.values[id+4])
end

figure(123)

subplot(1,5,1)
contourf(log.(M1 .+ 1e-10),100)

subplot(1,5,2)
contourf(log.(M2 .+ 1e-10),100)

subplot(1,5,3)
contourf(log.(M3 .+ 1e-10),100)

subplot(1,5,4)
contourf(log.(M4 .+ 1e-10),100)

subplot(1,5,5)
contourf(log.(M5 .+ 1e-10),100)
colorbar()

figure(124)
contourf(log.(M1 .+ 1e-10),100)
colorbar()

# #=
figure(321)
PyPlot.plot(epss,l2)
PyPlot.plot(epss,l3)
PyPlot.plot(epss,l4)
PyPlot.plot(epss,l5)
PyPlot.plot(epss,l6)
# =#

