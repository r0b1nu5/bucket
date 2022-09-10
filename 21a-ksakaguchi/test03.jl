using PyPlot, DelimitedFiles

include("tools.jl")

#L = readdlm("ntw_data/ntw10_L.csv",',')
L = readdlm("ntw_data/ntw9_L.csv",',')
b,w = L2B(L)

B = [b -b]
Bout = B.*(B .> 0.)

n,m = size(B)
Id = diagm(0 => ones(m))

D = diagm(0 => rand(m))
Δ = Id - diagm(0 => rand(m))
on = opnorm(Δ)
x = Array{Float64,1}()
y = Array{Float64,1}()

for i in 1:1000
	d = rand(Int(m/2))
	D = diagm(0 => [d;d])
P = Id - D*B'*pinv(Bout*D*B')*Bout

α = 1.
Q = (1+α)*P'*P + α*(Id - P - P')
Qh = sqrt(Q)
Qi = inv(Qh)

	push!(x,opnorm(Qh*Δ*Qi))
end

PyPlot.plot(c*ones(1000),x.-1,".")

c+=1

