using LinearAlgebra, DelimitedFiles

include("tools.jl")

ntw = "cyc5"
#ntw = "ntw20"

L = readdlm("data/"*ntw*"_L.csv",',')
n = size(L)[1]
λ2 = eigvals(L)[2]
δ = maximum(diag(L))
ϕ = .3
γbar = atan(λ2/(δ*tan(ϕ)))

ampl = π

B,w = L2B(L)
Bt = B'

Id = diagm(0 => ones(n))
Π = Id - ones(n,n)./n

R = gen_R(n)
Rt = R'

cohes = Vector{Float64}()
avcoh = Vector{Float64}()
μs = Vector{Float64}()
excepx = zeros(n,0)
excepi = Int64[]

iter = 0
max_iter = 500000
c = 0

while iter < max_iter
	global iter,c,cohes,avcoh,μs,Bt,ϕ,Id,R,Rt,excepx,excepi

	iter += 1

	x = ampl*[rand(n-1);0.]

	dx = mod.(Bt*x .+ π,2π) .- π

	push!(cohes,maximum(abs.(dx)))
	push!(avcoh,sum(abs.(dx)))

	J = -L.*cos.(x*ones(1,n) - ones(n)*x'.- ϕ).*(1 .- Id)
	J -= diagm(0 => vec(sum(J,dims=2)))
	Jr = R*J*Rt

	push!(μs,maximum(eigvals(Jr + Jr'))/2)
	
	if (cohes[end] < γbar) && (μs[end] > 0.)
		excepx = [excepx x]
		push!(excepi,iter)
	end

	if iter%1000 == 0
		c += 1
		writedlm("temp_data/"*ntw*"_$c.csv",[cohes avcoh μs],',')
		cohes = Vector{Float64}()
		avcoh = Vector{Float64}()
		μs = Vector{Float64}()
	end
	
end

c += 1
writedlm("temp_data/"*ntw*"_$c.csv",[cohes avcoh μs],',')
writedlm("temp_data/"*ntw*"_exceptionsx.csv",excepx,',')
writedlm("temp_data/"*ntw*"_exceptionsi.csv",excepi,',')




