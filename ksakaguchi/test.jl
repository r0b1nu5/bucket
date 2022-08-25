using LinearAlgebra, DelimitedFiles

n = 3
ϕ = .3
λ2 = 3.
δ = 2.

B = [1. 1. 0.;-1 0. 1.;0. -1. -1.]
Bt = [1. -1. 0.;1. 0. -1.;0. 1. -1.]
Id = diagm(0 => ones(n))
Π = Id - ones(n,n)./n
v1 = [1.,-1.,0.]./sqrt(2)
v2 = [1.,1.,-2.]./sqrt(6)
R = [v1'; v2']
Rd = pinv(R)

cohes = Vector{Float64}()
avcoh = Vector{Float64}()
μs = Vector{Float64}()

iter = 0
max_iter = 100
c = 0
test = true

γbar = atan(λ2/(δ*tan(ϕ)))

x = rand(n)
xs = Vector{Vector{Float64}}()

while test
	global x,xs,bar,iter,c,cohes,avcoh,μs,Bt,ϕ,Id,R,Rd

	iter += 1

	x = 2π*[rand(),rand(),0.]

	dx = mod.(Bt*x .+ π,2π) .- π

	cohes = maximum(abs.(dx))

	J = cos.(dx .- ϕ).*(1 .- Id)
	J -= diagm(0 => vec(sum(J,dims=2)))
	Jr = R*J*Rd

	μs = maximum((eigvals(Jr + Jr'))/2)
	
	if iter%1000 == 0
		@info "iter: $iter"
	end

	if cohes < γbar && μs > 0.
		push!(xs,x)
	end

	test = (cohes > γbar) || μs < 0.
	
end





