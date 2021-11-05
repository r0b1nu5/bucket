using LinearAlgebra,PyPlot,DelimitedFiles

include("cycle.jl")
include("L2B.jl")

n = 10
L = cycle(n)
B,w = L2B(L)
A = diagm(0 => diag(L)) - L
γ = 2/π

θ0 = 2π*rand(n) .- π
θ0 = 2π*(1:n)./n .- π + .1*rand(n)
θs = copy(θ0)
θ = copy(θ0)

x0 = exp.(im*θ0)
xs = copy(x0)
x = copy(x0)

T = 100000
h = .001
t = 0
c = 0

while t < T
	global θ,θs,x,xs,c
	global t += 1

	if t%1000 == 0
		@info "t = $t"
		c += 1
		writedlm("temp/th_$c.csv",θs[:,1:end-1],',')
		θs = θs[:,end]
		writedlm("temp/xr_$c.csv",real.(xs[:,1:end-1]),',')
		writedlm("temp/xi_$c.csv",imag.(xs[:,1:end-1]),',')
		xs = xs[:,end]
	end

	k1 = -B*sin.(B'*θ)
	k2 = -B*sin.(B'*(θ + h*k1/2))
	k3 = -B*sin.(B'*(θ + h*k2/2))
	k4 = -B*sin.(B'*(θ + h*k3))

	θ += h*(k1 + 2*k2 + 2*k3 + k4)/6
	θs = [θs (mod.(θ .+ π,2π) .- π)]

	k1 = γ*A*x
	k2 = γ*A*(x + h*k1/2)
	k3 = γ*A*(x + h*k2/2)
	k4 = γ*A*(x + h*k3)

	x += h*(k1 + 2*k2 + 2*k3 + k4)/6
	xs = [xs x]
end

Θs = Matrix{Float64}(undef,n,0)
Xs = Matrix{Complex{Float64}}(undef,n,0)
for d in 1:c
	global Θs,Xs
	Θs = [Θs readdlm("temp/th_$d.csv",',')]
	rm("temp/th_$d.csv")
	Xs = [Xs (readdlm("temp/xr_$d.csv",',') + im*readdlm("temp/xi_$d.csv",','))]
	rm("temp/xr_$d.csv")
	rm("temp/xi_$d.csv")
end
Θs = [Θs θs]
Xs = [Xs xs]

for i in 1:n
	subplot(1,2,1)
	PyPlot.plot((0:T)*h,Θs[i,:])
	subplot(1,2,2)
	PyPlot.plot((0:T)*h,angle.(Xs[i,:]))
end


