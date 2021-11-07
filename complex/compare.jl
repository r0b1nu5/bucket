using LinearAlgebra,PyPlot,DelimitedFiles

include("cycle.jl")
include("L2B.jl")

n = 10
L = cycle(n)
B,w = L2B(L)
A = diagm(0 => diag(L)) - L
γ = 2/π

#θ0 = 2π*(1:n)./n .- π .- π/10 + .3*(2*rand(n).-1)
θs = copy(θ0)
θ = copy(θ0)

z0 = copy(θ0) .+ im*0.
zs = copy(z0)
z = copy(z0)

x0 = exp.(im*θ0)
xs = copy(x0)
x = copy(x0)

T = 50000
h = .001
t = 0
c = 0

while t < T
	global θ,θs,z,zs,x,xs,c
	global t += 1

	if t%1000 == 0
		@info "t = $t"
		c += 1
		writedlm("temp/th_$c.csv",θs[:,1:end-1],',')
		θs = θs[:,end]
		writedlm("temp/zr_$c.csv",real.(zs[:,1:end-1]),',')
		writedlm("temp/zi_$c.csv",imag.(zs[:,1:end-1]),',')
		zs = zs[:,end]
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

	k1 = -γ*B*(sin.(B'*z) - im*cos.(B'*z))
	k2 = -γ*B*(sin.(B'*(z + h*k1/2)) - im*cos.(B'*(z + h*k1/2)))
	k1 = -γ*B*(sin.(B'*(z + h*k2/2)) - im*cos.(B'*(z + h*k2/2)))
	k1 = -γ*B*(sin.(B'*(z + h*k3)) - im*cos.(B'*(z + h*k3)))

	z += h*(k1 + 2*k2 + 2*k3 + k4)/6
	zs = [zs z]
	
	k1 = γ*A*x
	k2 = γ*A*(x + h*k1/2)
	k3 = γ*A*(x + h*k2/2)
	k4 = γ*A*(x + h*k3)

	x += h*(k1 + 2*k2 + 2*k3 + k4)/6
	xs = [xs x]
end

Θs = Matrix{Float64}(undef,n,0)
Zs = Matrix{Complex{Float64}}(undef,n,0)
Xs = Matrix{Complex{Float64}}(undef,n,0)
for d in 1:c
	global Θs,Zs,Xs
	Θs = [Θs readdlm("temp/th_$d.csv",',')]
	rm("temp/th_$d.csv")
	Zs = [Zs (readdlm("temp/zr_$d.csv",',') + im*readdlm("temp/zi_$d.csv",','))]
	rm("temp/zr_$d.csv")
	rm("temp/zi_$d.csv")
	Xs = [Xs (readdlm("temp/xr_$d.csv",',') + im*readdlm("temp/xi_$d.csv",','))]
	rm("temp/xr_$d.csv")
	rm("temp/xi_$d.csv")
end
Θs = [Θs θs]
Zs = [Zs zs]
Xs = [Xs xs]

for i in 1:n
	subplot(1,3,1)
	PyPlot.plot((0:T)*h,Θs[i,:])
	subplot(1,3,2)
	PyPlot.plot((0:T)*h,real.(Zs[i,:]))
	subplot(1,3,3)
	PyPlot.plot((0:T)*h,angle.(Xs[i,:]))
end

subplot(1,3,1)
title("Original Kuramoto")
xlabel("t")
ylabel("θ")
subplot(1,3,2)
title("Complex Kuramoto - 1st")
xlabel("t")
ylabel("real(z)")
subplot(1,3,3)
title("Complex Kuramoto - 2nd")
xlabel("t")
ylabel("arg(x)")


