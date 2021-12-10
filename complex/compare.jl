using LinearAlgebra,PyPlot,DelimitedFiles

include("cycle.jl")
include("L2B.jl")
include("kuramoto.jl")

n = 10
L = cycle(n)
B,w = L2B(L)
A = diagm(0 => diag(L)) - L
γ = 1.
win = 100
q = 4
α = .01

θ0 = 2π*q*(1:n)./n + α*(2*rand(n).-1)
#θ0 = vec(readdlm("th1.csv",','))
θ0 = vec(readdlm("th2.csv")) + α*(2*rand(n) .-1)
θs = copy(θ0)
θ = copy(θ0)

z0 = copy(θ0) .+ im*0.
zs = copy(z0)
z = copy(z0)

x0 = exp.(im*θ0)
xs = copy(x0)
x = copy(x0)

T1 = 10000
T2 = 50000
h = .0005
t = 0
c = 10

while t < T1
	global θ,θs,z,zs,x,xs,c
	global t += 1000

	θs = kuramoto(L,zeros(n),θs,false,h,1000,-1.)

	zs = complex_k1(L,zeros(n),zs,win,false,h,1000,-1.)

	xs = exp.(im*zs)
	
	if t%1000 == 0
		@info "t = $t"
		c += 1
		writedlm("temp/th_$c.csv",θs[:,1:end-1],',')
		θs = θs[:,end]
		writedlm("temp/zr_$c.csv",real.(zs[:,1:end-1]),',')
		writedlm("temp/zi_$c.csv",imag.(zs[:,1:end-1]),',')
		zs = real.(zs[:,end]) .+ im*0.
		writedlm("temp/xr_$c.csv",real.(xs[:,1:end-1]),',')
		writedlm("temp/xi_$c.csv",imag.(xs[:,1:end-1]),',')
		xs = xs[:,end]
	end

end

η = copy(real.(zs))
ηs = copy(η)
c2 = copy(c)

while t < T2
	global θ,θs,η,ηs,z,zs,x,xs,c
	global t += 1000

	θs = kuramoto(L,zeros(n),θs,false,h,1000,-1.)
	ηs = kuramoto(L,zeros(n),ηs,false,h,1000,-1.)

	zs = complex_k1(L,zeros(n),zs,win,false,h,1000,-1.)

	xs = exp.(im*zs)
	
	if t%1000 == 0
		@info "t = $t"
		c += 1
		writedlm("temp/th_$c.csv",θs[:,1:end-1],',')
		θs = θs[:,end]
		writedlm("temp/et_$c.csv",ηs[:,1:end-1],',')
		ηs = ηs[:,end]
		writedlm("temp/zr_$c.csv",real.(zs[:,1:end-1]),',')
		writedlm("temp/zi_$c.csv",imag.(zs[:,1:end-1]),',')
		zs = real.(zs[:,end]) .+ im*0.
		writedlm("temp/xr_$c.csv",real.(xs[:,1:end-1]),',')
		writedlm("temp/xi_$c.csv",imag.(xs[:,1:end-1]),',')
		xs = xs[:,end]
	end

end

Θs = Matrix{Float64}(undef,n,0)
Hs = Matrix{Float64}(undef,n,0)
Zs = Matrix{Complex{Float64}}(undef,n,0)
Xs = Matrix{Complex{Float64}}(undef,n,0)
for d in 11:c2
	global Θs,Zs,Xs
	Θs = [Θs (mod.(readdlm("temp/th_$d.csv",',') .+ π,2π) .- π)]
	rm("temp/th_$d.csv")
	Zs = [Zs (readdlm("temp/zr_$d.csv",',') + im*readdlm("temp/zi_$d.csv",','))]
	rm("temp/zr_$d.csv")
	rm("temp/zi_$d.csv")
	Xs = [Xs (readdlm("temp/xr_$d.csv",',') + im*readdlm("temp/xi_$d.csv",','))]
	rm("temp/xr_$d.csv")
	rm("temp/xi_$d.csv")
end
for d in c2+1:c
	global Θs,Hs,Zs,Xs
	Θs = [Θs (mod.(readdlm("temp/th_$d.csv",',') .+ π,2π) .- π)]
	rm("temp/th_$d.csv")
	Hs = [Hs (mod.(readdlm("temp/et_$d.csv",',') .+ π,2π) .- π)]
	rm("temp/et_$d.csv")
	Zs = [Zs (readdlm("temp/zr_$d.csv",',') + im*readdlm("temp/zi_$d.csv",','))]
	rm("temp/zr_$d.csv")
	rm("temp/zi_$d.csv")
	Xs = [Xs (readdlm("temp/xr_$d.csv",',') + im*readdlm("temp/xi_$d.csv",','))]
	rm("temp/xr_$d.csv")
	rm("temp/xi_$d.csv")
end
Θs = [Θs (mod.(θs .+ π,2π) .- π)]
Hs = [Hs (mod.(ηs .+ π,2π) .- π)]
Zs = [Zs zs]
Xs = [Xs xs]

for i in 1:n
	AA = Θs[i,:]
	AAA = Hs[i,:]
	BB = mod.(real.(Zs[i,:]) .+ π,2π) .- π
	CC = angle.(Xs[i,:])

	figure(1)
	subplot(2,3,1)
	PyPlot.plot((0:T2)*h,AA)
	subplot(2,3,2)
	PyPlot.plot((0:T2)*h,BB)
	subplot(2,3,3)
	PyPlot.plot((0:T2)*h,imag.(Zs[i,:]))
#	PyPlot.plot((0:T)*h,CC)
	subplot(2,3,4)
	PyPlot.plot((T1:T2)*h,AAA)
#	PyPlot.plot((0:T)*h,abs.(AA - BB))
	subplot(2,3,5)
	PyPlot.plot((0:T2)*h,abs.(AA - CC))
	subplot(2,3,6)
	PyPlot.plot((0:T2)*h,abs.(BB - CC))

	figure(2)
	PyPlot.plot((0:T2)*h,AA,color="C$(i-1)")
	PyPlot.plot((T1:T2)*h,AAA,"--",color="C$(i-1)")
	PyPlot.plot((0:T2)*h,BB,":",color="C$(i-1)")
end

figure(1)
subplot(2,3,1)
title("Original Kuramoto")
xlabel("t[s]")
ylabel("θ")
subplot(2,3,2)
title("Complex Kuramoto - 1st")
xlabel("t[s]")
ylabel("real(z)")
subplot(2,3,3)
title("Complex Kuramoto - 2nd")
xlabel("t[s]")
ylabel("arg(x)")
subplot(2,3,4)
xlabel("t[s]")
ylabel("θ - real(z)")
subplot(2,3,5)
xlabel("t[s]")
ylabel("θ - arg(x)")
subplot(2,3,6)
xlabel("t[s]")
ylabel("real(z) - arg(x)")




