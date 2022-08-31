using PyPlot, LinearAlgebra

res = 100
cm = "RdBu"

ϕ = .3
R = [1/sqrt(2) -1/sqrt(2) 0.;1/sqrt(6) 1/sqrt(6) -2/sqrt(6)]
λ2 = 3.
dmax = 2.
γbar = atan(λ2/(dmax*tan(ϕ)))

xmin = -π
xmax = π
xs = LinRange(xmin,xmax,res)
X = repeat(Vector(xs)',res,1)
Y = repeat(Vector(xs),1,res)

μ = zeros(res,res)
for i in 1:res
	for j in 1:res
		x1 = xs[i]
		x2 = xs[j]
		J = [(-cos(x1-x2-ϕ)-cos(x1-ϕ)) cos(x1-x2-ϕ) cos(x1-ϕ);
		     cos(x2-x1-ϕ) (-cos(x2-x1-ϕ)-cos(x2-ϕ)) cos(x2-ϕ);
		     cos(-x1-ϕ) cos(-x2-ϕ) (-cos(-x1-ϕ)-cos(-x2-ϕ))]
		μ[i,j] = maximum(eigvals(R*(J+J')*R'/2))
	end
end

figure("ϕ = $ϕ")

ks = 0:1
ls = 0:1
for k in ks
	for l in ls
		PyPlot.contourf(X .+ 2π*k,Y .+ 2π*l,μ',50,cmap=cm)
		PyPlot.plot([-π,0,-π,-π] .+ 2π*k,[0,π,π,0] .+ 2π*l,"k")
		PyPlot.plot([0,π,π,0] .+ 2π*k,[-π,-π,0,-π] .+ 2π*l,"k")
		#PyPlot.plot([minimum(xs),maximum(xs)] .+ 2π*k,[π,π] .+ 2π*l,"k")
		#PyPlot.plot([minimum(xs),maximum(xs)] .+ 2π*k,[-π,-π] .+ 2π*l,"k")
		#PyPlot.plot([π,π] .+ 2π*k,[minimum(xs),maximum(xs)] .+ 2π*l,"k")
		#PyPlot.plot([-π,-π] .+ 2π*k,[minimum(xs),maximum(xs)] .+ 2π*l,"k")
		PyPlot.plot([-π,0,π,π,0,-π,-π] .+ 2π*k,[-π,-π,0,π,π,0,-π] .+ 2π*l,"k")
		PyPlot.plot(γbar/π*[-π,0,π,π,0,-π,-π] .+ 2π*k,γbar/π*[-π,-π,0,π,π,0,-π] .+ 2π*l,"--k")
	end
end

colorbar(label="μ")
xlabel("x1")
ylabel("x2")




