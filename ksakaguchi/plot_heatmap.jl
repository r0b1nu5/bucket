using PyPlot

res = 100
cm = "RdBu"

ϕ = .3
R = [1/sqrt(2) -1/sqrt(2) 0.;1/sqrt(6) 1/sqrt(6) -2/sqrt(6)]
λ2 = 3.
dmax = 2.
γbar = atan(λ2/(dmax*tan(ϕ)))

xs = LinRange(-π-.1,π+.1,res)
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

PyPlot.contourf(X,Y,μ',50,cmap=cm)
PyPlot.plot([-π,0,-π,-π],[0,π,π,0],"k")
PyPlot.plot([0,π,π,0],[-π,-π,0,-π],"k")
#PyPlot.plot([minimum(xs),maximum(xs)],[π,π],"k")
#PyPlot.plot([minimum(xs),maximum(xs)],[-π,-π],"k")
#PyPlot.plot([π,π],[minimum(xs),maximum(xs)],"k")
#PyPlot.plot([-π,-π],[minimum(xs),maximum(xs)],"k")
PyPlot.plot([-π,0,π,π,0,-π,-π],[-π,-π,0,π,π,0,-π],"k")
PyPlot.plot(γbar/π*[-π,0,π,π,0,-π,-π],γbar/π*[-π,-π,0,π,π,0,-π],"--k")

colorbar(label="μ")
xlabel("x1")
ylabel("x2")




