using PyPlot, LinearAlgebra, DelimitedFiles

include("tools.jl")
include("L2B.jl")

cm = "RdBu"

n = 13
A = readdlm("ntw_data/ntw13_A.csv",',')
L = readdlm("ntw_data/ntw13_L.csv",',')
B,w = L2B(L)
R = zeros(n-1,n)
for i in 1:n-1
	for j in 1:i
		R[i,j] = 1/sqrt(i + i^2)
	end
	R[i,i+1] = -i/sqrt(i + i^2)
end
Σ = [[1,2,3,4,5,6,7,1],[7,8,9,10,11,12,13,1,7]]

#ϕ = .3
ϕ = .05
#ϕ = .01
λ2 = minimum(eigvals(R*L*R'))
dmax = maximum(diag(L))
γbar = atan(λ2/(dmax*tan(ϕ)))

res = 1000

θ0 = zeros(n)
θ1 = vec(readdlm("ntw_data/ntw13_$(ϕ)_th1.csv",','))
θ2 = vec(readdlm("ntw_data/ntw13_$(ϕ)_th2.csv",','))
θ3 = vec(readdlm("ntw_data/ntw13_$(ϕ)_th3.csv",','))
θ4 = vec(readdlm("ntw_data/ntw13_$(ϕ)_th4.csv",','))

v1 = θ1 - θ0
v1 ./= norm(v1)
v2 = θ2 - θ0
v2 -= dot(v1,v2)*v1
v2 ./= norm(v2)
V = [v1 v2]

R = gen_R(n)
P = pinv(R*V)*R

x1 = [v1';v2']*θ1
x2 = [v1';v2']*θ2
#x1 = P*θ1
#x2 = P*θ2
#x3 = P*θ3
#x4 = P*θ4

xmin = min(0.,x1[1],x2[1])
xmax = max(0.,x1[1],x2[1])
ymin = min(0.,x1[2],x2[2])
ymax = max(0.,x1[2],x2[2])
dx = xmax-xmin
dy = ymax-ymin

ampl = 3.
xs = Vector(LinRange(xmin-ampl*dx,xmax+ampl*dx,res))
X = repeat(xs',res,1)
ys = Vector(LinRange(ymin-ampl*dy,ymax+ampl*dy,res))
Y = repeat(ys,1,res)

μ = zeros(res,res)
p = zeros(res,res)
c = zeros(res,res)
for i in 1:res
	for j in 1:res
		θ = X[i,j]*v1 + Y[i,j]*v2
		J = A.*cos.(θ*ones(1,n) - ones(n)*θ' .- ϕ)
		J -= diagm(0 => vec(sum(J,dims=2)))

		μ[i,j] = maximum(eigvals(R*(J+J')*R'/2))

		q = winding(θ,Σ)
		p[i,j] = (2.)^(q[1]+1) * (3.)^(q[2]+1)

		c[i,j] = (maximum(abs.(dcc(B'*θ))) < γbar)
	end
end
ma = maximum(μ)
mi = minimum(μ)

μtilde = zeros(res,res)
for i in 1:res
	for j in 1:res
		if μ[i,j] < 0.
			μtilde[i,j] = -μ[i,j]/mi
		else
			μtilde[i,j] = μ[i,j]/ma
		end
	end
end

figure("xxx")
#PyPlot.contourf(X,Y,μ,50,cmap=cm,vmin=-maximum(abs.(μ)),vmax=maximum(abs.(μ)))
PyPlot.contourf(X,Y,μ,150,cmap=cm,vmin=minimum(μ),vmax=-minimum(μ))
PyPlot.plot([0.,x1[1],x2[1]],[0.,x1[2],x2[2]],"ok")

colorbar(label="μ")
xlabel("x1")
ylabel("x2")

figure("yyy")
PyPlot.contourf(X,Y,μtilde,50,cmap=cm)
PyPlot.plot([0.,x1[1],x2[1]],[0.,x1[2],x2[2]],"ok")
#PyPlot.plot([0.,x1[1],x2[1],x3[1],x4[1]],[0.,x1[2],x2[2],x3[2],x4[2]],"ok")

title("μmin = $mi, μmax = $ma")
colorbar(label="μ")
xlabel("x1")
ylabel("x2")


figure("zzz")
PyPlot.contour(X,Y,log.(p),150)
#PyPlot.plot([0.,x1[1],x2[1]],[0.,x1[2],x2[2]],"ok")
PyPlot.plot(0.,0.,"ok")
PyPlot.text(0.,0.,"$(winding(θ0,Σ))")
PyPlot.plot(x1[1],x1[2],"ok")
PyPlot.text(x1[1],x1[2],"$(winding(θ1,Σ))")
PyPlot.plot(x2[1],x2[2],"ok")
PyPlot.text(x2[1],x2[2],"$(winding(θ2,Σ))")
colorbar()


figure("aaa")
PyPlot.contour(X,Y,c)
PyPlot.plot([0.,x1[1],x2[1]],[0.,x1[2],x2[2]],"ok")
colorbar()



