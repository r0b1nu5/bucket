using LinearAlgebra,PyPlot

n = 10
a = pi/5
#th = pi/2*rand(n)
dth = th*ones(1,n) - ones(n,1)*transpose(th)
D = Array{Float64,1}()
for i in 1:n
	push!(D,sum(cos.(th.-th[i])))
end
Dv = D*ones(1,n)
Dh = ones(n,1)*transpose(D)

al1 = .5*sum(sin.(dth).^2)
be1 = -.5*sum(sin.(dth).^2 .*(Dv+Dh))
ga1 = .5*sum(sin(dth).^2 .*Dv.*Dh)
de1 = -n*cos(a)
ep1 = cos(a)*sum(D)
ph1 = 1

al2 = -(n-2)/2*sum(sin.(dth).^2)
be2 = 0.
ga2 = 0.
de2 = 0.
for i in 1:n
	for j in 1:n
		for k in 1:n
			if i!=j && i!=k && j!=k
				global be2,ga2,de2
				be2 += .5*sin(th[i]-th[j]).^2*(D[i]+D[j]+D[k])
				ga2 += -.5*sin(th[i]-th[j]).^2*(D[i]*D[j]+D[j]*D[k]+D[k]*D[i])
				de2 += .5*sin(th[i]-th[j]).^2*D[i]*D[j]*D[k]
			end
		end
	end
end
ep2 = (n^2-n)*cos(a)
ph2 = 0.
mu2 = 0.
for i in 1:n
	for j in 1:n
		if i!=j
			global ph2,mu2
			ph2 += -cos(a)*(D[i]+D[j])
			mu2 += cos(a)*D[i]*D[j]
		end
	end
end
nu2 = -n
et2 = sum(D)

ls = LinRange(-300,300,1000)
rs = LinRange(-5,5,500)

D1 = zeros(1000,500)
D2 = zeros(1000,500)

for i in 1:1000
	for j in 1:500
		D1[i,j] = (al1*ls[i]^2+be1*ls[i]+ga1)*rs[j]^2+(de1*ls[i]+ep1)*rs[j]+ph1
		D2[i,j] = (al2*ls[i]^3*be2*ls[i]^2+ga2*ls[i]+de2)*rs[j]^2+(ep2*ls[i]^2+ph2*ls[i]+mu2)*rs[j]+nu2*ls[j]+et2
	end
end

a1 = al1*ls.^2 + be1*ls .+ ga1
b1 = de1*ls .+ ep1
c1 = ones(1000)

matshow(D1)
matshow(sign.(D1))
PyPlot.plot([250,250],[0,1000],"--r")
PyPlot.plot([375,375],[0,1000],"--r")
matshow(D2)
matshow(sign.(D2))
PyPlot.plot([250,250],[0,1000],"--r")
PyPlot.plot([375,375],[0,1000],"--r")
matshow(sign.(D1.*D2))
PyPlot.plot([250,250],[0,1000],"--r")
PyPlot.plot([375,375],[0,1000],"--r")

figure()
surf(D1)

figure()
PyPlot.plot(l,a1,label="a1")
PyPlot.plot(l,b1,label="b1")
PyPlot.plot(l,c1,label="c1")
legend()



