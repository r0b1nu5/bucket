using PyPlot, LinearAlgebra

nm = 4
nn = 2
ns = 3

A = 2*rand(nm+nn,nm+nn) .- 1

Axx = A[1:nm,1:nm]
Axy = A[1:nm,nm+1:nm+nn]
Ayx = A[nm+1:nm+nn,1:nm]
Ayy = A[nm+1:nm+nn,nm+1:nm+nn]

B = Symmetric(2*rand(ns,ns).-1)
us = eigvecs(B)

T = us[:,1:nn]
TT = [diagm(0 => ones(nm)) zeros(nm,nn);zeros(ns,nm) T]

#At = [Axx Axy*T';T*Ayx T*Ayy*T']
At = TT*A*TT'

z0 = 2*rand(nm+nn) .- 1
x0 = z0[1:nm]
y0 = z0[nm+1:nm+nn]
#zt0 = [x0;T*y0]
zt0 = TT*z0

ts = LinRange(0,1,200)

z = zeros(nm+nn,0)
zt = zeros(nm+ns,0)

for t in ts
	global z = [z exp(A*t)*z0]
	global zt = [zt exp(At*t)*zt0]
end

figure()

for i in 1:nm
	subplot(2,2,1)
	PyPlot.plot(ts,z[i,:])
	subplot(2,2,3)
	PyPlot.plot(ts,zt[i,:])
end

for i in nm+1:nm+nn
	subplot(2,2,2)
	PyPlot.plot(ts,z[i,:])
end

for i in nm+1:nm+ns
	subplot(2,2,4)
	PyPlot.plot(ts,zt[i,:])
end





