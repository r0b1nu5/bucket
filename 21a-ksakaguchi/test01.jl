using PyPlot, Statistics

include("ksakaguchi.jl")
include("tools.jl")

n = 12

th0 = Array(2pi/n*(1:n))
th1 = copy(th0)
th2 = copy(th0)

L = cyqle(n)

a = .5

om = [-1;1;zeros(n-2)]

bs = LinRange(.55,.8,100)

X = Array{Float64,2}(undef,n,0)
Y = Array{Float64,2}(undef,n,0)
ls1 = Array{Float64,2}(undef,n,0)
ls2 = Array{Float64,2}(undef,n,0)

for b in bs
	global th1,th2,X,Y,ls1,ls2

	x = ksakaguchi(L,b*om,th1,a,false)
	y = ksakaguchi(L,b*om,th2,0.,false)
	
	th1 = x[1][:,end]
	X = [X th1]
	th2 = y[1][:,end]
	Y = [Y th2]

	ls1 = [ls1 eigvals(jacobian(L,th1,a))]
	ls2 = [ls2 eigvals(jacobian(L,th2,0.))]
end

nn,T1 = size(X)
dX = mod.(X - [X[2:end,:];X[[1,],:]] .+ pi,2pi) .- pi
q1 = [winding(X[:,i],Array(1:n)) for i in 1:T1]
nn,T2 = size(Y)
dY = mod.(Y - [Y[2:end,:];Y[[1,],:]] .+ pi,2pi) .- pi
q2 = [winding(Y[:,i],Array(1:n)) for i in 1:T2]

figure()
subplot(2,2,1)
PyPlot.plot([bs[1],bs[end]],[pi/2+a,pi/2+a],"--b")
PyPlot.plot([bs[1],bs[end]],[pi/2-a,pi/2-a],"--r")
PyPlot.plot([bs[1],bs[end]],[-pi/2+a,-pi/2+a],"--b")
PyPlot.plot([bs[1],bs[end]],[-pi/2-a,-pi/2-a],"--r")
PyPlot.plot([bs[1],bs[end]],[pi/2,pi/2],"--k")
PyPlot.plot([bs[1],bs[end]],[-pi/2,-pi/2],"--k")
subplot(2,2,2)
PyPlot.plot([bs[1],bs[end]],[pi/2,pi/2],"--k")
PyPlot.plot([bs[1],bs[end]],[-pi/2,-pi/2],"--k")
for i in 1:n
	subplot(2,2,1)
	PyPlot.plot(bs,dX[i,:])
	subplot(2,2,2)
	PyPlot.plot(bs,dY[i,:])
	subplot(2,2,3)
	PyPlot.plot(bs,ls1[i,:])
	subplot(2,2,4)
	PyPlot.plot(bs,ls2[i,:])
end
#=
subplot(2,2,3)
PyPlot.plot(bs,q1)
subplot(2,2,4)
PyPlot.plot(bs,q2)
=#





