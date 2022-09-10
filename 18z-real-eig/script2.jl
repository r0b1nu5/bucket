using LinearAlgebra,PyPlot

n = 10
a = pi/6

th = (pi/2-a)*rand(n)
dth = th*ones(1,n) - ones(n,1)*transpose(th)
J = cos.(dth .- a).*(ones(n,n) - I)
J = J - diagm(0 => vec(sum(J,dims=2)))

th0 = th[1:n-1]
dth0 = th0*ones(1,n-1) - ones(n-1,1)*transpose(th0)
J0 = cos.(dth0 .- a).*(ones(n-1,n-1) - I)
J0 = J0 - diagm(0 => vec(sum(J0,dims=2)))
J0 = [J0 zeros(n-1);zeros(1,n)]

dJ = J - J0

#= 
J1 = copy(J0)
sxy = Array{Float64,2}(undef,n,0)

for i in n-1:-1:1
	global J1,dth,sxy
	x = [zeros(i-1);-cos(dth[i,n]-a);zeros(n-1-i);cos(dth[n,i]-a)]
	y = [zeros(i-1);1;zeros(n-1-i);-1]
	U = eigvecs(J1)
	Ux = inv(U)*x
	yU = vec(transpose(y)*U)
	
	sxy = [sxy sign.(Ux.*yU)]
	
	J1 = J1 + x*transpose(y)
end

matshow(sxy)
colorbar()
=#
#=
ls = Array{Float64,2}(undef,n,0)

for t in LinRange(0,1,1000)
	global ls
	J1 = J0 + t*dJ
	ls = [ls sort(eigvals(J1))]
end

for i in 1:n
	PyPlot.plot(LinRange(0,1,1000),vec(ls[i,:]))
end
=#

ls = Array{Float64,2}(undef,n,0)
ls = [ls eigvals(J0)]
J1 = J0 + diagm(0 => diag(dJ))
ls = [ls eigvals(J1)]
J1 = J1 + [zeros(n-1,n);transpose([dJ[end,1:n-1];0])]
ls = [ls eigvals(J1)]
J1 = J1 + [zeros(n,n-1) [dJ[1:n-1,end];0]]
ls = [ls eigvals(J1)]
ls = [ls eigvals(J)]



