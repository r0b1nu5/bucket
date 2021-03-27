using PyPlot

include("ksakaguchi.jl")
include("tools.jl")

n = 10

L = cyqle(n)
ω = [.12,.0,.675,.195,-.21,-.66,.75,-.015,-.015,-.84]
θ0 = zeros(n)
θ1 = -2π*Array(1:n)/n
α = 1.

T00 = 15.
T01 = 100.
T10 = 15.
T11 = 15.

t00,x00 = ksakaguchi_ND(L,ω,θ0,0.,(0.,T00))
t01,x01 = ksakaguchi_ND(L,ω,θ0,1.,(0.,T01))
t10,x10 = ksakaguchi_ND(L,ω,θ1,0.,(0.,T10))
t11,x11 = ksakaguchi_ND(L,ω,θ1,1.,(0.,T11))

for i in 1:n
	subplot(2,2,1)
	PyPlot.plot(t00,x00[i,:])
	subplot(2,2,2)
	PyPlot.plot(t01,x01[i,:])
	subplot(2,2,3)
	PyPlot.plot(t10,x10[i,:])
	subplot(2,2,4)
	PyPlot.plot(t11,x11[i,:])
end

subplot(2,2,1)
title("w = 0, α = 0.0")
ylabel("θ")
subplot(2,2,2)
title("w = 0, α = 1.0")
subplot(2,2,3)
title("w = 1, α = 0.0")
xlabel("t")
ylabel("θ")
subplot(2,2,4)
title("w = 1, α = 1.0")
xlabel("t")

