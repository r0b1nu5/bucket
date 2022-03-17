using LinearAlgebra, PyPlot

n = 10
ϕs = [.1,.2,.5]
γ = LinRange(0,π/2,200)

for i in 1:length(ϕs)
	ϕ = ϕs[i]

	subplot(1,3,1)
	PyPlot.plot(γ,2*(1-cos(π/n))*cos(ϕ)*cos.(γ),color="C$(i-1)",label="ϕ = $ϕ")
	PyPlot.plot(γ,2*sin(ϕ)*sin.(γ),"--",color="C$(i-1)")

	subplot(1,3,2)
	PyPlot.plot(γ,2*(1-cos(2π/n))*cos(ϕ)*cos.(γ),color="C$(i-1)",label="ϕ = $ϕ")
	PyPlot.plot(γ,2*sin(ϕ)*sin.(γ),"--",color="C$(i-1)")

	subplot(1,3,3)
	PyPlot.plot(γ,n*cos(ϕ)*cos.(γ),color="C$(i-1)",label="ϕ = $ϕ")
	PyPlot.plot(γ,(n-1)*sin(ϕ)*sin.(γ),"--",color="C$(i-1)")
end

subplot(1,3,1)
xlabel("γ")
ylabel("λ2, D")
legend()
title("Path graph")

subplot(1,3,2)
xlabel("γ")
legend()
title("Cycle graph")

subplot(1,3,3)
xlabel("γ")
legend()
title("Complete graph")


