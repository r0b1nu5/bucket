using PyPlot

function plot_osc(Vs::Array{Float64,2}, θs::Array{Float64,1}, vs::Array{Float64,1}, nc::Int64=3)
	nn,T = size(Vs)
	n = Int64(nn/2)

	nl = ceil(Int64,n/nc)

	ϕ = atan(Vs[2,T]/Vs[1,T])
	θ = ϕ .- θs

	figure()
	for i in 1:n
		subplot(nl,nc,i)
		PyPlot.plot(vs[i]*cos.(LinRange(0,2π,200)),vs[i]*sin.(LinRange(0,2π,200)),"--",color="C1")
		PyPlot.plot(vs[i]*cos(θ[i]),vs[i]*sin(θ[i]),"o",color="C1")

		PyPlot.plot(Vs[2*i-1,:],Vs[2*i,:],color="C0")
		PyPlot.plot(Vs[2*i-1,1],Vs[2*i,1],"o",color="C3")
		PyPlot.plot(Vs[2*i-1,T],Vs[2*i,T],"o",color="C2")


	end

	return nothing
end



