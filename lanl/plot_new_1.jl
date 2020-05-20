using PyPlot, DelimitedFiles

include("final_new.jl")

figure()

for i in 1:3
	Xs = readdlm("data_marc/Xs$i.csv",',')
	Ns = round.(Int,LinRange(49,50,30)*1000)

	L_l0 = Array{Float64,1}()
	L_l2 = Array{Float64,1}()
	g0 = Array{Float64,1}()
	g2 = Array{Float64,2}(undef,3,0)
	k0 = Array{Int64,1}()
	k2 = Array{Int64,1}()

	for N in Ns
		LLL = run_new(Xs[:,1:N],1e-3,50)
		push!(L_l0,LLL[2][1])
		push!(L_l2,LLL[4][1])
		push!(g0,LLL[2][4])
		g2 = [g2 LLL[4][4]]
		push!(k0,LLL[2][5])
		push!(k2,LLL[4][5])
	end

	subplot(3,3,(i-1)*3 + 1)
	PyPlot.plot(Ns*1e-3,L_l0,"o")
	PyPlot.plot(Ns*1e-3,L_l2,"-")
	xlabel("T")
	ylabel("L_min")

	subplot(3,3,(i-1)*3 + 2)
	PyPlot.plot(Ns*1e-3,g0,"o")
	PyPlot.plot(Ns*1e-3,g2[1,:],"-")
	PyPlot.plot(Ns*1e-3,g2[2,:],"--")
	PyPlot.plot(Ns*1e-3,g2[3,:],"--")
	PyPlot.plot(Ns*1e-3,ones(length(Ns)),"--k")
	xlabel("T")
	ylabel("γ")

	subplot(3,3,(i-1)*3 + 3)
	PyPlot.plot(Ns*1e-3,2*pi*(k0 .- 1)./(Ns*1e-3),"o")
	PyPlot.plot(Ns*1e-3,2*pi*(k2 .- 1)./(Ns*1e-3),"-")
	PyPlot.plot(Ns*1e-3,ones(length(Ns)),"--k")
	xlabel("T")
	ylabel("ω = 2πν")
end


