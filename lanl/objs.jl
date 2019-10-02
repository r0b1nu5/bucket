include("objective.jl")

fs = LinRange(0.,.05,201)

o1 = Array{Float64,1}()
o2 = Array{Float64,1}()
o3 = Array{Float64,1}()
o4 = Array{Float64,1}()
o5 = Array{Float64,1}()

T1 = Int(1e5)
T2 = Int(5e4)
T3 = Int(1e4)
T4 = Int(5e3)
T5 = Int(1e3)

for f in fs
	@info "$(round(f/maximum(fs)*100))%"

	push!(o1,objective1(Xs, diagm(0 => 1 ./m)*L, d./m, [.2;zeros(9)], [f;zeros(9)], zeros(10), dt, 0., 0.))
	oo = 0.
	for i in 0:1
		oo += objective1(Xs[:,i*T2+1:(i+1)*T2], diagm(0 => 1 ./m)*L, d./m, [.2;zeros(9)], [f;zeros(9)], zeros(10), dt, 0., 0.)
	end
	push!(o2,oo)

	ooo = 0.
	for i in 0:9
		ooo += objective1(Xs[:,i*T3+1:(i+1)*T3], diagm(0 => 1 ./m)*L, d./m, [.2;zeros(9)], [f;zeros(9)], zeros(10), dt, 0., 0.)
	end
	push!(o3,ooo)

	oooo = 0.
	for i in 0:19
		oooo += objective1(Xs[:,i*T4+1:(i+1)*T4], diagm(0 => 1 ./m)*L, d./m, [.2;zeros(9)], [f;zeros(9)], zeros(10), dt, 0., 0.)
	end
	push!(o4,oooo)

	ooooo = 0.
	for i in 0:99
		ooooo += objective1(Xs[:,i*T5+1:(i+1)*T5], diagm(0 => 1 ./m)*L, d./m, [.2;zeros(9)], [f;zeros(9)], zeros(10), dt, 0., 0.)
	end
	push!(o5,ooooo)
end

PyPlot.plot(fs,o5,label="T=1e3")
PyPlot.plot(fs,o4,label="T=5e3")
PyPlot.plot(fs,o3,label="T=1e4")
PyPlot.plot(fs,o2,label="T=5e4")
PyPlot.plot(fs,o1,label="T=1e5")

legend()





