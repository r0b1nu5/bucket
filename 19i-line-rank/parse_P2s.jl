using PyPlot,DelimitedFiles

P0 = .11
taus = [.025,.075,.125,.175,.225,.275,.325,.375,.425,.475]

P2s = Array{Float64,2}(undef,165,0)
losses = Array{Float64,2}(undef,165,0)

for tau in taus
	global P2s,losses
	
	P2s = [P2s Array{Float64,1}(vec(readdlm("data/P2s_$(P0)_$(tau).csv",',')))]
	losses = [losses Array{Float64,1}(vec(readdlm("data/losses_$(P0)_$(tau).csv",',')))]
end

PL = Array{Float64,2}(undef,0,10)
for i in 1:165
	global PL
	if !isinf(P2s[i,1]/losses[i,1])
		PL = [PL;P2s[[i,],:]./losses[[i,],:]]
	end
end
miPL = minimum(PL,dims=1)
maPL = maximum(PL,dims=1)
l = size(PL)[1]
PLn = (PL - repeat(miPL,l))./(repeat(maPL,l) - repeat(miPL,l))
PyPlot.subplot(1,3,1)
for i in 1:l
		PyPlot.plot(taus,vec(PLn[i,:]),"-o")
end
xlabel("τ")
ylabel("P2/losses")
for i in 1:165
	if !isinf(P2s[i,1]./losses[i,1])
		PyPlot.subplot(1,3,2)
		PyPlot.plot(taus,vec(P2s[i,:]./losses[i,:]),"-o")
		PyPlot.subplot(1,3,3)
		PyPlot.loglog(taus,vec(P2s[i,:]./losses[i,:]),"-o")
	end
end
PyPlot.subplot(1,3,2)
xlabel("τ")
ylabel("P2/losses")
PyPlot.subplot(1,3,3)
xlabel("τ")
ylabel("P2/losses")



