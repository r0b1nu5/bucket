using Statistics, DelimitedFiles, PyPlot

include("tools.jl")
include("volumes.jl")

n = 83
qmax = floor(Int64,n/4)

#ns = [43,83,123,163]

type = "undir"
n_intra = 1000

ids = Int64.(vec(readdlm("data/ids-"*type*"-$n-$(n_intra).csv",',')))
Qs = load_Qs("data/",ids,type,n,n_intra)

cmap = get_cmap("plasma")

qil = [length(Qs[qi]) for qi in keys(Qs)]
qis = [qi for qi in keys(Qs)][qil .> 10]

figure("Distribtions")

subplot(1,2,1)
sf = .16
PyPlot.plot(-qmax:qmax,exp.((-(-qmax:qmax).^2)./(2*n*sf^2)),"--k")
ccol = -1
for qi in sort(qis)
	global ccol += 1
	
	qf = collect(keys(Qs[qi]))
	nqf = [Qs[qi][q] for q in qf]
	qfnqf = sortslices([qf nqf],dims=1)
	qf = qfnqf[:,1]
	nqf = qfnqf[:,2]
	
	μqf = sum(qf.*nqf)/sum(nqf)
	ηqf = nqf./maximum(nqf)
	
	PyPlot.plot(qf .- μqf,ηqf,color=cmap(ccol/(length(qis)-1)),label="qi = $qi",alpha=.4)
end
xlabel("qf - μqf")
ylabel("P(qf|qi)")
legend()
title("n = $n")

subplot(1,2,2)
ccol = -1
si = .25
PyPlot.plot(-qmax:qmax,exp.((-(-qmax:qmax).^2)./(2*n*si^2)),"--k")
for qf in -qmax:qmax
	global ccol += 1

	qi = [haskey(Qs[q],qf) ? q : NaN for q in sort(collect(keys(Qs)))]
	nqi = [haskey(Qs[q],qf) ? Qs[q][qf] : 0 for q in sort(collect(keys(Qs)))]
	qi = qi[nqi .> 10]
	nqi = nqi[nqi .> 10]

	if length(nqi) > 0
		μqi = sum(qi.*nqi)/sum(nqi)
		ηqi = nqi./maximum(nqi)
		
		PyPlot.plot(qi .- μqi,ηqi,color=cmap(ccol/(2*qmax)),label="qf = $qf",alpha=.4)
	end
end
xlabel("qi - μqi")
ylabel("P(qi|qf)")
legend()
title("n = $n")


