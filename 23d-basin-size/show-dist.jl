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


