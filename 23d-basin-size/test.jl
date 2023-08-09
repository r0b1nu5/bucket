using Statistics

include("dir-kuramoto.jl")
include("cycle-kuramoto.jl")
include("tools.jl")

 #= 
n = 83
qmax = floor(Int64,n/4)
qs = -qmax:qmax
nq = length(qs)
Qs = zeros(nq,nq)
Qi = Dict{Int64,Vector{Int64}}(q => Int64[] for q in -qmax:qmax)
Qf = Dict{Int64,Vector{Int64}}(q => Int64[] for q in -qmax:qmax)

#L = gen_cycle_dir(n)
L = gen_cycle_undir(n)
σ = collect(1:n)

n_iter = 1000

for iter in 1:n_iter
	if iter%100 == 0
		@info "############### iter: $iter/$n_iter"
	end

	θi = 2π*rand(n)
	qi = winding(θi,σ)
	θf = cycle_kuramoto(θi,zeros(n),ones(n))
	qf = winding(θf,σ)

	Qs[qi+qmax+1,qf+qmax+1] += 1
	push!(Qi[qf],qi)
	push!(Qf[qi],qf)
end

# =#

QQf = collect(values(Qf))[.!isempty.(values(Qf))]
qfmax = maximum([maximum(QQf[i]) for i in 1:length(QQf)])
qfmin = minimum([minimum(QQf[i]) for i in 1:length(QQf)])
QQi = collect(values(Qi))[.!isempty.(values(Qi))]
qimax = maximum([maximum(QQi[i]) for i in 1:length(QQi)])
qimin = minimum([minimum(QQi[i]) for i in 1:length(QQi)])
qifmax = maximum([qfmax,-qfmin,qimax,-qimin])

# #=
matshow(Qs)
xticks(collect(0:nq-1),collect(-qmax:qmax))
yticks(collect(0:nq-1),collect(-qmax:qmax))
# =#

Mi = Float64[]
Mf = Float64[]
Σi = Float64[]
Σf = Float64[]
Pni = Float64[]
Pnf = Float64[]
Pei = Float64[]
Pef = Float64[]

min_stat = 10

figure("Histograms")
for i in 1:nq
	qi = i-1-qmax
	si = sum(Qs[i,:])
#	μi = sum(Qs[i,:].*(-qmax:qmax))/si
	μi = mean(Qf[qi])
	push!(Mi,μi)
	σi = std(Qf[qi])
	push!(Σi,σi)
	if si >= min_stat
		subplot(2,2,1)
		PyPlot.plot(-qmax:qmax,Qs[i,:]./si,"-o")
		subplot(2,2,2)
		PyPlot.plot(-qmax-μi:qmax-μi,Qs[i,:]./si,"-o")
	
		push!(Pni,gof_normal(Qf[qi],qmax,50))
		push!(Pei,gof_exp(Qf[qi],qmax,50))
	else
		push!(Pni,NaN)
		push!(Pei,NaN)
	end

	qf = i-1-qmax
	sf = sum(Qs[:,i])
#	μf = sum(Qs[:,i].*(-qmax:qmax))/sf
	μf = mean(Qi[qf])
	push!(Mf,μf)
	σf = std(Qi[qf])
	push!(Σf,σf)
	if sf > min_stat
		subplot(2,2,3)
		PyPlot.plot(-qmax:qmax,Qs[:,i]./sf,"-o")
		subplot(2,2,4)
		PyPlot.plot(-qmax-μf:qmax-μf,Qs[:,i]./sf,"-o")
	
		push!(Pnf,gof_normal(Qi[qf],qmax,50))
		push!(Pef,gof_exp(Qi[qf],qmax,50))
	else
		push!(Pnf,NaN)
		push!(Pef,NaN)
	end
end

subplot(2,2,1)
xlabel("qf")
ylabel("P(qf|qi)")
subplot(2,2,2)
PyPlot.plot(-qmax:qmax,exp.(-(((-qmax:qmax) .- Mi[qmax+1]).^2)./(2Σi[qmax+1]^2))./(sqrt(2π)*Σi[qmax+1]),"--k")
xlabel("qf-μ")
ylabel("P(qf|qi)")

subplot(2,2,3)
xlabel("qi")
ylabel("P(qi|qf)")
subplot(2,2,4)
PyPlot.plot(-qmax:qmax,exp.(-(((-qmax:qmax) .- Mf[qmax+1]).^2)./(2Σf[qmax+1]^2))./(sqrt(2π)*Σf[qmax+1]),"--k")
xlabel("qi-μ")
ylabel("P(qi|qf)")

figure("Gaussian parameters")
PyPlot.plot(-qmax:qmax,Mi,color="C0")
PyPlot.plot(-qmax:qmax,Σi,"--",color="C0")
PyPlot.plot(-qmax:qmax,Mf,color="C1")
PyPlot.plot(-qmax:qmax,Σf,"--",color="C1")
xlabel("qi or qf")

figure("GoFs Normal")
PyPlot.bar(-qmax:qmax,Pni,color="C0",align="edge",width=-.4)
PyPlot.bar(-qmax:qmax,Pnf,color="C1",align="edge",width=.4)
PyPlot.plot([-qmax,qmax],[.05,.05],"--k")
axis([-qifmax-1,qifmax+1,0.,1.])
xlabel("qi or qf")
ylabel("p-value")
xticks(-qifmax:qifmax)

figure("GoFs Exponential")
PyPlot.bar(-qmax:qmax,Pei,color="C0",align="edge",width=-.4)
PyPlot.bar(-qmax:qmax,Pef,color="C1",align="edge",width=.4)
PyPlot.plot([-qmax,qmax],[.05,.05],"--k")
axis([-qifmax-1,qifmax+1,0.,1.])
xlabel("qi or qf")
ylabel("p-value")
xticks(-qifmax:qifmax)





