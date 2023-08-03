include("dir-kuramoto.jl")

n = 83
qmax = floor(Int64,n/4)
qs = -qmax:qmax
nq = length(qs)
Qs = zeros(nq,nq)

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
	θf = dir_kuramoto(L,θi,zeros(n))
	qf = winding(θf,σ)

	Qs[qi+qmax+1,qf+qmax+1] += 1
end

# #=
matshow(Qs)
xticks(collect(0:nq-1),collect(-qmax:qmax))
yticks(collect(0:nq-1),collect(-qmax:qmax))
# =#

figure("Histogram")
for i in 1:nq
	qi = i-1-qmax
	s = sum(Qs[i,:])
	μ = sum(Qs[i,:].*(-qmax:qmax))/s
	if s > 10
		subplot(1,2,1)
		PyPlot.plot(-qmax:qmax,Qs[i,:]./s,"-o")
		subplot(1,2,2)
		PyPlot.plot(-qmax-μ:qmax-μ,Qs[i,:]./s,"-o")
	end
end

subplot(1,2,1)
xlabel("qf")
ylabel("p")
subplot(1,2,2)
xlabel("qf-μ")
ylabel("p")



