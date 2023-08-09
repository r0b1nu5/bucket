using SparseArrays, Dates

include("dir-kuramoto.jl")
include("cycle-kuramoto.jl")
include("tools.jl")

n_iter = 100000
n_intra = 1000
n_loop = Int64(n_iter/n_intra)

n_perc = ceil(Int64,n_iter/100)

n = 23

L = gen_cycle_undir(n); type = "undir"
#L = get_cycle_dir(n); type = "dir"
σ = collect(1:n)

ω = zeros(n)

qmax = floor(Int64,n/2)
nq = 2*qmax+1

iter = 0
ids = Int64[]
for str in readdir("data/")
	if str[5:end] == "-k-"*type*"-$n-$(n_intra)-Qs.csv"
		push!(ids,parse(Int64,str[1:4]))
	end
end

t0 = now()
for l in 1:n_loop
	id = rand(1111:9999)
	while id in ids
		id = rand(1111:9999)
	end

	global Qs = init_Qs(n)
	for i in 1:n_intra
		global iter += 1

		θi = 2π*rand(n)
		qi = winding(θi,σ)
#		θf = dir_kuramoto(L,θi,ω)
		θf = cycle_kuramoto(θi,ω)
		qf = winding(θf,σ)

		Qs[qi][qf] += 1

		if iter%n_perc == 0
			dt = Dates.value(now() - t0)
			prop = iter/n_iter
			dtfin = canonicalize(Millisecond(round(Int64,dt/prop*(1-prop))))
			
			@info "Achieved: $(iter/n_iter*100)%, remaining time: $dtfin"
		end
	end

	save_Qs(Qs,"data/$id-k-"*type*"-$n-$(n_intra)-")
	push!(ids,id)
	io = open("data/ids-"*type*"-$n-$(n_intra).csv","a")
	write(io,"$id\n")
	close(io)
end





