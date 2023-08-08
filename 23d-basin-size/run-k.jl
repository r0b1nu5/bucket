using SparseArrays

include("dir-kuramoto.jl")
include("tools.jl")

n_iter = 36000
n_intra = 1000
n_loop = Int64(n_iter/n_intra)

n_perc = ceil(Int64,n_iter/100)

n = 83

L = gen_cycle_undir(n); type = "undir"
#L = get_cycle_dir(n); type = "dir"
σ = collect(1:n)

ω = zeros(n)

qmax = floor(Int64,n/4)
nq = floor(Int64,n/2+1)

iter = 0
ids = Int64[]
for str in readdir("data/")
	if str[5:end] == "-k-"*type*"-$n-$(n_intra)-Qs.csv"
		push!(ids,parse(Int64,str[1:4]))
	end
end

for l in 1:n_loop
	id = rand(1111:9999)
	while id in ids
		id = rand(1111:9999)
	end

	Qs = init_Qs(n)
	for i in 1:n_intra
		global iter += 1

		θi = 2π*rand(n)
		qi = winding(θi,σ)
		θf = dir_kuramoto(L,θi,ω)
		qf = winding(θf,σ)

		Qs[qi][qf] += 1

		if iter%n_perc == 0
			@info "Achieved: $(iter/n_iter*100)%"
		end
	end

	save_Qs(Qs,"data/$id-k-"*type*"-$n-$(n_intra)-")
	push!(ids,id)
	io = open("data/ids-"*type*"-$n-$(n_intra).csv","a")
	write(io,"$id\n")
	close(io)
end





