include("iterations.jl")

n_iter = 1000

n,m = size(Bout)

ρs = Array{Float64,1}()
f1s = Array{Float64,2}(undef,m,0)
f2s = Array{Float64,2}(undef,m,0)

iter = 0
while length(ρs) < 1
	global iter += 1
#	if iter%1000 == 0
		@info "iter = $iter"
#	end

	f01 = rand_init(ω,Bout,hγ)
	f02 = rand_init(ω,Bout,hγ)
	df0 = f01-f02
	
	f1 = Tu_dir(f01,u,P,C,diagm(0 => Lmin))
	f2 = Tu_dir(f02,u,P,C,diagm(0 => Lmin))
	df = f1-f2

	ρ = norm(df)/norm(df0)

	if ρ > 1.
		push!(ρs,ρ)
		f1s = [f1s f01]
		f2s = [f2s f02]
	end
end







