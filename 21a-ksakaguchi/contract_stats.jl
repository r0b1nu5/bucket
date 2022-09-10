using PyPlot

include("ksakaguchi.jl")
include("iterations.jl")

iterate = false
read = true
plot3 = false
plot1 = true
write = false
verb = true

# #=
ntw = "ntw10"
δs = [.208,.204,.202,.2]
# =#

 #=
ntw = "rts96"
δs = [2.,1.,.5,.25]
# =#

n_iter = 10

L = readdlm("ntw_data/"*ntw*"_L.csv",',')
include("ntw_data/"*ntw*"_cycles.jl")
ω = vec(readdlm("ntw_data/"*ntw*"_om1.csv",','))
ths = ["th0","th1","th2"]
ths = ["th1",]

b,w = L2B(L)
B,Bout,Bin = L2B_bidir(L)

n,m = size(b)
m2 = 2*m

α = .1
γ = π/2 - α
hγ1 = h(-γ,α)
hγ2 = h(γ,α)

#xxx = ksakaguchi(L,ω,θ0,α,false,false,.01,1e-8)
#θf = xxx[1][:,end]

for th in ths
	θf = vec(readdlm("ntw_data/"*ntw*"_"*th*".csv",','))
	uf = winding(θf,Σ)
	
	Δf = dcc(b'*θf)
	Ff = h([Δf;-Δf])

	u = uf
	T = 1000
	
	Lmin = ones(m2)
	P = dir_cycle_proj(B,Lmin)

	if iterate
		global d1f,d2f,d12
	
		global min_decrease = Dict{Tuple{String,Float64},Vector{Float64}}()
		for δ in δs
			min_decrease[("d1f",δ)] = Vector{Float64}()
			min_decrease[("d2f",δ)] = Vector{Float64}()
			min_decrease[("d12",δ)] = Vector{Float64}()
			
#			global X1 = Vector{Matrix{Float64}}()
#			global X2 = Vector{Matrix{Float64}}()
			
			for i in 1:n_iter
				if verb
					@info "δ = $δ, iter = $i"
				end
			
				Δ01 = γ*(2*rand(m) .- 1)
				Δ1 = iterations5(Δ01,θf,Bout,B,C,ω,α,u,δ,T,false)
#				f1 = h([Δ1;-Δ1])
#				push!(X1,Δ1)
		
				Δ02 = γ*(2*rand(m) .- 1)
				Δ2 = iterations5(Δ02,θf,Bout,B,C,ω,α,u,δ,T,false)
#				f2 = h([Δ2;-Δ2])
#				push!(X2,Δ2)
	
				d1f = [norm(dcc(Δ1[:,t]-Δf)) for t in 1:T]
				d2f = [norm(dcc(Δ2[:,t]-Δf)) for t in 1:T]
				d12 = [norm(dcc(Δ1[:,t] - Δ2[:,t])) for t in 1:T]
			
				push!(min_decrease[("d1f",δ)],maximum(d1f[2:T] - d1f[1:T-1]))
				push!(min_decrease[("d2f",δ)],maximum(d2f[2:T] - d2f[1:T-1]))
				push!(min_decrease[("d12",δ)],maximum(d12[2:T] - d12[1:T-1]))
			end
		end
	end

	if read
		global min_decrease = Dict{Tuple{String,Float64},Vector{Float64}}()
		for d in ["d1f","d2f","d12"]
			for δ in δs
				min_decrease[(d,δ)] = vec(readdlm("temp_data/mindec_"*ntw*"_"*th*"_"*d*"_$(δ).csv",','))
			end
		end
	end
	
	if plot3
		figure()
		n_bins = min(n_iter,30)
	
		M = maximum([maximum(min_decrease[(v,δ)]) for v in ["d1f","d2f","d12"] for δ in δs])
		binw = M/n_bins
		bins = Vector(-binw:binw:M) .+ 1e-8
		
		for δ in δs
			subplot(1,3,1)
			PyPlot.hist(min_decrease[("d1f",δ)],bins,label="δ = $δ")
		
			subplot(1,3,2)
			PyPlot.hist(min_decrease[("d2f",δ)],bins,label="δ = $δ")
		
			subplot(1,3,3)
			PyPlot.hist(min_decrease[("d12",δ)],bins,label="δ = $δ")
		end
		
		subplot(1,3,1)
		xlabel("max. decrease of ||Δ1 - Δf||")
		ylabel("p")
		legend()
		
		subplot(1,3,2)
		xlabel("max. decrease of ||Δ2 - Δf||")
		ylabel("p")
		legend()
		
		subplot(1,3,3)
		xlabel("max. decrease of ||Δ1 - Δ2||")
		ylabel("p")
		legend()
	end
	
	if plot1
		figure()
		n_bins = min(n_iter,30)
	
		M = maximum([maximum(min_decrease[("d12",δ)]) for δ in δs])
		binw = M/n_bins
		bins = Vector(-binw:binw:M) .+ 1e-8
		
		PyPlot.hist([min_decrease[("d12",δ)] for δ in δs],bins,label=["δ = $δ" for δ in δs])
		
		xlabel("max. increase of ||Δ1 - Δ2||")
		ylabel("p")
		legend()
		title(th)
	end

	if write && !read
		for d in ["d1f","d2f","d12"]
			for δ in δs
				writedlm("temp_data/mindec_"*ntw*"_"*th*"_"*d*"_$(δ).csv",min_decrease[(d,δ)],',')
			end
		end
	end
end



