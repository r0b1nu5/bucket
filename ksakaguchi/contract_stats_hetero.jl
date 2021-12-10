using PyPlot

include("ksakaguchi.jl")
include("iterations.jl")

iterate = true
read = false
plot3 = false
plot1 = true
write = false
verb = true

 #=
ntw = "ntw10"
δs = [.2,.15,.1,.05]
α = .2*rand(m2) .+ .1
w = rand(m2) .+ .5
# =#

# #=
ntw = "rts96_w"
δs = [.5,.25,.1,.05]
Lb = readdlm("ntw_data/"*ntw*"_Lb.csv",',')
Bb,wb = L2B(Lb)
Lg = readdlm("ntw_data/"*ntw*"_Lg.csv",',')
Bg,wg = L2B(Lg)
α = atan.(wg./wb)
α = [α;copy(α).+1e-8]
w = sqrt.(wb.^2 + wg.^2)
w = [w;copy(w).+1e-8]
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

hh = Vector{Function}()
for i in 1:m
	push!(hh,(x -> w[i]*(sin(x - α[i]) + sin(α[i]))))
end
for i in 1:m
	push!(hh,(x -> w[i]*(sin(x - α[i]) + sin(α[i]))))
end
#hh = [(x -> w[i]*(sin(x - α[i]) + sin(α[i]))) for i in 1:m2]
γs = [(α[i] - π/2,π/2 - α[i]) for i in 1:m2]
@info "$(γs[109])"
@info "$(hh[109](0.))"
hγs = [(hh[i](γs[i][1]),hh[i](γs[i][2])) for i in 1:m2]
λs = ([max(γs[i][1],-γs[i+m][2]) for i in 1:m],[min(γs[i][2],-γs[i+m][1]) for i in 1:m])

#xxx = ksakaguchi(L,ω,θ0,α,false,false,.01,1e-8)
#θf = xxx[1][:,end]

for th in ths
	θf = vec(readdlm("ntw_data/"*ntw*"_"*th*".csv",','))
	uf = winding(θf,Σ)
	
	Δf = dcc(b'*θf)
	Δff = [Δf;-Δf]
	Ff = [hh[i](Δff[i]) for i in 1:m2]

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
			
				Δ01 = (λs[2] - λs[1]).*rand(m) + λs[1]
				Δ1 = iterations5_hetero(Δ01,θf,Bout,B,C,ω,hh,γs,u,δ,T,false)
#				f1 = h([Δ1;-Δ1])
#				push!(X1,Δ1)

				Δ02 = (λs[2] - λs[1]).*rand(m) + λs[1]
				Δ2 = iterations5_hetero(Δ02,θf,Bout,B,C,ω,hh,γs,u,δ,T,false)
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



