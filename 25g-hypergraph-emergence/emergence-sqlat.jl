using Random, Dates, DelimitedFiles

include("kuramoto.jl")
include("this.jl")
include("this-filter-distance.jl")
include("graph-tools.jl")
include("coarse-grain.jl")
include("gen-lattice.jl")

@info "############# START: $(now())"

n = 8
T = 500
ks = [1,2]

save = true # Saving the inference and the ground truth?

A,B,coord = gen_square_lattice(n)
m = Int64(nnz(A)/2)

# ========================================================================

amplitude = 1.
ξ0 = 0.0005
X = amplitude*(rand(n^2,T) .- .5)
Y = f_kuramoto(X,zeros(n^2),B,1.,π/4)

 #= # Normalizing the states
for i in 1:size(X)[1]
	X[i,:] = X[i,:]./mean(X[i,:])
end
# =#

#################### INFERENCE ###################################

ooi = [2,3,4]
dmax = 3
zer0 = 1e-4
λ = 0.1

#Ainf,coeff,relerr = this(X,Y,ooi,dmax,λ)

dist_keep = 2.1
Ainf,coeff,relerr = this_filter_distance(X,Y,ooi,dmax,dist_keep,coord,λ)

k2Ainf = Dict{Int64,Any}(0 => Ainf)

m2max = 1.
m3max = 1.
m4max = 1.
#m2max = (n*(n-1))
#m3max = (n*(n-1)*(n-2))
m2 = [sum(abs.(Ainf[2][:,3]) .> zer0)/m2max,]
m3 = [sum(abs.(Ainf[3][:,4]) .> zer0)/m3max,]
m3 = [sum(abs.(Ainf[4][:,5]) .> zer0)/m4max,]
m2true = [2*(n-1)^2+2*(n-1),]

list_3edges = Vector{Vector{Tuple{Int64,Int64,Int64}}}()
push!(list_3edges,Tuple{Int64,Int64,Int64}[])
for i in 1:size(Ainf[3])[1]
	if abs(Ainf[3][i,4]) > zer0
		push!(list_3edges[1],(Ainf[3][i,1],Ainf[3][i,2],Ainf[3][i,3]))
	end
end

# Contribution of 2-edges to the dynamics
global mag2 = zeros(size(X)[2])
for i in 1:size(Ainf[2])[1]
	a = Int64(Ainf[2][i,2])
	global mag2 += abs.(Ainf[2][i,3]*X[a,:])
end
contribution2 = [median(mag2),]
# Contribution of 3-edges to the dynamics
global mag3 = zeros(size(X)[2])
for i in 1:size(Ainf[3])[1]
	a = Int64(Ainf[3][i,2])
	b = Int64(Ainf[3][i,3])
	global mag3 += abs.(Ainf[3][i,4]*(X[a,:].*X[b,:]))
end
contribution3 = [median(mag3),]
global mag4 = zeros(size(X)[2])
for i in 1:size(Ainf[4])[1]
	a = Int64(Ainf[4][i,2])
	b = Int64(Ainf[4][i,3])
	c = Int64(Ainf[4][i,4])
	global mag4 += abs.(Ainf[4][i,5]*(X[a,:].*X[b,:].*X[c,:]))
end
contribution4 = [median(mag4),]

distances = Vector{Int64}[]

for k in ks
	@info "k/kmax = $k/$(maximum(ks))"
	
	global A2,B2,X = coarse_grain_sqlat_4(n,X)
	global A2,B2,Y = coarse_grain_sqlat_4(n,Y)
	X2 = X
	Y2 = Y
	global n = Int64(n/2)
	push!(m2true,2*(n-1)^2 + 2*(n-1))

	A,B,coord = gen_square_lattice(n)
#	Ainf2,coeff2,relerr2 = this(X2,Y2,ooi,dmax,λ)
	Ainf2,coeff2,relerr2 = this_filter_distance(X2,Y2,ooi,dmax,dist_keep,coord,λ)

	global k2Ainf[k] = Ainf2

	push!(m2,sum(abs.(Ainf2[2][:,3]) .> zer0))
	push!(m3,sum(abs.(Ainf2[3][:,4]) .> zer0))
	#push!(m2,sum(abs.(Ainf2[2][:,3]) .> zer0)/(nc*(nc-1)))
	#push!(m3,sum(abs.(Ainf2[3][:,4]) .> zer0)/(nc*(nc-1)*(nc-2)))
	
	push!(list_3edges,Tuple{Int64,Int64,Int64}[])
	for i in 1:size(Ainf2[3])[1]
		if abs(Ainf2[3][i,4]) > zer0
			push!(list_3edges[end],(Ainf2[3][i,1],Ainf2[3][i,2],Ainf2[3][i,3]))
		end
	end

	# Contribution of 2-edges to the dynamics
	mag2 = zeros(size(X2)[2])
	for i in 1:size(Ainf2[2])[1]
		a = Int64(Ainf2[2][i,2])
		mag2 += abs.(Ainf2[2][i,3]*X2[a,:])
	end
	push!(contribution2,median(mag2))
	# Contribution of 3-edges to the dynamics
	mag3 = zeros(size(X2)[2])
	for i in 1:size(Ainf2[3])[1]
		a = Int64(Ainf2[3][i,2])
		b = Int64(Ainf2[3][i,3])
		mag3 += abs.(Ainf2[3][i,4]*(X2[a,:].*X2[b,:]))
	end
	push!(contribution3,median(mag3))
	mag4 = zeros(size(X)[2])
	for i in 1:size(Ainf2[4])[1]
		a = Int64(Ainf2[4][i,2])
		b = Int64(Ainf2[4][i,3])
		c = Int64(Ainf2[4][i,4])
		mag4 += abs.(Ainf[4][i,5]*(X[a,:].*X[b,:].*X[c,:]))
	end
	push!(contribution4,median(mag4))
end


# Plot number of edges (normalized)
fig1, (ax11,ax21) = subplots(1,2,figsize=(15,5))

ax11.plot([0;ks], m2, color="C0")
ax11.plot([0;ks], 2*m2true, "ok")
ax11.set_xlabel("k")
ax11.set_ylabel("#2-edges", color="C0")
ax11.set_ylim(-maximum(m2)*0.05,maximum(m2)*1.05)
ax11.tick_params(axis="y", labelcolor="C0")

ax12 = ax11.twinx()   # Share the same X-axis
ax12.plot([0;ks], m3, color="C1")
ax12.set_ylabel("#3-edges", color="C1")
ax12.set_ylim(-maximum(m3)*0.05,maximum(m3)*1.05)
ax12.tick_params(axis="y", labelcolor="C1")

title("λ = $λ")

tight_layout()
show()

# Plot contribution of edges
#fig2, ax21 = subplots(1,2,2)

ax21.plot([0;ks], contribution2, color="C0")
ax21.set_xlabel("k")
ax21.set_ylabel("2-edges contribution", color="C0")
ax21.set_ylim(-maximum(contribution2)*0.05,maximum(contribution2)*1.05)
ax21.tick_params(axis="y", labelcolor="C0")

ax22 = ax21.twinx()
ax22.plot([0;ks], contribution3, color="C1")
ax22.set_ylabel("3-edges contribution", color="C1")
ax22.set_ylim(-maximum(contribution3)*0.05,maximum(contribution3)*1.05)
ax22.tick_params(axis="y", labelcolor="C1")

tight_layout()
show()







