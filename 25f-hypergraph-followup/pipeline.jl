using PyPlot

include("gen-hypergraph.jl")
include("hyper-kuramoto.jl")
include("tools.jl")

include("../../THIS/this.jl")

n = 10
p = .1
ω0 = .0
ξ0 = 5.
ρ = 5.
λ = .005

h = .01
nstep = 5
δt = nstep*h
niter1 = 10000
#niter2 = 5000; niter3 = 5000
#niter2 = 2000; niter3 = 8000
niter2 = 1000; niter3 = 2000
niter2 = 500; niter2bis = 50; niter3 = 1450

cm1 = get_cmap("Greens")
cm2 = get_cmap("Purples")
cm3 = get_cmap("Oranges")

A2 = zeros(0,3)
A3,B,B2,E = rand_3_digraph(n,p)
l_ref = [A3[i,1:3] for i in 1:length(A3[:,1])]
idx_ref = get_loose_ends(A2,A3,n)
b_ref = zeros(n)
d_ref = zeros(n)
if length(idx_ref) == 0
	b_ref[1] = 1.
	d_ref[1] = .8
else
	b_ref[idx_ref] = ones(length(idx_ref)) 
	d_ref[idx_ref] = .8*ones(length(idx_ref))
end
connected = (length(idx_ref) == 0)
ω1 = ω0*rand(n); ω1 .-= mean(ω1)
ϕ2 = .5
ϕ3 = .5

# 0. Reach the sync state
Θ1,dΘ1,iter1 = hyper_k(A2,A3,ω1,zeros(n),ϕ2,ϕ3,h,niter1,1e-5)
θstar = Θ1[:,end]

# 1. Run the system without control
ω2 = ω1 .- mean(dΘ1[:,end])
Θ2,dΘ2,iter2 = hyper_k_gaussian_noise(A2,A3,ω2,Θ1[:,end],ξ0,δt,ϕ2,ϕ3,h,niter2+niter2bis+niter3,-1.)

figure("fig",(12,4))
subplot(2,1,1)
for i in 1:n
        PyPlot.plot(h*(1:(niter2+niter2bis+niter3)),Θ2[i,:],color=cm1((i+n/2)/(1.5*n)))
end
subplot(2,1,2)
for i in 1:n
	PyPlot.plot(h*(1:(niter2+niter2bis)),Θ2[i,1:(niter2+niter2bis)],color=cm1((i+n/2)/(1.5*n)))
end
#=
figure("fig (mod 2π)")
subplot(2,1,1)
for i in 1:n
        PyPlot.plot(h*(1:length(Θ2[i,:])),mod.(Θ2[i,:] .+ π,2π) .- π,color=cm1((i+n/2)/(1.5*n)))
end
subplot(2,1,2)
for i in 1:n
        PyPlot.plot(h*(1:length(Θ2[i,:])),mod.(Θ2[i,:] .+ π,2π) .- π,color=cm1((i+n/2)/(1.5*n)))
end
=#

# 2. Run THIS
X = Θ2[:,1:nstep:niter2]
Xm = mean(X,dims=2)
X .-= Xm
X ./= mean(abs.(X),dims=2)
Y = dΘ2[:,1:nstep:niter2]
ooi = [3,]
dmax = 2

Ainf,coeff,relerr = this(X,Y,ooi,dmax,λ)
l = [Ainf[3][i,1:3] for i in 1:length(Ainf[3][:,1])]
idx = get_loose_ends(Ainf[2],Ainf[3],n)

#=
b = ones(n)
d = .8*ones(n)
idx = get_loose_ends(A2,A3)
=#

b = zeros(n)
d = zeros(n)
if length(idx) == 0
        @info "System is connected ($(connected))"
	b[1] = 1.
	d[1] = .8
else
        @info "System is not connected ($(!connected))"
	b[idx] = ones(length(idx)) 
	d[idx] = .8*ones(length(idx))
end

# 3. Run the damped system
Θ3,dΘ3,iter3 = hyper_k_drooped_gaussian_noise(A2,A3,ω2,Θ2[:,niter2+niter2bis],ρ*b,θstar,ξ0,δt,ϕ2,ϕ3,h,niter3,-1.)
#Θ3,dΘ3,iter3 = hyper_k_damped_gaussian_noise(A2,A3,ω2,Θ2[:,end],d,ξ0,ϕ2,ϕ3,h,1000,-1.)

figure("fig")
subplot(2,1,2)
for i in 1:n-1
        PyPlot.plot(h*(niter2+niter2bis .+ (1:niter3)),Θ3[i,:],color=cm2((i+n/2)/(1.5*n)))
end
PyPlot.plot(h*(niter2+niter2bis .+ (1:niter3)),Θ3[n,:],color=cm2((n+n/2)/(1.5*n)),label="$(Int64(sum(b))) controlled nodes")
#=
figure("fig (mod 2π)")
subplot(2,1,1)
for i in 1:n
        PyPlot.plot(h*(length(Θ2[i,:]) .+ (1:length(Θ3[i,:]))),mod.(Θ3[i,:] .+ π,2π) .- π,color=cm2((i+n/2)/(1.5*n)))
end
=#

# 3bis. Run the reference-damped system
Θ4,dΘ4,iter4 = hyper_k_drooped_gaussian_noise(A2,A3,ω2,Θ2[:,niter2],ρ*b_ref,θstar,ξ0,δt,ϕ2,ϕ3,h,niter3,-1.)
#Θ3,dΘ3,iter3 = hyper_k_damped_gaussian_noise(A2,A3,ω2,Θ2[:,end],d,ξ0,ϕ2,ϕ3,h,1000,-1.)

#=
figure("fig")
subplot(2,1,2)
for i in 1:n-1
        PyPlot.plot(h*(length(Θ2[i,:]) .+ (1:length(Θ4[i,:]))),Θ4[i,:],color=cm3((i+n/2)/(1.5*n)))
end
PyPlot.plot(h*(length(Θ2[n,:]) .+ (1:length(Θ4[n,:]))),Θ4[n,:],color=cm3((n+n/2)/(1.5*n)),label="$(Int64(sum(b_ref))) controlled nodes")
=#
#=
figure("fig (mod 2π)")
subplot(2,1,2)
for i in 1:n
        PyPlot.plot(h*(length(Θ2[i,:]) .+ (1:length(Θ4[i,:]))),mod.(Θ4[i,:] .+ π,2π) .- π,color=cm3((i+n/2)/(1.5*n)))
end
=#

@info "$(length(idx)) controlled nodes out of $(length(idx_ref)) needed"

H2_1 = sum((Θ2 .- θstar).^2)*h/niter2
H2_2 = sum((Θ3 .- θstar).^2)*h/niter3
H2_ref = sum((Θ4 .- θstar).^2)*h/niter3

@info "H2-norm without control: $(H2_1)"
@info "H2-norm with our control: $(H2_2)"
@info "H2-norm with ideal control: $(H2_ref)"

xmin = 0
xmax = h*(niter2+niter3)
ymin = min(minimum(Θ2),minimum(Θ3))#,minimum(Θ4))
ymax = max(maximum(Θ2),maximum(Θ3))#,maximum(Θ4))
dy = ymax-ymin

figure("fig")
subplot(2,1,1)
#PyPlot.plot([h*niter2,h*niter2],[ymin-.05*dy,ymax+.05*dy],"--k")
ylabel("x - x*")
axis([xmin,xmax,ymin-.05*dy,ymax+.05*dy])
#legend()
subplot(2,1,2)
PyPlot.plot([h*niter2,h*niter2],[ymin-.05*dy,ymax+.05*dy],"--k")
PyPlot.plot([h*(niter2+niter2bis),h*(niter2+niter2bis)],[ymin-.05*dy,ymax+.05*dy],"--k")
xlabel("t [a.u.]")
ylabel("x - x*")
axis([xmin,xmax,ymin-.05*dy,ymax+.05*dy])
legend()

ainf2 = Vector{Float64}[]
if length(Ainf[2]) > 0
	ainf2 = [Ainf[2][i,1:2] for i in 1:length(Ainf[2][:,1])]
end
a2 = Vector{Float64}[]
if length(A2) > 0
	a2 = [A2[i,1:2] for i in 1:length(A2[:,1])]
end
ainf3 = Vector{Float64}[]
if length(Ainf[3]) > 0
	ainf3 = [Ainf[3][i,1:3] for i in 1:length(Ainf[3][:,1])]
end
a3 = Vector{Float64}[]
if length(A3) > 0
	a3 = [A3[i,1:3] for i in 1:length(A3[:,1])]
end

tp_2 = 0
fp_2 = 0
for e in ainf2
	if e in a2
		global tp_2 += 1
	else
		global fp_2 += 1
	end
end
tpr_2 = tp_2/max(1,length(a2))
fpr_2 = fp_2/max(1,n*(n-1) - length(a2))
@info "Order 2: TPR = $(tpr_2), FPR = $(fpr_2)"

tp_3 = 0
fp_3 = 0
for e in ainf3
	if e in a3
		global tp_3 += 1
	else
		global fp_3 += 1
	end
end
tpr_3 = tp_3/max(1,length(a3))
fpr_3 = fp_3/max(1,n*(n-1)*(n-2) - length(a3))
@info "Order 3: TPR = $(tpr_3), FPR = $(fpr_3)"

tpr_b = sum((b.*b_ref))/sum(b_ref)
fpr_b = sum(b.*(1 .- b_ref))/(n - sum(b_ref))
@info "Controlled nodes: TPR = $(tpr_b), FPR = $(fpr_b)"

figure("fig")
subplot(2,1,1)
title("TPR_2 = $(round(tpr_2,digits=2)), FPR_2 = $(round(fpr_2,digits=2)), TPR_3 = $(round(tpr_3,digits=2)), FPR_3 = $(round(fpr_3,digits=2)), TPR_b = $(round(tpr_b,digits=2)), FPR_b = $(round(fpr_b,digits=2))")



