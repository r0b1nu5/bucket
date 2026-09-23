using PyPlot

include("gen-hypergraph.jl")
include("hyper-lv.jl")
include("tools.jl")

include("../../THIS/this.jl")

n = 10
p = .1
r0 = 1.
l0 = 1.
ξ0 = 2.
ρ = 5.
λ = .005
λ = 2e-5
zer0 = 1e-10

h = .01
nstep = 5
δt = nstep*h
niter1 = 10000
niter1 = 100
#niter2 = 5000; niter3 = 5000
#niter2 = 2000; niter3 = 8000
niter2 = 1000; niter3 = 2000
niter2 = 30000; niter2bis = 500; niter3 = 30000
#niter2 = 1500; niter2bis = 150; niter3 = 4350

cm1 = get_cmap("Greens")
cm2 = get_cmap("Purples")
cm3 = get_cmap("Oranges")

A2 = zeros(0,3)
A3,B,B2,E = rand_3_digraph(n,p)
le = size(A3)[1]
#A3[:,4] = sign.(3*rand(le) .- 2) .* .003 .* rand(le)
A3[:,4] = -.003 .* rand(le)
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
r1 = r0*rand(n); r1 .-= mean(r1)
l1 = ones(n)

# 0. Reach the steady state
X1,dX1,iter1 = hyper_lv(A2,A3,r1,l1,rand(n),h,niter1,-1e-5)
xstar = X1[:,end]

# 1. Run the system without control
X2,dX2,iter2 = hyper_lv_gaussian_noise(A2,A3,r1,l1,X1[:,end],ξ0,δt,h,niter2+niter2bis+niter3,-1.,zer0,true)

figure("fig",(12,4))
subplot(2,1,1)
for i in 1:n
#        PyPlot.plot(h*(1:(niter2+niter2bis+niter3)),X2[i,:],color=cm1((i+n/2)/(1.5*n)))
        PyPlot.plot(h*(1:iter2),X2[i,:],color=cm1((i+n/2)/(1.5*n)))
end

subplot(2,1,2)
for i in 1:n
	PyPlot.plot(h*(1:(niter2+niter2bis)),X2[i,1:(niter2+niter2bis)],color=cm1((i+n/2)/(1.5*n)))
end

#=
subplot(3,1,3)
PyPlot.plot(h*(1000:(niter2+niter2bis+niter3)),[sum(minimum(X2[:,j-999:j],dims=2) .< zer0) for j in 1000:(niter2+niter2bis+niter3)],color=cm1(.8))
=#

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
XX = X2[:,1:nstep:niter2]
XXm = mean(XX,dims=2)
XX .-= XXm
mXX = mean(abs.(XX),dims=2)
XX ./= mXX
YY = dX2[:,1:nstep:niter2]
YY ./= mXX
ooi = [3,]
dmax = 2

Ainf,coeff,relerr = this(XX,YY,ooi,dmax,λ)
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
X3,dX3,iter3 = hyper_lv_drooped_gaussian_noise(A2,A3,r1,l1,X2[:,niter2+niter2bis],ρ*b,xstar,ξ0,δt,h,niter3,-1.)
#Θ3,dΘ3,iter3 = hyper_k_damped_gaussian_noise(A2,A3,ω2,Θ2[:,end],d,ξ0,ϕ2,ϕ3,h,1000,-1.)

figure("fig")
subplot(2,1,2)
for i in 1:n-1
        PyPlot.plot(h*(niter2+niter2bis .+ (1:niter3)),X3[i,:],color=cm2((i+n/2)/(1.5*n)))
end
PyPlot.plot(h*(niter2+niter2bis .+ (1:niter3)),X3[n,:],color=cm2((n+n/2)/(1.5*n)),label="$(Int64(sum(b))) controlled nodes")

#=
subplot(3,1,3)
PyPlot.plot(h*(niter2+niter2bis .+ (1000:niter3)),[sum(minimum(X3[:,j-999:j],dims=2) .< zer0) for j in 1000:niter3],color=cm2(.8))
=#

#=
figure("fig (mod 2π)")
subplot(2,1,1)
for i in 1:n
        PyPlot.plot(h*(length(Θ2[i,:]) .+ (1:length(Θ3[i,:]))),mod.(Θ3[i,:] .+ π,2π) .- π,color=cm2((i+n/2)/(1.5*n)))
end
=#

# 3bis. Run the reference-damped system
X4,dX4,iter4 = hyper_lv_drooped_gaussian_noise(A2,A3,r1,l1,X2[:,niter2],ρ*b_ref,xstar,ξ0,δt,h,niter3,-1.)
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

H2_1 = sum((X2 .- xstar).^2)*h/niter2
H2_2 = sum((X3 .- xstar).^2)*h/niter3
H2_ref = sum((X4 .- xstar).^2)*h/niter3

@info "H2-norm without control: $(H2_1)"
@info "H2-norm with our control: $(H2_2)"
@info "H2-norm with ideal control: $(H2_ref)"

xmin = 0
xmax = h*(niter2+niter3)
ymin = min(minimum(X2),minimum(X3))#,minimum(Θ4))
ymax = min(max(maximum(X2),maximum(X3)),2*maximum(X3))#,maximum(Θ4))
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
#=
subplot(3,1,3)
ylabel("# surviving species")
xlabel("t [a.u.]")
xlim(xmin,xmax)
=#

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



