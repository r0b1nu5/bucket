using LinearAlgebra, PyPlot, DelimitedFiles, Statistics

include("line_dist.jl")
# #=
#L = readdlm("ntws_data/uk_lap_mat.csv",',')
Lsp = readdlm("ntws_data/hamsterster2_lap_mat_sp.csv",',')
L = sparse(Int.(Lsp[:,1]),Int.(Lsp[:,2]),Lsp[:,3])

#h = .01
h = .001

B,w,Bt = L2B(L)
W = diagm(0 => w)

n = size(L)[1]

qs = [1.,2.,3.]
Ks = [1.,1.,1.]

# #=
om = rand(n)
om .-= mean(om)
# =#

th0 = rand(n)
#th0 = pinv(L)*om

ls = [51,20]
#as = [5.,1.]
as = [5.,1.]
Ti = [1200*h,1000*h]
Tf = [1200*h,2000*h]

#sig = .0005
sig = .0005

ths0,dh1,ds = kuramoto_q(L,qs,Ks,.2*om,th0,ls,zeros(length(as)),Ti,Tf,0.,true,20000,1e-4,h)
th1 = ths0[:,end]

ths,dhs,ds = kuramoto_q(L,qs,Ks,.2*om,th1,ls,as,Ti,Tf,sig,true,3000,-1.,h)
psi = L*ths

T = size(ths)[2]
# =#
t = Array{Float64,2}(undef,2*length(ls),0)
d = Array{Float64,2}(undef,2*length(ls),0)
p = Array{Float64,2}(undef,2*length(ls),0)
ts = Array{Float64,2}(undef,n-2*length(ls),0)
ds = Array{Float64,2}(undef,n-2*length(ls),0)
ps = Array{Float64,2}(undef,n-2*length(ls),0)

ii = Array{Int64,1}()
jj = Array{Int64,1}()

for l in ls
	push!(ii,Int(maximum(B[:,l].*(1:n))))
	push!(jj,-Int(minimum(B[:,l].*(1:n))))
end
ids = setdiff(Array(1:n),[ii;jj])

for i in 1:T
	global t = [t ths[[ii;jj],i]-ths[[ii;jj],1]]
	global ts = [ts sort(ths[ids,i]-ths[ids,1])]
	global d = [d dhs[[ii;jj],i]-dhs[[ii;jj],1]]
	global ds = [ds sort(dhs[ids,i]-dhs[ids,1])]
	global p = [p psi[[ii;jj],i]-psi[[ii;jj],1]]
	global ps = [ps sort(psi[ids,i]-psi[ids,1])]
end

@info "==============================================="
@info "$(sin.(Bt[ls,:]*th1*qs')*Ks)"


# =#

#=
mit = [minimum(t[ids,i]) for i in 1:T]
met = [median(t[ids,i]) for i in 1:T]
mat = [maximum(t[ids,i]) for i in 1:T]
mid = [minimum(d[ids,i]) for i in 1:T]
med = [median(d[ids,i]) for i in 1:T]
mad = [maximum(d[ids,i]) for i in 1:T]
mip = [minimum(p[ids,i]) for i in 1:T]
mep = [median(p[ids,i]) for i in 1:T]
map = [maximum(p[ids,i]) for i in 1:T]
=#
mit = ts[1,:]
met = ts[round(Int,n/2-1),:]
mat = ts[end,:]
mid = ds[1,:]
med = ds[round(Int,n/2-1),:]
mad = ds[end,:]
mip = ps[1,:]
mep = ps[round(Int,n/2-1),:]
map = ps[end,:]


T1 = 1
T2 = T
# #= Reduced time span
T1 = 500
T2 = 3000
# =#

figure()


subplot(3,1,1)
PyPlot.fill(h*[T1:T2;T2:-1:T1],[mit[T1:T2];mat[T2:-1:T1]],"k",alpha=.3)
PyPlot.plot(h*(T1:T2),met[T1:T2],"k")
for i in [ii;jj]
	PyPlot.plot(h*(T1:T2),t[i,T1:T2])
end
ylabel("θ")
subplot(3,1,2)
PyPlot.fill(h*[T1:T2;T2:-1:T1],[mid[T1:T2];mad[T2:-1:T1]],"k",alpha=.3)
PyPlot.plot(h*(T1:T2),med[T1:T2],"k")
for i in [ii;jj]
	PyPlot.plot(h*(T1:T2),d[i,T1:T2])
end
ylabel("θ'")
subplot(3,1,3)
PyPlot.fill(h*[T1:T2;T2:-1:T1],[mip[T1:T2];map[T2:-1:T1]],"k",alpha=.3)
PyPlot.plot(h*(T1:T2),mep[T1:T2],"k")
for i in [ii;jj]
	PyPlot.plot(h*(T1:T2),p[i,T1:T2])
end
xlabel("t")
ylabel("ψ")





