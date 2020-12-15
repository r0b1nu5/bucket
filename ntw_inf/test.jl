using LinearAlgebra, PyPlot, DelimitedFiles, Statistics, Dates

include("line_dist.jl")
# #=
#L = readdlm("ntws_data/uk_lap_mat.csv",',')
#Lsp = readdlm("ntws_data/hamsterster2_lap_mat_sp.csv",',')
#Lsp = readdlm("ntws_data/euroroad_red_lap_mat_sp.csv",',')
Lsp = readdlm("ntws_data/polblogs_red_lap_mat_sp.csv",',')
L = sparse(Int.(Lsp[:,1]),Int.(Lsp[:,2]),Lsp[:,3])

#h = .01
h = .001

B,w,Bt = L2B(L)
W = diagm(0 => w)

n = size(L)[1]

qs = [1.,2.,3.]
Ks = [1.,1.,1.]

 #=
om0 = rand(n)
om0 .-= mean(om0)
# =#
om = 0.5 * om0

th0 = rand(n)
#th0 = pinv(L)*om

ls = [rand(1:16000),rand(1:16000)]
ls = [9840,1047]
as = [1.,1.5]
#as = [1.,1.]
Ti = [1200*h,1000*h]
Tf = [1200*h,2000*h]

tau = [.5,.2]

sig = .001

ths0,dh1,ds = kuramoto_q(L,qs,Ks,.2*om,th0,ls,zeros(length(as)),Ti,Tf,0.,true,20000,1e-4,h)
#ths0,dh1,prt = linear_noise(L,om,th0,ls,zeros(length(as)),tau,sig,true,20000,1e-4,h)
th1 = ths0[:,end]

ths,dhs,ds = kuramoto_q(L,qs,Ks,.2*om,th1,ls,as,Ti,Tf,sig,true,3000,-1.,h)
#ths,dhs,prt = linear_noise(L,om,th1,ls,as,tau,sig,true,3000,-1.,h)

@info "Compute ψ"
psi = L*ths

T = size(ths)[2]

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



@info "Sorted time series"
io_t = open("data1/t.csv","a")
io_ts = open("data1/ts.csv","a")
io_d = open("data1/d.csv","a")
io_ds = open("data1/ds.csv","a")
io_p = open("data1/p.csv","a")
io_ps = open("data1/ps.csv","a")
for i in 1:T
	if i%1000 == 0
		@info "$(now()) -- t = $i/$T"
	end
	dths = ths[:,i]-ths[:,1]
	writedlm(io_t,dths[[ii;jj]]',",")
	writedlm(io_ts,[minimum(dths[ids]) median(dths[ids]) maximum(dths[ids])],",")
	ddhs = dhs[:,i]-dhs[:,1]
	writedlm(io_d,ddhs[[ii;jj]]',",")
	writedlm(io_ds,[minimum(ddhs[ids]) median(ddhs[ids]) maximum(ddhs[ids])],",")
	dps = psi[:,i]-psi[:,1]
	writedlm(io_p,dps[[ii;jj]]',",")
	writedlm(io_ps,[minimum(dps[ids]) median(dps[ids]) maximum(dps[ids])],",")
#	global t = [t ths[[ii;jj],i]-ths[[ii;jj],1]]
#	global ts = [ts sort(ths[ids,i]-ths[ids,1])]
#	global d = [d dhs[[ii;jj],i]-dhs[[ii;jj],1]]
#	global ds = [ds sort(dhs[ids,i]-dhs[ids,1])]
#	global p = [p psi[[ii;jj],i]-psi[[ii;jj],1]]
#	global ps = [ps sort(psi[ids,i]-psi[ids,1])]
end

close(io_t)
close(io_ts)
close(io_d)
close(io_ds)
close(io_p)
close(io_ps)

t = Array(readdlm("data1/t.csv",',')')
ts = Array(readdlm("data1/ts.csv",',')')
d = Array(readdlm("data1/d.csv",',')')
ds = Array(readdlm("data1/ds.csv",',')')
p = Array(readdlm("data1/p.csv",',')')
ps = Array(readdlm("data1/ps.csv",',')')

rm("data1/t.csv")
rm("data1/ts.csv")
rm("data1/d.csv")
rm("data1/ds.csv")
rm("data1/p.csv")
rm("data1/ps.csv")


@info "Time series to be displayed"
mit = ts[1,:]
#met = ts[round(Int,n/2-1),:]
met = ts[2,:]
mat = ts[end,:]
mid = ds[1,:]
#med = ds[round(Int,n/2-1),:]
med = ds[2,:]
mad = ds[end,:]
mip = ps[1,:]
#mep = ps[round(Int,n/2-1),:]
mep = ps[2,:]
map = ps[end,:]


T1 = 1
T2 = T
# #= Reduced time span
T1 = 500
T2 = 3000
# =#

# =#
figure()

subplot(3,1,1)
PyPlot.fill(h*[T1:T2;T2:-1:T1],[mit[T1:T2];mat[T2:-1:T1]],"k",alpha=.3)
PyPlot.plot(h*(T1:T2),met[T1:T2],"k")
for i in 1:length([ii;jj])
	PyPlot.plot(h*(T1:T2),t[i,T1:T2])
end
ylabel("θ")
subplot(3,1,2)
PyPlot.fill(h*[T1:T2;T2:-1:T1],[mid[T1:T2];mad[T2:-1:T1]],"k",alpha=.3)
PyPlot.plot(h*(T1:T2),med[T1:T2],"k")
for i in 1:length([ii;jj])
	PyPlot.plot(h*(T1:T2),d[i,T1:T2])
end
ylabel("θ'")
subplot(3,1,3)
PyPlot.fill(h*[T1:T2;T2:-1:T1],[mip[T1:T2];map[T2:-1:T1]],"k",alpha=.3)
PyPlot.plot(h*(T1:T2),mep[T1:T2],"k")
for i in 1:length([ii;jj])
	PyPlot.plot(h*(T1:T2),p[i,T1:T2])
end
xlabel("t")
ylabel("ψ")





