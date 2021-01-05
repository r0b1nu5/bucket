using LinearAlgebra, PyPlot, DelimitedFiles, Statistics, Dates

include("line_dist.jl")
# #=
#L = readdlm("ntws_data/uk_lap_mat.csv",',')
#Lsp = readdlm("ntws_data/hamsterster2_lap_mat_sp.csv",',')

 #= airports
Lsp = readdlm("ntws_data/usairports_lap_mat_sp.csv",',')
om0 = vec(readdlm("line_dist_data/usairports_om1.csv",','))
om = om0
#ls = rand(1:17214,2)
ls = [9850,8885]
as = [.9,.8]
sig = .0005
# =#

 #= polblogs
Lsp = readdlm("ntws_data/polblogs_redred_lap_mat_sp.csv",',')
om0 = vec(readdlm("line_dist_data/polblogs_om2.csv",','))
om = .5 * om0
ls = [rand(1:16000),rand(1:16000)]
#ls = [9838,1047]
as = [.3,.4]
sig = .0005
# =#

# #= euroroad
Lsp = readdlm("ntws_data/euroroad_red_lap_mat_sp.csv",',')
om0 = vec(readdlm("line_dist_data/euroroad_om1.csv",','))
om = .1 * om0
#ls = [rand(1:1305),rand(1:1305)]
ls = [454,613]
as = [.5,.5]
sig = .005
# =#

L = sparse(Int.(Lsp[:,1]),Int.(Lsp[:,2]),Lsp[:,3])

#h = .01
h = .001

B,w,Bt = L2B(L)
W = diagm(0 => w)

n = size(L)[1]

qs = [1.,2.,3.]
Ks = [1.,1.,1.]

th0 = rand(n)
#th0 = pinv(L)*om

Ti = [1200*h,1000*h]
Tf = [1200*h,2000*h]

tau = [.8,1.]


ths0,dh1,ds = kuramoto_q(L,qs,Ks,om,th0,ls,zeros(length(as)),Ti,Tf,0.,true,20000,1e-4,h) 
#ths0,dh1,prt = linear_noise(L,om,th0,ls,zeros(length(as)),tau,0.,true,20000,1e-4,h)
th1 = ths0[:,end]

ths,dhs,ds = kuramoto_q(L,qs,Ks,om,th1,ls,as,Ti,Tf,sig,true,3000,-1.,h) 
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
	writedlm(io_ts,[minimum(dths[ids]) quantile(dths[ids],.05) quantile(dths[ids],.25) median(dths[ids]) quantile(dths[ids],.75) quantile(dths[ids],.95) maximum(dths[ids])],",")
	ddhs = dhs[:,i]-dhs[:,1]
	writedlm(io_d,ddhs[[ii;jj]]',",")
	writedlm(io_ds,[minimum(ddhs[ids]) quantile(ddhs[ids],.05) quantile(ddhs[ids],.25) median(ddhs[ids]) quantile(ddhs[ids],.75) quantile(ddhs[ids],.95) maximum(ddhs[ids])],",")
	dps = psi[:,i]-psi[:,1]
	writedlm(io_p,dps[[ii;jj]]',",")
	writedlm(io_ps,[minimum(dps[ids]) quantile(dps[ids],.05) quantile(dps[ids],.25) median(dps[ids]) quantile(dps[ids],.75) quantile(dps[ids],.95) maximum(dps[ids])],",")
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

writedlm("data1/t_end.csv",ths[:,end]-ths[:,1],',')
writedlm("data1/d_end.csv",dhs[:,end]-dhs[:,1],',')
writedlm("data1/p_end.csv",psi[:,end]-psi[:,1],',')

t = Array(readdlm("data1/t.csv",',')')
ts = Array(readdlm("data1/ts.csv",',')')
te = ths[:,end]-ths[:,1]
d = Array(readdlm("data1/d.csv",',')')
ds = Array(readdlm("data1/ds.csv",',')')
de = dhs[:,end]-dhs[:,1]
p = Array(readdlm("data1/p.csv",',')')
ps = Array(readdlm("data1/ps.csv",',')')
pe = psi[:,end]-psi[:,1]

# #=
rm("data1/t.csv")
rm("data1/ts.csv")
rm("data1/t_end.csv")
rm("data1/d.csv")
rm("data1/ds.csv")
rm("data1/d_end.csv")
rm("data1/p.csv")
rm("data1/ps.csv")
rm("data1/p_end.csv")
# =#


@info "Time series to be displayed"
mit = ts[1,:]
d05t = ts[2,:]
d25t = ts[3,:]
met = ts[4,:]
d75t = ts[5,:]
d95t = ts[6,:]
mat = ts[7,:]
mid = ds[1,:]
d05d = ds[2,:]
d25d = ds[3,:]
med = ds[4,:]
d75d = ds[5,:]
d95d = ds[6,:]
mad = ds[7,:]
mip = ps[1,:]
d05p = ps[2,:]
d25p = ps[3,:]
mep = ps[4,:]
d75p = ps[5,:]
d95p = ps[6,:]
map = ps[7,:]


T1 = 1
T2 = T
# #= Reduced time span
T1 = 500
T2 = 3000
# =#

# =#
figure("time series")

subplot(3,1,1)
PyPlot.fill(h*[T1:T2;T2:-1:T1],[mit[T1:T2];mat[T2:-1:T1]],"k",alpha=.2)
PyPlot.fill(h*[T1:T2;T2:-1:T1],[d05t[T1:T2];d95t[T2:-1:T1]],"k",alpha=.3)
PyPlot.fill(h*[T1:T2;T2:-1:T1],[d25t[T1:T2];d75t[T2:-1:T1]],"k",alpha=.3)
PyPlot.plot(h*(T1:T2),met[T1:T2],"k")
for i in 1:length([ii;jj])
	PyPlot.plot(h*(T1:T2),t[i,T1:T2])
end
ylabel("θ")
subplot(3,1,2)
PyPlot.fill(h*[T1:T2;T2:-1:T1],[mid[T1:T2];mad[T2:-1:T1]],"k",alpha=.2)
PyPlot.fill(h*[T1:T2;T2:-1:T1],[d05d[T1:T2];d95d[T2:-1:T1]],"k",alpha=.3)
PyPlot.fill(h*[T1:T2;T2:-1:T1],[d25d[T1:T2];d75d[T2:-1:T1]],"k",alpha=.3)
PyPlot.plot(h*(T1:T2),med[T1:T2],"k")
for i in 1:length([ii;jj])
	PyPlot.plot(h*(T1:T2),d[i,T1:T2])
end
ylabel("θ'")
subplot(3,1,3)
PyPlot.fill(h*[T1:T2;T2:-1:T1],[mip[T1:T2];map[T2:-1:T1]],"k",alpha=.2)
PyPlot.fill(h*[T1:T2;T2:-1:T1],[d05p[T1:T2];d95p[T2:-1:T1]],"k",alpha=.3)
PyPlot.fill(h*[T1:T2;T2:-1:T1],[d25p[T1:T2];d75p[T2:-1:T1]],"k",alpha=.3)
PyPlot.plot(h*(T1:T2),mep[T1:T2],"k")
for i in 1:length([ii;jj])
	PyPlot.plot(h*(T1:T2),p[i,T1:T2])
end
xlabel("t")
ylabel("ψ")

figure("snapshot")

subplot(3,1,1)
PyPlot.plot(ids,te[ids],".",color="C7")
for i in [ii;jj]
	PyPlot.plot(i,te[i],"x")
end
ylabel("θ")
subplot(3,1,2)
PyPlot.plot(ids,de[ids],".",color="C7")
for i in [ii;jj]
	PyPlot.plot(i,de[i],"x")
end
ylabel("θ'")
subplot(3,1,3)
PyPlot.plot(ids,pe[ids],".",color="C7")
for i in [ii;jj]
	PyPlot.plot(i,pe[i],"x")
end
ylabel("ψ")
xlabel("i")



