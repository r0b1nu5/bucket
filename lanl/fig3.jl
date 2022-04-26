using PyPlot, DelimitedFiles, RecipesBase, Shapefile

cmap = get_cmap("plasma")
colshift = .5
cols = [cmap(1-(i+colshift)/(2+colshift)) for i in 0:2]

# #=
ntw = "utk"
n = 99
ls = 1:n
ks = 1:1500
T = 300.
file = "mysterious_forcing"
fs1 = 51
fs2 = 66
fs3 = 85
# =#

L0 = zeros(length(ls),length(ks))
for i in 1:length(ls)
	for j in 1:length(ks)
		L0[i,j] = readdlm("data/"*file*"_l0_$(ls[i]).$(ks[j])_obj.csv",',')[1]
	end
end
nL0 = (L0 .- maximum(L0))./(maximum(L0) - minimum(L0))

figure("fig 4",(19,6))

subplot(1,2,1)

states = Shapefile.Table("data_utk/states_shapes/tl_2021_us_state.shp")
skip = ["02","15","60","66","69","72","78"]

ss = empty(Shapefile.shapes(states))
for row in states
	if !(row.STATEFP in skip)
		push!(ss, Shapefile.shape(row))
	end
end

s_xy = Vector{Tuple{Vector{Float64},Vector{Float64}}}()
for s in ss
	x = Vector{Float64}()
	y = Vector{Float64}()
	for p in s.points
		push!(x,p.x)
		push!(y,p.y)
	end
	push!(s_xy,(x,y))
end

for s in s_xy
	PyPlot.plot(s[1],s[2],"-k",linewidth=.5)
end

xy = readdlm("data_utk/utk1_coord.csv",',')
ids = (1:size(xy)[1]).*(xy[:,1] .< -65.)
ids = setdiff(ids,[0,])
PyPlot.plot(xy[ids,1],xy[ids,2],"ok",markersize=5.)
PyPlot.plot(xy[fs1,1],xy[fs1,2],"o",color=cols[1],markersize=10.)
PyPlot.plot(xy[fs2,1],xy[fs2,2],"o",color=cols[3],markersize=10.)
PyPlot.plot(xy[fs3,1],xy[fs3,2],"o",color=cols[2],markersize=10.)
axis([-105.,-65.,24.,48.])
PyPlot.xticks([])
PyPlot.yticks([])

subplot(1,2,2)

L0red = nL0[[1:fs1-1;fs1+1:fs2-1;fs2+1:fs3-1;fs3+1:size(nL0)[1]],:]
L0max = [maximum(L0red[:,i]) for i in 1:size(L0red)[2]]
L0min = [minimum(L0red[:,i]) for i in 1:size(L0red)[2]]
PyPlot.fill([ks;ks[end:-1:1]]/T,[L0max;L0min[end:-1:1]],color="C7")
PyPlot.plot(ks/T,nL0[fs3,:],"-",color=cols[2],linewidth=2.,label="$fs3")
PyPlot.plot(ks/T,nL0[fs1,:],"-",color=cols[1],linewidth=2.,label="$fs1")
PyPlot.plot(ks/T,nL0[fs2,:],"-",color=cols[3],linewidth=2.,label="$fs2")
axis([0.,maximum(ks)/T,-1.1,.1])
xlabel("freq")
ylabel("normalized inverse log-likelihodd")
legend()



#=
subplot2grid((2,4),(0,2),colspan=1,rowspan=2)
for i in 1:length(ls)
	PyPlot.plot(ks/T,L0[i,:],"-o",color=cmap((i-1)/(length(ls)-1)))
end
xlabel("freq")
ylabel("obj")

subplot2grid((2,4),(0,3),colspan=1,rowspan=1)
PyPlot.plot(ks/T,L1,"-o",color=cmap(.4))
xlabel("freq")
ylabel("obj")

subplot2grid((2,4),(1,3),colspan=1,rowspan=1)
PyPlot.plot(ls,γ1[iii],"-o",color=cmap(0.))
xlabel("node id")
ylabel("amplitude")

