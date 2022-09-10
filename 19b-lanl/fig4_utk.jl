using PyPlot, DelimitedFiles, RecipesBase, Shapefile, Colors

include("processing_tools.jl")

 #=
cmap = get_cmap("plasma")
colshift1 = .5
colshift2 = .5
cols = [cmap(1-(i+colshift1)/(2+colshift1+colshift2)) for i in 0:2]
# =#
# #=
cols = [(243,111,33)./255,(0,145,194)./255,(161,0,202)./255]
# =#
gra = (.8,.8,.8,1.)

ntw = "utk"
n = 98
ls = 1:n
ks = 1:1500
#nok = [(3,684),(5,635),(15,771),(15,775),(16,1016),(16,1044),(16,1045),(17,793),(25,451),(25,784),(29,1223),(34,705),(35,141),(35,980),(35,1171),(36,1332),(36,1437),(37,1053),(37,1187),(38,20),(38,53),(39,1445),(42,666),(42,1300),(43,119),(44,1104),(44,1107),(46,923),(54,969),(66,470),(66,617),(66,965),(67,624),(68,76),(69,389),(72,479),(73,1173),(73,1177),(74,1344),(78,12),(80,253),(80,260),(80,267),(86,30),(87,1450),(88,659),(94,377),(94,484),(94,660),(94,1347),(94,1348),(97,825),(98,376),(98,415),(98,1189)]
T = 300.
τ = .1
file = "mysterious_forcingDrift"
id0 = Int64.(vec(readdlm("data_utk/utk2_ids.csv",',')))
fs1 = 51
fs1id = id0[fs1]
#fs1 = 66

# #=
L0 = zeros(length(ls),length(ks))
for i in 1:length(ls)
	@info "i = $i"
	for j in 1:length(ks)
		L0[i,j] = readdlm("data/utk/"*file*"_l0_$(ls[i]).$(ks[j])_obj.csv",',')[1]
	end
end
nL0 = (L0 .- maximum(L0))./(maximum(L0) - minimum(L0))

@info "L0 is loaded."

# =#
##############################################################

figure("UTK",(14,5))

subplot2grid((2,6),(0,0),colspan=3,rowspan=2)

# #=
s_xy, red_i, name, l_xy, lids, lakes = utk_preprocess_boundaries()

for i in 1:length(s_xy)
	s = s_xy[i]
	for j in 1:length(red_i[i])
		PyPlot.fill(s[1][red_i[i][j]],s[2][red_i[i][j]],color=(.95,.95,.95,1.))
		PyPlot.plot(s[1][red_i[i][j]],s[2][red_i[i][j]],"-k",linewidth=.5)
	end
end

for i in 1:length(l_xy)
	l = l_xy[i]
	ids = lids[i]
	PyPlot.fill(l[1][ids],l[2][ids],color=(1.,1.,1.,1.))
	PyPlot.plot(l[1][ids],l[2][ids],"-k",linewidth=.5)
end

# =#

xy = readdlm("data_utk/utk1_coord.csv",',')[id0,:]
ids = (1:size(xy)[1]).*(xy[:,1] .< -55.)
ids = setdiff(ids,[0,])
PyPlot.plot(xy[ids,1],xy[ids,2],"ok",markersize=5.)
PyPlot.plot(xy[fs1,1],xy[fs1,2],"o",color=cols[1],markersize=10.)
axis([-105.,-63.,24.,48.])
PyPlot.xticks([])
PyPlot.yticks([])


subplot2grid((2,6),(0,3),colspan=3,rowspan=2)

L0red = nL0[[1:fs1-1;fs1+1:size(nL0)[1]],:]
L0max = [maximum(L0red[:,i]) for i in 1:size(L0red)[2]]
L0min = [minimum(L0red[:,i]) for i in 1:size(L0red)[2]]
PyPlot.fill([ks;ks[end:-1:1]]/T,-[L0max;L0min[end:-1:1]],color=gra)
PyPlot.plot(ks/T,-nL0[fs1,:],"-",color=cols[1],linewidth=2.,label="$fs1id")
axis([0.,maximum(ks)/T,-.1,1.1])
xlabel("freq")
ylabel("normalized log-likelihodd")
legend()

