using PyPlot, DelimitedFiles

dates = ["2013-01-15_00", "2013-03-10_04", "2013-04-03_02", "2013-04-03_03", "2013-04-03_07", "2013-07-30_01", "2013-07-30_04", "2013-07-30_09"]

cols = ["C0","C1","C2","C3","C4","C5","C6","C7","C8","C9","C10"]
c = 0

figure()
mah = zeros(230)
k = zeros(230)

for date in dates
	global c,mah,k
	c += 1
	
	ah = vec(readdlm("data_PEN/pen_$(date)_ah.csv",','))
	ids = Int.(vec(readdlm("data_PEN/pen_$(date)_final_ids.csv",',')))

	p = abs.(ah)./sum(abs.(ah))

	mah[ids] += p
	k[ids] .+= 1

	subplot(3,3,c)
	PyPlot.plot(ids,p,"o",color=cols[c])
	title(date)
end

iids = Int.(setdiff((k .!= 0).*(1:230),[0,]))

mp = mah[iids]./k[iids]

subplot(3,3,c+1)
PyPlot.plot(iids,mp,"o",color=cols[c+1])
title("Average confidence")
	



	



