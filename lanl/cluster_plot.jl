using PyPlot, DelimitedFiles

dates = ["2013-01-15_00", "2013-03-10_04", "2013-04-03_02", "2013-04-03_03", "2013-04-03_07", "2013-07-30_01", "2013-07-30_04", "2013-07-30_09"]

ntw = "uk" # ["ntw5", "ntw10", "ntw20", "uk10", "pen"] 
n = 120 # [5, 10, 20, 120, ???]

if ntw == "pen"
	AA = Array{Float64,2}(undef,3,0)
	i = 0
	for date in dates
		i += 1

		ids = Int.(vec(readdlm("data_PEN/pen_"*date*"_ids.csv",',')))
		n = length(ids)

		ah = vec(readdlm("data_PEN/pen_"*date*"_ah.csv",','))

		AA = sortslices([abs.(ah) 1:n],dims=1,rev=true)
		PP = AA[:,1]./sum(AA[:,1])
		id1 = Int(AA[1,2])
		id2 = Int(AA[2,2])
		id3 = Int(AA[3,2])

		figure(666)
		PyPlot.plot([i,i],[0,PP[1]],color="C0")
		PyPlot.plot([i,i],[0,PP[2]],color="C1")
		PyPlot.plot([i,i],[0,PP[3]],color="C2")

		PyPlot.text(i,PP[1]+.05,"($id1,$id2,$id3)",rotation="270")
	end
	xticks(1:i,dates,rotation=20)

	xlabel("Dates")
	ylabel("Confidence")
else
	reL = Array{Float64,1}()
	rea = Array{Float64,1}()
	ref = Array{Float64,1}()
	rep = Array{Float64,1}()
	
	ah = Array{Float64,2}(undef,120,0)
	
	for i in 1:n
		push!(reL,readdlm("data/"*ntw*"_reL_$(i).csv",',')[1])
		push!(rea,readdlm("data/"*ntw*"_rea_$(i).csv",',')[1])
		push!(ref,readdlm("data/"*ntw*"_ref_$(i).csv",',')[1])
		push!(rep,readdlm("data/"*ntw*"_rep_$(i).csv",',')[1])
	
		ah = [ah vec(readdlm("data/"*ntw*"ah_$(i).csv",','))]
	
		AA = sortslices([abs.(ah[:,end]) 1:n],dims=1,rev=true)
	        PP = AA[:,1]./sum(AA[:,1])
	        id1 = Int(AA[1,2])
	        id2 = Int(AA[2,2])
	        id3 = Int(AA[3,2])
	
	        figure(666)
	        subplot(1,3,3)
	        PyPlot.plot([i,i],[0,abs(ah[id1])/sum(abs.(ah))],color="C0")
	        PyPlot.plot([i,i],[0,abs(ah[id2])/sum(abs.(ah))],color="C1")
	        PyPlot.plot([i,i],[0,abs(ah[id3])/sum(abs.(ah))],color="C2")
	        PyPlot.plot(i,abs(ah[i])/sum(abs.(ah)),"o",color="C3",markersize=4.)
	end
	
	figure(666)
	
	subplot(1,3,1)
	plot_ntw(ntw)
	
	subplot(1,3,2)
	PyPlot.semilogy(1:n,reL,"o")
	PyPlot.semilogy(n .+ (1:n),rea,"o")
	PyPlot.semilogy(2*n .+ (1:n),ref,"o")
	PyPlot.semilogy(3*n .+ (1:n),rep,"o")
end




