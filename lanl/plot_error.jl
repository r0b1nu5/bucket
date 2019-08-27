using PyPlot

function plot_error(Lm::Array{Float64,2}, Lh::Array{Float64,2}, dm::Array{Float64,1}, dh::Array{Float64,1}, a::Array{Float64,1}, ah::Array{Float64,1}, f::Array{Float64,1}, fh::Array{Float64,1})
	n = length(dm)

	figure()
	
	subplot(1,4,1)
	PyPlot.plot([min(minimum(Lm),minimum(Lh)),max(maximum(Lm),maximum(Lh))],[min(minimum(Lm),minimum(Lh)),max(maximum(Lm),maximum(Lh))],"--k")
	subplot(1,4,2)
	PyPlot.plot([min(minimum(dm),minimum(dh)),max(maximum(dm),maximum(dh))],[min(minimum(dm),minimum(dh)),max(maximum(dm),maximum(dh))],"--k")
	subplot(1,4,3)
	PyPlot.plot([min(minimum(a),minimum(ah)),max(maximum(a),maximum(ah))],[min(minimum(a),minimum(ah)),max(maximum(a),maximum(ah))],"--k")
	subplot(1,4,4)
	PyPlot.plot([min(minimum(f),minimum(fh)),max(maximum(f),maximum(fh))],[min(minimum(f),minimum(fh)),max(maximum(f),maximum(fh))],"--k")
	for i in 1:n
		subplot(1,4,1)
		for j in 1:n
			PyPlot.plot(Lm[i,j],Lh[i,j],"o",color="C0")
			PyPlot.text(Lm[i,j],Lh[i,j],"($i,$j)")
		end
		xlabel("L")
		ylabel("Lh")

		subplot(1,4,2)
		PyPlot.plot(dm[i],dh[i],"o",color="C1")
		PyPlot.text(dm[i],dh[i],"$i")
		xlabel("d")
		ylabel("dh")

		subplot(1,4,3)
		PyPlot.plot(a[i],ah[i],"o",color="C2")
		PyPlot.text(a[i],ah[i],"$i")
		xlabel("a")
		ylabel("ah")

		subplot(1,4,4)
		PyPlot.plot(f[i],fh[i],"o",color="C3")
		PyPlot.text(f[i],fh[i],"$i")
		xlabel("f")
		ylabel("fh")
	end
end




