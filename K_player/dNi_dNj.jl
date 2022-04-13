using PyPlot, DelimitedFiles, Statistics, Dates

include("journals.jl")

gen_data = false # otherwise, load data

ref_year = "00"

#journal = "science"; col = "C2"; lw = 2.
#journal = "prl"; col = "C3"; lw = 1.
journal = "lancet"; col = "C0"; lw = 2.

yi = 99
yf = 8

Nmin = 1
Nmax = 200000

if gen_data
	xi = readdlm("data/time_split/"*journal*"_$(ref_year)_$yi.txt",'\t')[2:end-2,1:2]
	auti = xi[:,1]
	numi = Int64.(xi[:,2])
	
	xf = readdlm("data/time_split/"*journal*"_$(ref_year)_$yf.txt",'\t')[2:end-2,1:2]
	autf = xf[:,1]
	numf = Int64.(xf[:,2])
	
	as = Vector{String}()
	ni = Vector{Int64}()
	nf = Vector{Int64}()
	t = now()
	for k in 1:length(auti)
		if k%100 == 0
			global t
			@info "$k/$(length(auti)), t = $(now() - t)"
			t = now()
		end
		if Nmin .<= numi[k] .<= Nmax
			idx = findmax(autf .== auti[k])[2]
			if numi[k] < numf[idx]
				push!(as,auti[k])
				push!(ni,numi[k])
				push!(nf,numf[idx])
			end
		end
	end

	writedlm("data/time_split/"*journal*"_cadv_ni.csv",ni,',')
	writedlm("data/time_split/"*journal*"_cadv_nf.csv",nf,',')
else
	ni = vec(readdlm("data/time_split/"*journal*"_cadv_ni.csv",','))
	nf = vec(readdlm("data/time_split/"*journal*"_cadv_nf.csv",','))
end

mi = union(ni)
Rs = Dict{Int64,Vector{Float64}}(n => Float64[] for n in mi)
ρs = Vector{Float64}()
for i in 1:length(ni)
	ρ = ni[i]/nf[i]
	push!(ρs,ρ)
	push!(Rs[ni[i]],ρ)
end
Qs = Dict{Int64,Vector{Float64}}()
for k in keys(Rs)
	Qs[k] = quantile(Rs[k],[0.,.1,.25,.5,.75,.9,1.])
end

# #=
figure(journal*" 1")
PyPlot.fill([mi;mi[end:-1:1]],[[Qs[k][2] for k in mi];[Qs[k][6] for k in mi[end:-1:1]]],color=col,alpha=.3)
PyPlot.fill([mi;mi[end:-1:1]],[[Qs[k][3] for k in mi];[Qs[k][5] for k in mi[end:-1:1]]],color=col,alpha=.6)
PyPlot.plot(mi,[Qs[k][4] for k in mi],color=col,linewidth=2)
PyPlot.plot(mi,[Qs[k][1] for k in mi],"--",color=col)
PyPlot.plot(mi,[Qs[k][7] for k in mi],"--",color=col)
xlabel("n")
ylabel("ratio (1999-2008)")
# =#

# #=
figure(journal*" 2")
for k in mi
	PyPlot.plot([k,k],[Qs[k][2],Qs[k][6]],color=col,linewidth=1)
	PyPlot.plot([k,k],[Qs[k][3],Qs[k][5]],color=col,linewidth=5)
end
PyPlot.plot(mi,[Qs[k][4] for k in mi],"o",color=col)
PyPlot.plot(mi,[Qs[k][1] for k in mi],"x",color=col)
PyPlot.plot(mi,[Qs[k][7] for k in mi],"x",color=col)
xlabel("n")
ylabel("ratio (1999-2008)")
# =#

# #=
figure(journal*" 3")
for k in mi
	PyPlot.plot([k,k],[Qs[k][3],Qs[k][5]],color=col,linewidth=lw)
end
PyPlot.plot(mi,[Qs[k][4] for k in mi],"o",color=col)
PyPlot.plot(mi,[Qs[k][1] for k in mi],"x",color=col)
PyPlot.plot(mi,[Qs[k][7] for k in mi],"x",color=col)
xlabel("n")
ylabel("ratio (1999-2008)")
# =#



