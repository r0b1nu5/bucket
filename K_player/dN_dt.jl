using PyPlot, Logging, DelimitedFiles
# #=
include("journals.jl")

#ref_year = "00"
ref_year = "80"

journal = "science"
yi = 99
ny = 10
hor = 3

nums = Array{Array{Float64,1},1}()
push!(nums,collect(Float64,readdlm("data/time_split/"*journal*"_"*ref_year*"_$yi.txt",'\t')[2:end-2,2]))
for i in 1:ny
	push!(nums,collect(Float64,readdlm("data/time_split/"*journal*"_"*ref_year*"_$((yi+i)%100).txt",'\t')[2:end-2,2]))
end

ma = maximum(nums[end])

Ns = Array{Array{Float64,1},1}()
N = Array{Float64,1}()
for i in 1:Int(ma)
	push!(N,sum(nums[1] .== i))
end
push!(Ns,N)

dNs = Array{Array{Float64,1},1}()
mtot = Array{Float64,1}()

for num in nums[2:end]
	N = Array{Float64,1}()
	for i in 1:Int(ma)
		push!(N,sum(num .== i))
	end
	push!(Ns,N)
	push!(dNs,Ns[end]-Ns[end-1])
	push!(mtot,sum(dNs)[1])
end

ms = Array{Array{Float64,1},1}()

aut = readdlm("data/time_split/"*journal*"_"*ref_year*"_$yi.txt",'\t')[2:end-2,1]


for i in 1:(ny-hor+1)
	numi = collect(Float64,readdlm("data/time_split/"*journal*"_"*ref_year*"_$((yi+i-1)%100).txt",'\t')[2:end-2,2])
	auti = readdlm("data/time_split/"*journal*"_"*ref_year*"_$((yi+i-1)%100).txt",'\t')[2:end-2,1]
	dnumi = collect(Float64,readdlm("data/time_split/"*journal*"_$((yi+i)%100).txt",'\t')[2:end-2,2])
	dauti = readdlm("data/time_split/"*journal*"_$((yi+i)%100).txt",'\t')[2:end-2,1]
	
	for k in 2:hor
		x = collect(Float64,readdlm("data/time_split/"*journal*"_$((yi+i+k-1)%100).txt",'\t')[2:end-2,2])
		y = readdlm("data/time_split/"*journal*"_$((yi+i+k-1)%100).txt",'\t')[2:end-2,1]
		for l in 1:length(x)
			idx = (dauti .== y[l]).*(1:length(dauti))
			os = intersect(idx,(1:length(dauti)))

			if length(os) > 0
				dnumi[os] .+= x[l]/length(os)
			else
				push!(dauti,y[l])
				push!(dnumi,x[l])
			end
		end
	end
	dnumi /= hor

	m = zeros(Int(10*ma))
	
	for j in 1:length(dnumi)
		if j%1000 == 0
		@info(" $i/$ny : $j/$(length(dnumi))")
		end
		
		nu = dnumi[j]
		au = dauti[j]
		
		idx = (auti .== au).*(1:length(auti))
		ks = intersect(idx,(1:length(auti)))
		
		if length(ks) > 0
			mean_n = round(nu/length(ks))
			for k in ks
				nold = numi[k]
				nnew = nold + mean_n
				m[Int(nold)] += mean_n
			end
		end
	end
	
	push!(ms,m)
end

# =#

x = Array{Float64,1}()
y = Array{Float64,1}()

for i in 1:(ny-hor+1)
	for j in 1:length(Ns[i])
		if Ns[i][j] > .1 && ms[i][j] > 0
			push!(x,j)
 			push!(y,ms[i][j]./Ns[i][j])
		end
	end
end

mx = sum(x)/length(x)
my = sum(y)/length(y)
mxy = sum(x.*y)/length(x)
sx = sqrt(sum((x.-mx).^2)/length(x))
sy = sqrt(sum((y.-my).^2)/length(y))
r = (mxy - mx*my)/(sx*sy)

figure()
PyPlot.plot(x,y,"o",color=journals_colors[journal][1])
xlabel("Number of articles (k)")
ylabel("Number of new articles per author with k articles")
title(journal*": Correlation = $(r-r%.0001)")


