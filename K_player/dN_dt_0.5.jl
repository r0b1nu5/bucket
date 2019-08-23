using PyPlot, Logging
@Logging.configure(level=INFO)
include("journals.jl")

# #=
journals = ["prl","science","lancet"]
journal = journals[1]
spans = ["00","80"]
span = spans[1]
# =#

yi = 99
ny = 10
hor = 1

# each column of aut_year is a year, and each line is an author, each element is the number of articles published by author i in year j.
aut_fin = sort(readdlm("data/time_split/"*journal*"_"*span*"_$((yi+ny)%100).txt",'\t')[2:end-2,1])
aut_year = zeros(length(aut_fin),ny+1)
@info(" Year: $yi")
an = sortrows(readdlm("data/time_split/"*journal*"_"*span*"_$yi.txt",'\t')[2:end-2,1:2])
j = 1
k = 1
while j <= length(aut_fin) && k <= size(an)[1]
	if an[k,1] == aut_fin[j]
		aut_year[j,1] = an[k,2]
		k += 1
		j += 1
	elseif an[k,1] < aut_fin[j]
		k += 1
	else
		j += 1
	end
end

for i in 1:ny
	@info(" Year: $((yi+i)%100)")
#	an = sortrows(readdlm("data/time_split/"*journal*"_$((yi+i)%100).txt",'\t')[2:end-2,1:2])	
	an = sortrows(readdlm("data/time_split/"*journal*"_"*span*"_$((yi+i)%100).txt",'\t')[2:end-2,1:2])	
	j = 1
	k = 1
	while j <= length(aut_fin) && k <= size(an)[1]
		if an[k,1] == aut_fin[j]
			aut_year[j,i+1] = an[k,2]
			k += 1
			j += 1
		elseif an[k,1] < aut_fin[j]
			k += 1
		else
			j += 1
		end	
	end
end


x = Array{Float64,1}()
y = Array{Float64,1}()

for i in 1:ny-hor
	@info(" Year: $i/$(ny-hor)")
#	num = collect(sum(aut_year[:,1:i],2))
	num = collect(aut_year[:,i])
	for j in 1:maximum(num)
		idx = (num .== j).*(1:length(num))
		ids = intersect(idx,1:length(num))
#		nnew_papers = sum(aut_year[ids,i+1:i+hor])/hor
		nnew_papers = sum(aut_year[ids,i+hor]-aut_year[ids,i])/hor
		if length(ids) > 0 && nnew_papers > 0
#		if length(ids) > 0
			new_p_aut = nnew_papers/length(ids)
			push!(x,j)
			push!(y,new_p_aut)
		end
	end
end



sx = sqrt(mean((x-mean(x)).^2))
sy = sqrt(mean((y-mean(y)).^2))
mx = mean(x)
my = mean(y)
mxy = mean(x.*y)
r = (mxy - mx*my)/(sx*sy)

figure(journal)
PyPlot.plot(x,y,"o",color=journals_colors[journal][1],markersize=8,label="horizon = $hor, r = $(round(r,4))")
xlabel("Number of articles (k)")
ylabel("Number of new articles per author with k articles")
legend()
title(journals_code[journal])

