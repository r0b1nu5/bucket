using DelimitedFiles,LinearAlgebra,PyPlot

include("res_dist.jl")

Asp = readdlm("uk_adj_mat.csv",',') .+ 1.0
n = Int(maximum(Asp))
m = Int(size(Asp)[1]/2)
mm = copy(m)

A = zeros(Int,n,n)
for i in 1:330
	A[Int(Asp[i,1]),Int(Asp[i,2])] = 1.0
end

L = Array{Float64,2}(diagm(0 => vec(sum(A, dims = 2))) - A)

Om = res_dist(L)

line_list = Array{Array{Int,1},1}()
dist = Array{Float64,1}()
for i in 1:m
	push!(line_list,Array{Int,1}(vec([sort(Asp[2*i-1,[1,2]]);i])))
	push!(dist,Om[line_list[end][1],line_list[end][2]])
end

rank1 = Array{Int64,}(vec(sortslices([dist 1:m],dims=1)[:,2]))

L2 = copy(L)
L3 = copy(L)
iter = 0
i2 = 0
i3 = 0
Om2 = copy(Om)
Om3 = copy(Om)
ll2 = copy(line_list)
ll3 = copy(line_list)

rank2 = [Int(rank1[1]),]
rank3 = [Int(rank1[1]),]
rankt = copy(rank1)

while m > n-1
	global n,m,L2,L3,iter,i2,i3,Om2,Om3,line_list,ll2,ll3,rank1,rank2,rank3,rankt
	iter += 1
	m -= 1
	
	@info("iter = $iter, n = $n, m = $m")
	
	l2 = line_list[rank1[iter+i2]][1:2]
	while abs(-L2[l2[1],l2[2]] - Om2[l2[1],l2[2]]) < 1e-8
		i2 += 1
		l2 = line_list[rank1[iter+i2]][1:2]
	end
	b2 = -L2[l2[1],l2[2]]
	L2[l2,l2] -= b2*[1 -1;-1 1]
	
	i3 = 0
	l3 = line_list[rank3[end]][1:2]
#	while abs(-L3[l3[1],l3[2]] - Om3[l3[1],l3[2]]) < 1e-8
#		i3 += 1
#		l3 = line_list[rankt[1+i3]]
#	end
	b3 = -L3[l3[1],l3[2]]
	L3[l3,l3] -= b3*[1 -1;-1 1]
	
	ll2 = Array{Array{Int,1},1}()
	ll3 = Array{Array{Int,1},1}()
	for i in 1:n-1
		for j in i+1:n
			if (abs(L2[i,j]) > 1e-6) && (abs(-L2[i,j]-Om2[i,j]) > 1e-8)
				push!(ll2,[i,j])
			end
			if (abs(L3[i,j]) > 1e-6) && (abs(-L3[i,j]-Om3[i,j]) > 1e-8)
				push!(ll3,[i,j])
			end
		end
	end
	
	Om2 = res_dist(L2)
	Om3 = res_dist(L3)
	
	dist2 = Array{Float64,1}()
	dist3 = Array{Float64,1}()
	for i in 1:length(ll2)
		push!(dist2,Om2[ll2[i][1],ll2[i][2]])
	end
	for i in 1:length(ll3)
		push!(dist3,Om3[ll3[i][1],ll3[i][2]])
	end
	
	w2 = ll2[Int(sortslices([dist2 1:length(dist2)],dims=1)[1,2])]
	w3 = ll3[Int(sortslices([dist3 1:length(dist3)],dims=1)[1,2])]
	
	for l in line_list
		if w2 == l[1:2]
			push!(rank2,l[3])
		end
		if w3 == l[1:2]
			push!(rank3,l[3])
		end
	end
	
#	rankt = Array{Int64,1}(vec(sortslices([dist3 1:m],dims=1)[:,2]))
end

figure(1,(10,7))
PyPlot.subplot(1,3,1)
PyPlot.plot(rank1[1:mm-n+2],rank2,"o")

PyPlot.subplot(1,3,2)
PyPlot.plot(rank1[1:mm-n+2],rank3,"o")

PyPlot.subplot(1,3,3)
PyPlot.plot(rank2,rank3,"o")


