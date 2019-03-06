using PyPlot,LinearAlgebra

n = 10
th = pi/2*rand(n)
a = pi/5

D = Array{Float64,1}()
for i in 1:n
	push!(D,sum(cos.(th.-th[i])))
end

ls = LinRange(6.5,10.5,10000)

khi = Array{Float64,1}()
khi0 = Array{Float64,1}()
khin = Array{Float64,1}()

for t in 1:10000
	global khi,khi0,khin
	temp = prod(D.-ls[t])
	tempn = 1.
	push!(khi0,temp)
	for i in 1:n
		temp2 = cos(a)
		for j in 1:n
			if i!=j
				temp *= (D[j]-ls[t])
			end
		end
		temp += temp2
		tempn += cos(a)/(D[i]-ls[t])
	end
	for i in 1:n-1
		for j in i+1:n
			temp2 = sin(th[i]-th[j])^2
			for k in 1:n
				if k!=i && k!=j
					temp2 *= (D[k]-ls[t])
				end
			end
			temp += temp2
			tempn += sin(th[i]-th[j])^2/((D[i]-ls[t])*(D[j]-ls[t]))
		end
	end
	push!(khi,temp)
	push!(khin,tempn)
end

PyPlot.plot(ls,khi0,label="X0")
PyPlot.plot(ls,khi,label="X")
PyPlot.plot(ls,khin,label="X/X0")
legend()


