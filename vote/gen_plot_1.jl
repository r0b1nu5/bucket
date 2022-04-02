using PyPlot, NPZ

x = npzread("data/NUM_VOTES.npy")

N = 10
ϵ1 = .25
ϵ2 = 1.
ne = 20
ϵs = LinRange(ϵ1,ϵ2,ne)
de = (ϵs[2]-ϵs[1])/3
niter = 200

win = zeros(Int64,niter,ne)
for i in 1:niter
	for j in 1:ne
		a,b = findmax(x[i,j,:])
		win[i,j] = b
	end
end

p = zeros(N,ne)
for i in 1:niter
	for j in 1:ne
		p[win[i,j],j] += 1/niter
	end
end

for i in 1:N
	subplot(1,3,1)
	PyPlot.plot(ϵs,100*p[i,:],label="party $i")
	s = 100*vec(sum(p[1:N+1-i,:],dims=1))
	subplot(1,3,2)
	PyPlot.fill([ϵs;ϵ2;ϵ1],[s;0.;0.])
	subplot(1,3,3)
	for j in 1:ne
		ϵ = ϵs[j]
		PyPlot.fill([ϵ-de,ϵ+de,ϵ+de,ϵ-de],[0.,0.,s[j],s[j]],color="C$(mod(i-1,10))")
	end
end

subplot(1,3,1)
legend()
xlabel("ϵ")
ylabel("% of win")

subplot(1,3,2)
xlabel("ϵ")

subplot(1,3,3)
xlabel("ϵ")





