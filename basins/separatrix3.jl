using PyPlot, LinearAlgebra

n = 3
res = 200

A = ones(n,n) - diagm(0 => ones(n))

θs = LinRange(-π,π,res)
θ∞ = 0.

V = zeros(res,res)
Λ2 = zeros(res,res)
Λ3 = zeros(res,res)
X = zeros(res,res)

for i in 1:res
	for j in 1:res
		θ = [θs[i],θs[j],θ∞]
		
		V[i,j] = 3 - cos(θ[1] - θ[2]) - cos(θ[1] - θ[3]) - cos(θ[2] - θ[3])
		
		preJ = A.*cos.(θ*ones(1,n) - ones(n)*θ')
		J = preJ - diagm(0 => vec(sum(preJ,dims=1)))
		λs = eigvals(J)
		m,i1 = findmin(abs.(λs))
		m,i2 = findmax([λs[1:i1-1];-Inf;λs[i1+1:n]])
		j1 = min(i1,i2)
		j2 = max(i1,i2)
		Λ2[i,j] = m
		
		m,i3 = findmax([λs[1:j1-1];-Inf;λs[j1+1:j2-1];-Inf;λs[j2+1:n]])
		Λ3[i,j] = m

		X[i,j] = cos(θ[1]-θ[2])*cos(θ[2]-θ[3]) + cos(θ[2]-θ[3])*cos(θ[3]-θ[1]) + cos(θ[3]-θ[1])*cos(θ[1]-θ[2])
	end
end

#subplot(2,2,1)
subplot(1,4,1)
imshow(V,origin="lower",extent=(-π,π,-π,π))
colorbar(label="V")
ylabel("θ2")
xlabel("θ1")

#subplot(2,2,2)
subplot(1,4,2)
#contour(abs.(Λ2),50)
imshow(log.(abs.(Λ2)),origin="lower",extent=(-π,π,-π,π))
colorbar(label="log|λ2|")
xlabel("θ1")

#subplot(2,2,4)
subplot(1,4,3)
#contour(abs.(Λ3),50)
imshow(log.(abs.(Λ3)),origin="lower",extent=(-π,π,-π,π))
colorbar(label="log|λ3|")
xlabel("θ1")

#subplot(2,2,3)
subplot(1,4,4)
imshow(log.(abs.(X)),origin="lower",extent=(-π,π,-π,π))
colorbar(label="log|span. trees condition|")
xlabel("θ1")


