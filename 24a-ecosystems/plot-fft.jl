using DelimitedFiles

x = readdlm("data/FFT.dat")

F = zeros(49,199)
f = zeros(49,199)
σs = LinRange(2,4,199)
ks = 2:50
for l in 1:size(x)[1]
	j = ceil(Int64,l/49)
	i = mod(l-1,49) + 1

	if !isnan(x[l,3])
		F[i,j] = log.(x[l,3] + 1e-6)
	end
	f[i,j] = x[l,3]
end

figure("2d")
PyPlot.contourf(σs[80:149],ks,F[:,80:149],100)

figure("3d")
PyPlot.surf(σs[80:149],ks,F[:,80:149])


