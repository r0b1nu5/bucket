using PyPlot, DelimitedFiles

ids = [1,2,3,4]
n = 83

res1 = 50
res2 = 200

Qf = Matrix{Float64}(undef,0,res2)
Cst = Matrix{Float64}(undef,0,res2)
Λs = Vector{Matrix{Float64}}()
for i in 1:n-1
	push!(Λs,Matrix{Float64}(undef,0,res2))
end

for id in ids
	global Qf,Cst,Λs

# #=
	Qf = [Qf;readdlm("data/Qf_$(id).csv",',')]
	Cst = [Cst;readdlm("data/Cs_$(id).csv",',')]
	for i in 1:n-1
		Λs[i] = [Λs[i];readdlm("data/L$(i)_$(id).csv",',')]
	end
	global αmin = -π
	global αmax = π
# =#
 #=
	Qf = [Qf;readdlm("data/Qfbis_$(id).csv",',')]
	Cst = [Cst;readdlm("data/Csbis_$(id).csv",',')]
	for i in 1:n-1
		Λs[i] = [Λs[i];readdlm("data/L$(i)bis_$(id).csv",',')]
	end
	global αmin = -1.
	global αmax = 1.
# =#
end

figure("Qf")
imshow(Qf,origin="lower",extent=(αmin,αmax,αmin,αmax))
colorbar()
xlabel("α2")
ylabel("α1")

figure("Cst")
imshow(log.(abs.(Cst)),origin="lower",extent=(αmin,αmax,αmin,αmax))
colorbar()
xlabel("α2")
ylabel("α1")

figure("λ2")
imshow(Λs[1],origin="lower",extent=(αmin,αmax,αmin,αmax))
colorbar()
xlabel("α2")
ylabel("α1")

figure("λ3")
imshow(Λs[2],origin="lower",extent=(αmin,αmax,αmin,αmax))
colorbar()
xlabel("α2")
ylabel("α1")

figure("λn")
imshow(Λs[end],origin="lower",extent=(αmin,αmax,αmin,αmax))
colorbar()
xlabel("α2")
ylabel("α1")


