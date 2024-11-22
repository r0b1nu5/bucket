using Random, Dates, DelimitedFiles

include("hyper_inf.jl")
include("hyper_kuramoto.jl")
include("hyper_ktanh.jl")
include("gen_rand_hyperg.jl")
include("tools_hyper.jl")

@info "############# START: $(now())"

n = 4; T = 200
n = 8; T = 600
n = 16; T = 1000
n = 32; T = 1200
n = 64; T = 1500
n = 128; T = 2000
n = 256; T = 2000

save = true
A2 = zeros(n,n)
A3 = zeros(n,n,n)
A2l = zeros(0,3)
A3l = zeros(0,4)


ntw = "Simplicial-ER-py"
run = "999"
p1 = .001
p2 = .01
# =#
 #=
ntw = "Hyper-ER-py"
run = "002"
p1 = .001
p2 = .01
# =#
 #=
ntw = "Simplicial-ER-corr"
p1 = n > 30 ? .01 : .05
p2 = .4
# =#

el = readdlm("data/edgelist-n$n-"*run*".csv",',')
for l in 1:size(el)[1]
	global i,j,k,A2,A2l,A3,A3l
	i,j,k = el[l,:]
	if k == ""
		i,j = sort([i,j] .+ 1)
		A2[i,j] = A2[j,i] = 1.
		A2l = vcat(A2l,[i j 1.;j i 1.])
	else
		i,j,k = sort([i,j,k] .+ 1)
		A3[i,j,k] = A3[i,k,j] = A3[j,i,k] = A3[j,k,i] = A3[k,i,j] = A3[k,j,i] = 1.
		A3l = vcat(A3l,[i j k 1.;j i k 1.;k i j 1.])
	end
end


adj = cat_As(A2,A3)
amplitude = 1.
ξ0 = 0.0005
X = amplitude*(rand(n,T) .- .5)
Y = f_kuramoto_3rd(X,A2,A3,zeros(n),π/4,π/4) + ξ0*randn(size(X))

ooi = [2,3]
dmax = 2
c = 0
nkeep = min(n^2,10,ceil(Int64,n^2*0.05))

t0 = time()
xxx = hyper_inf_filter(X[:,1:iter],Y[:,1:iter],ooi,dmax,nkeep,1e-1) 
A2us = xxx[1][2]
A3us = xxx[1][3]
t1 = time()

if save
	writedlm("data/kuramoto-"*ntw*"-n$n-timetest-A2.csv",A2,',')
	writedlm("data/kuramoto-"*ntw*"-n$n-timetest-A3.csv",A3,',')
	writedlm("data/kuramoto-"*ntw*"-n$n-timetest-A2this.csv",A2us,',')
	writedlm("data/kuramoto-"*ntw*"-n$n-timetest-A3this.csv",A3us,',')
	writedlm("data/kuramoto-"*ntw*"-n$n-timetest.csv",t1-t0,',')
end

@info "============= WE ARE DONE: $(now()) ================"

tpr,fpr = my_ROC(abs.(A2us),A2l,abs.(A3us),A3l,n)
tpr2,fpr2 = my_ROC(abs.(A2us),A2l,n)
tpr3,fpr3 = my_ROC(abs.(A3us),A3l,n)

@info "############# FINISHED: $(now())"

