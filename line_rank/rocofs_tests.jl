using DelimitedFiles,Statistics,PyPlot,SparseArrays,LinearAlgebra

include("kuramoto.jl")
include("res_dist.jl")

ntw = "ieee118"
Lsp = readdlm(ntw*"_data/"*ntw*"_lap_mat.csv",',')
L = sparse(Lsp[:,1],Lsp[:,2],Lsp[:,3])

Om = res_dist(L)

n = size(L)[1]
m = round(Int,size(Lsp)[1]/2)

P0 = .001
ratedP = vec(readdlm(ntw*"_data/P_"*ntw*".csv",','))
P = P0*ratedP
P .-= mean(P)

g = .5
# For inertias, see Laurent's Plos One, Appendix 2, Eq. (S2)
H = 2*pi .+ 2*rand(n) .- 1.
omega0 = 50*2*pi
rP = 2*max.(abs.(ratedP),maximum(abs.(ratedP))/100*ones(n))
M = 2*H.*rP./omega0
D = g*M
LD = diagm(0 => D)^(-1/2)*L*diagm(0 => D)^(-1/2)

# Lines to cut: largest bij/mi, median bij/mi, smallest bij/mi,...


lines2cut = [(5,4),(37,34),(12,16),(96,94),(89,85),(69,49)]


# Distribution P: large variance on slow modes, large variance on fast modes,...

l4
  92.0
 103.0
  83.0
 104.0
  91.0
  89.0
 106.0
  90.0
 105.0
 107.0
  88.0
  84.0
  85.0
 108.0
 109.0
 110.0
 112.0
 111.0
  86.0
  87.0

l3
  60.0
 109.0
  47.0
  49.0
 110.0
  45.0
  59.0
  48.0
  46.0
 112.0
 111.0
  50.0
  54.0
  55.0
  56.0
  57.0
  51.0
  58.0
  53.0
  52.0

l2
 100.0
  92.0
  88.0
  89.0
 101.0
 102.0
  90.0
  91.0
  86.0
 103.0
  87.0
 104.0
 106.0
 105.0
 107.0
 108.0
 109.0
 110.0
 112.0
 111.0

l118
  78.0
  67.0
  60.0
  79.0
  37.0
  47.0
  70.0
  75.0
  77.0
  61.0
  63.0
  38.0
  80.0
  66.0
  64.0
  69.0
  81.0
  65.0
 116.0
  68.0

l117
   2.0
 113.0
  34.0
  26.0
  37.0
  13.0
  38.0
  17.0
   1.0
  10.0
  12.0
   7.0
  30.0
  11.0
   9.0
   3.0
   6.0
   8.0
   4.0
   5.0

l116
   4.0
   8.0
   5.0
  66.0
  41.0
  17.0
 116.0
  64.0
  65.0
  30.0
  19.0
  43.0
  40.0
  33.0
  39.0
  38.0
  35.0
  36.0
  34.0
  37.0





