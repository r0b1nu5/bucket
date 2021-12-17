include("acyclic_algorithm_repo.jl")
include("cyclic_iterations_repo.jl")
include("toolbox_repo.jl")


@info "==========================================================="
@info "Running the acyclic algorithm on a spanning tree of RTS-96."
@info "==========================================================="

Lst = readdlm("ntw_data/rts96_spantree_L.csv",',')
Yst = readdlm("ntw_data/rts96_spantree_Lg.csv",',') + im*readdlm("ntw_data/rts96_spantree_Lb.csv",',')
id_st = Int64.(vec(readdlm("ntw_data/rts96_ids_spantree.csv",',')))

ω = vec(readdlm("ntw_data/rts96_om1.csv",','))

h,hi,γ,Bst = load_ksakaguchi(Yst)

exists,θ,ff,φ = run_acyclic_algorithm(Bst,ω,h,hi,γ)

if exists
	dθ = cohesiveness_inc(θ,Bst)
	@info "The unique solution has been found. With maximal angular difference Δ = $(dθ)."
else
	@info "No solution as been found."
end




@info "================================================================"
@info "Running the cyclic iterations on RTS-96, in three winding cells."
@info "================================================================"

L = readdlm("ntw_data/rts96_w2_L.csv",',')
B,w,Bt = L2B(L)
ϕs = vec(readdlm("ntw_data/rts96_w2_phi.csv",','))
Bstd = pinv(Bst)

include("ntw_data/rts96_w_cycles.jl")

u1 = Int64.(vec(readdlm("ntw_data/rts96_w_u1.csv",',')))
u2 = Int64.(vec(readdlm("ntw_data/rts96_w_u2.csv",',')))
u3 = Int64.(vec(readdlm("ntw_data/rts96_w_u3.csv",',')))

ω = vec(readdlm("ntw_data/rts96_w2_om.csv",','))

h,hi,γ = load_ksakaguchi([w;w],[ϕs;ϕs])
Bb = [B -B]
Bout = Bb.*(Bb .> 1e-2)
n,m = size(B)

δ = .01

s = 1.

max_iter = 1000
tol = 1e-5

Δ0 = [((γ[i][2]-γ[i][1])*rand() + γ[i][1]) for i in 1:m]


@info "First case:"
Δ1,Δ1s = iterations(Δ0,B,C,u1,ω,h,γ,δ,s,max_iter,tol)
θ1 = Bstd'*Δ1[id_st]
q1 = winding(θ1,Σ)
f1 = [[h[i](Δ1[i]) for i in 1:m]; [h[i+m](-Δ1[i]) for i in 1:m]]
res1 = ω - Bout*f1

if q1 == u1
	@info "Winding vector matches."
else
	@info "Winding vector does not match."
end
if maximum(res1) - minimum(res1) > 1e-2
	@info "Did not find a solution."
else
	@info "Found a solution."
end
@info "-------------------------------------------"


@info "Second case:"
Δ2,Δ2s = iterations(Δ0,B,C,u2,ω,h,γ,δ,s,max_iter,tol)
θ2 = Bstd'*Δ2[id_st]
q2 = winding(θ2,Σ)
f2 = [[h[i](Δ2[i]) for i in 1:m]; [h[i+m](-Δ2[i]) for i in 1:m]]
res2 = ω - Bout*f2

if q2 == u2
	@info "Winding vector matches."
else
	@info "Winding vector does not match."
end
if maximum(res2) - minimum(res2) > 1e-2
	@info "Did not find a solution."
else
	@info "Found a solution."
end
@info "-------------------------------------------"

@info "Third case:"
Δ3,Δ3s = iterations(Δ0,B,C,u3,ω,h,γ,δ,s,max_iter,tol)
θ3 = Bstd'*Δ3[id_st]
q3 = winding(θ3,Σ)
f3 = [[h[i](Δ3[i]) for i in 1:m]; [h[i+m](-Δ3[i]) for i in 1:m]]
res3 = ω - Bout*f3

if q3 == u3
	@info "Winding vector matches."
else
	@info "Winding vector does not match."
end
if maximum(res3) - minimum(res3) > 1e-2
	@info "Did not find a solution."
else
	@info "Found a solution."
end
@info "-------------------------------------------"




