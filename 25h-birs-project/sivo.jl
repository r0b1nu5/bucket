using PyPlot, LinearAlgebra

function sivo(s0,x0,v0,o0,A,B,γ,ρ,h=.01,max_iter=1000,tol=1e-5)
    S = copy(s0)
    X = copy(x0)
    V = copy(v0)
    O = copy(o0)
    s = copy(s0)
    x = copy(x0)
    v = copy(v0)
    o = copy(o0)

    n = length(A[:,1])

    a = abs.(A)*ones(n)
    AA = A - diagm(0 => a)

    iter = 0
    err = 1000

    while iter < max_iter && err > tol
        iter += 1 
#        @info "iter = $iter, err = $err"

        k1 = f_sivo(s,x,v,o,AA,B,γ,ρ)
        k2 = f_sivo(s+k1[1]*h/2,x+k1[2]*h/2,v+k1[3]*h/2,o+k1[4]*h/2,AA,B,γ,ρ)
        k3 = f_sivo(s+k2[1]*h/2,x+k2[2]*h/2,v+k2[3]*h/2,o+k2[4]*h/2,AA,B,γ,ρ)
        k4 = f_sivo(s+k3[1]*h,x+k3[2]*h,v+k3[3]*h,o+k3[4]*h,AA,B,γ,ρ)

        ds = (k1[1]+2*k2[1]+2*k3[1]+k4[1])/6
        s += ds*h
        dx = (k1[2]+2*k2[2]+2*k3[2]+k4[2])/6
        x += dx*h
        dv = (k1[3]+2*k2[3]+2*k3[3]+k4[3])/6
        v += dv*h
        d_o = (k1[4]+2*k2[4]+2*k3[4]+k4[4])/6
        o += d_o*h

        S = [S s]
        X = [X x]
        V = [V v]
        O = [O o]
    end

    return S,X,V,O
end


function f_sivo(s,x,v,o,A,B,γ,ρ)
    κ = (1 .- 2*o)./(1 .+ 2*o)
    fs = γ*x - s.*(B*x + ρ*(1 .- v./κ))
    fx = s.*(B*x) .- γ*x
    fv = ρ*(1 .- v./κ).*s
    fo = (x - o .- .5) + A*o

    return fs,fx,fv,fo
end


function get_Jxx(s,x,v,o,A,B,γ)
    n = length(x)

    κ = (1 .- 2*o)./(1 .+ 2*o)

    return ((diagm(0 => ones(n)) - diagm(0 => κ) - diagm(0 => x))*B - diagm(0 => (B*x)) - γ*diagm(0 => ones(n)))
end






