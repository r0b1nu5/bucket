using LinearAlgebra, PyPlot

Gbr = [3.8,5.2]
Bbr = [19.1,25.8]
Ybr = Gbr - im*Bbr

E = [1 0;
     -1 1;
     0 -1]

Ybu = E*diagm(0 => Ybr)*E'
Y = Ybu + diagm(0 => rand(Complex{Float64},3)) # add shunt
Yi = inv(Y)

vs = [1.,.98,.97].*exp.(im*[0.,-.1,-.15])
ws = [vs[1]*conj(vs[2]),vs[2]*conj(vs[3])]

Sbar = conj(vs).*Y*vs
P = real(Sbar)
Q = -imag(Sbar)
S = P + im*Q

v0 = Complex{Float64}.(ones(3))
dv0 = v0 - vs

# Step 1
dvi0 = diagm(0 => 1 ./v0)
v1 = Yi*dvi0*S
dv1 = v1 - vs

# Step 2
dvi1 = diagm(0 => 1 ./v1)
v2 = Yi*dvi1*S
dv2 = v2 - vs

# Step 3
dvi2 = diagm(0 => 1 ./v2)
v3 = Yi*dvi2*S
dv3 = v3 - vs

# Step 4
dvi3 = diagm(0 => 1 ./v3)
v4 = Yi*dvi3*S
dv4 = v4 - vs

# Show err
dv = [dv0 dv1 dv2 dv3 dv4]
figure()
subplot(1,3,1)
PyPlot.plot(0:4,abs.(dv[2,:]),"o")
PyPlot.plot(0:4,abs.(dv[3,:]),"o")
subplot(1,3,2)
PyPlot.plot(0:4,abs.(real.(dv[2,:])),"o")
PyPlot.plot(0:4,abs.(real.(dv[3,:])),"o")
subplot(1,3,3)
PyPlot.plot(0:4,abs.(imag.(dv[2,:])),"o")
PyPlot.plot(0:4,abs.(imag.(dv[3,:])),"o")



#=
# Initialize 
w0 = Complex{Float64}.([1.,1.]) # [v1v2,v2v3]
dw0 = w0 - ws

# Step 1
w1 = copy(w0)
w1[2] = (conj(S[3]) - Y[3,3]*(abs(w0[2])/abs(w0[1]))^2)/Y[2,3]
w1[1] = (conj(S[2]) - Y[2,2]*(abs(w0[1]))^2)/Y[1,2] - Y[2,3]*w1[2]
dw1 = w1 - ws

# Step 2
w2 = copy(w1)
w2[2] = (conj(S[3]) - Y[3,3]*(abs(w1[2])/abs(w1[1]))^2)/Y[2,3]
w2[1] = (conj(S[2]) - Y[2,2]*(abs(w1[1]))^2)/Y[1,2] - Y[2,3]*w2[2]
dw2 = w2 - ws

# Step 3
w3 = copy(w2)
w3[2] = (conj(S[3]) - Y[3,3]*(abs(w2[2])/abs(w2[1]))^2)/Y[2,3]
w3[1] = (conj(S[2]) - Y[2,2]*(abs(w2[1]))^2)/Y[1,2] - Y[2,3]*w3[2]
dw3 = w3 - ws

# Step 4
w4 = copy(w3)
w4[2] = (conj(S[3]) - Y[3,3]*(abs(w3[2])/abs(w3[1]))^2)/Y[2,3]
w4[1] = (conj(S[2]) - Y[2,2]*(abs(w3[1]))^2)/Y[1,2] - Y[2,3]*w4[2]
dw4 = w4 - ws

# Show error
dw = [dw0 dw1 dw2 dw3 dw4]
figure()
subplot(1,3,1)
PyPlot.plot(0:4,abs.(dw[1,:]),"o")
PyPlot.plot(0:4,abs.(dw[2,:]),"o")
subplot(1,3,2)
PyPlot.plot(0:4,abs.(real.(dw[1,:])),"o")
PyPlot.plot(0:4,abs.(real.(dw[2,:])),"o")
subplot(1,3,3)
PyPlot.plot(0:4,abs.(imag.(dw[1,:])),"o")
PyPlot.plot(0:4,abs.(imag.(dw[2,:])),"o")




# Step 1
v1 = copy(v0)
v1[2] = (S[2] - conj(Y[2,2])*abs(v0[2])^2)/Y[1,2] - (conj(S[3]) - Y[3,3]*abs(v0[3])^2)*conj(Y[2,3])/(conj(Y[1,2])*Y[2,3])
v1[3] = ((S[3] - conj(Y[3,3])*abs(v0[3])^2)/conj(Y[2,3]))/conj(v1[2])
dv1 = v1 - vs

# Step 2
v2 = copy(v1)
v2[2] = (S[2] - conj(Y[2,2])*abs(v1[2])^2)/Y[1,2] - (conj(S[3]) - Y[3,3]*abs(v1[3])^2)*conj(Y[2,3])/(conj(Y[1,2])*Y[2,3])
v2[3] = ((S[3] - conj(Y[3,3])*abs(v1[3])^2)/conj(Y[2,3]))/conj(v2[2])
dv2 = v2 - vs

# Step 3
v3 = copy(v2)
v3[2] = (S[2] - conj(Y[2,2])*abs(v2[2])^2)/Y[1,2] - (conj(S[3]) - Y[3,3]*abs(v2[3])^2)*conj(Y[2,3])/(conj(Y[1,2])*Y[2,3])
v3[3] = ((S[3] - conj(Y[3,3])*abs(v2[3])^2)/conj(Y[2,3]))/conj(v3[2])
dv3 = v3 - vs

# Step 4
v4 = copy(v3)
v4[2] = (S[2] - conj(Y[2,2])*abs(v3[2])^2)/Y[1,2] - (conj(S[3]) - Y[3,3]*abs(v3[3])^2)*conj(Y[2,3])/(conj(Y[1,2])*Y[2,3])
v4[3] = ((S[3] - conj(Y[3,3])*abs(v3[3])^2)/conj(Y[2,3]))/conj(v4[2])
dv4 = v4 - vs


# Show error
dv = [dv0 dv1 dv2 dv3 dv4]
figure()
subplot(1,3,1)
PyPlot.plot(0:4,abs.(dv[2,:]),"o")
PyPlot.plot(0:4,abs.(dv[3,:]),"o")
subplot(1,3,2)
PyPlot.plot(0:4,abs.(real.(dv[2,:])),"o")
PyPlot.plot(0:4,abs.(real.(dv[3,:])),"o")
subplot(1,3,3)
PyPlot.plot(0:4,abs.(imag.(dv[2,:])),"o")
PyPlot.plot(0:4,abs.(imag.(dv[3,:])),"o")

# =#




