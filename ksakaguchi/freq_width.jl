include("ksakaguchi.jl")
include("tools.jl")
include("plots.jl")
include("cyqle.jl")

n = 10

L = cyqle(n)

a1 = .1
a0 = 0.

x = ksakaguchi(L,zeros(n),zeros(n),a1,true,false,.01,1e-7)
th0 = x[1][:,end]
q0 = winding(th0,Array(1:n))
x = ksakaguchi(L,zeros(n),Array((2pi/n)*(1:n)),a1,true,false,.01,1e-7)
th1 = x[1][:,end]
q1 = winding(th1,Array(1:n))
x = ksakaguchi(L,zeros(n),Array((2pi/n)*(n:-1:1)),a1,true,false,.01,1e-7)
th2 = x[1][:,end]
q2 = winding(th2,Array(1:n))

 #=
om = rand(n)
om .-= mean(om)
om /= norm(om)
# =# 

b = 0.
db = 1.
th = copy(th0)

while db > .0005
	global b,db,th,om,L,a1,x
	@info "b = $b, db = $db"
	q = q0
	it = 0
	while q == q0 && it < 100000
		b += db
		x = ksakaguchi(L,b*om,th,a1,true,false,.01,1e-6)
		q = winding(x[1][:,end],Array(1:n))
		it = x[4]
		@info "q = $q, it = $it"
	end
	th = x[1][:,1]
	b -= db
	db /= 10
end

@info "$b"

		



