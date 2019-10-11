using DynamicalBilliards, PyPlot

function animate_stadium()
	bd = billiard_stadium()
	
	N = 20

	cs = [(i/N,0,1-i/N,.5) for i in 1:N]
	ps = [Particle(1,.6+.0005*i,0) for i in 1:N]

	animate_evolution(ps,bd,20.;colors = cs, tailtime = 1.5)

end


function animate_sinai()
	bd = billiard_sinai()
	
	N = 20

	cs = [(i/N,0,1-i/N,.5) for i in 1:N]
	ps = [Particle(0.01,.6+.0005*i,0) for i in 1:N]

	animate_evolution(ps,bd,10.;colors = cs, tailtime = 1.5)

end

function animate_circle()
	ob = Obstacle{Float64}[]
	push!(ob,Semicircle([0.,0.],2.,[1.,0],"Left"))
	push!(ob,Semicircle([0.,0.],2.,[-1.,0],"Right"))

	bd = Billiard(ob)

	N = 10
	cs = [(i/N,0,1-i/N,.5) for i in 1:N]
	ps = [Particle(1.1*cos(i/N),1.1*sin(i/N),i/N+pi/2) for i in 1:N]

	animate_evolution(ps,bd,30.;colors=cs,tailtime = 30.)
end




