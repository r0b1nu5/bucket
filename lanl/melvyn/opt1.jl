

using JuMP, LinearAlgebra, GLPK, Test, DelimitedFiles, Ipopt, PyPlot

function example_basic(; verbose = true)
	model = Model(with_optimizer(Ipopt.Optimizer))
	ns=20
	@variable(model, A[1:ns*2,1:ns*2])
	@variable(model, c[1:2*ns])
	@variable(model, f[1:2*ns])
	@variable(model, phi[1:2*ns])
	for i in ns+1:2*ns
		@constraint(model, c[i]==0)
		@constraint(model, f[i]==0)
		@constraint(model, phi[i]==0)
		@constraint(model, -4.0<=c[i]<=4.0)

	end
	X=readdlm("data/C20_50000_0.2.csv",',')'
	X1=(A*X')'


	println("OK")
	@objective(model, Min, sum(sum((X[i+1,k]-X1[i,k])^2 for k in 1:2*ns) for i in 1:10000))

	for i in 1:ns
		for j in 1:ns  
			if(i!=j)
				 #@constraint(model, A[i,i] == 1.0)
			#else
				 @constraint(model, A[i,j] == 0.0)	
			end
		end
	end
	for i in 1:ns

		for j in ns+1:2*ns   
			 @constraint(model, A[i,j] == 0.0)
		end
	end
	for i in ns+1:2*ns
		for j in ns+1:2*ns   

			if(i!=j)
				 @constraint(model, A[i,j] == 0.0)	
			end
		end
	end
    if verbose
        print(model)
    end

    JuMP.optimize!(model)

    obj_value = JuMP.objective_value(model)
    A_value = JuMP.value.(A)

    if verbose
        println("Objective value: ", obj_value)
        println("x = ", A_value)
	matshow(A_value[ns+1:2*ns,1:ns])
	colorbar()
	show()
    end

end

example_basic(verbose = true)


