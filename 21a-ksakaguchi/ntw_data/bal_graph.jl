using LinearAlgebra 

L = [1 -1 0 0 0 0;
     0 1 -1 0 0 0;
     0 0 1 0 -1 0;
     0 -1 0 1 0 0;
     0 0 0 -1 2 -1;
     -1 0 0 0 0 1.]

B = [1 0 0 0 0 0 -1;
     -1 1 -1 0 0 0 0;
     0 -1 0 1 0 0 0;
     0 0 1 0 -1 0 0;
     0 0 0 -1 1 1 0;
     0 0 0 0 0 -1 1.]

Bout = [1 0 0 0 0 0 0;
	0 1 0 0 0 0 0;
	0 0 0 1 0 0 0;
	0 0 1 0 0 0 0;
	0 0 0 0 1 1 0;
	0 0 0 0 0 0 1.]

Bin = [0 0 0 0 0 0 1;
       1 0 1 0 0 0 0;
       0 1 0 0 0 0 0;
       0 0 0 0 1 0 0;
       0 0 0 1 0 0 0;
       0 0 0 0 0 1 0.]

u = abs.(real.(eigvecs(L')[:,1]))
w = [u[1],u[2],u[4],u[3],u[5],u[5],u[6]]
W = diagm(0 => w)



