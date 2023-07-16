MANUAL FOR THE CODE ACCOMPANYING THE MANUSCRIPT: 
On the Robustness of Democratic Electoral Processes to Computational Propaganda

Last edited: July 12, 2023

Model_execution.py
        outcome
        type_natural_opinion

Nat_opn_generator.py
        admissible_summits
        gen_rand_perm_vect
        Nat_opn_2
        rand_opinion
        rand_opinion_pos
        Random_generation
        Random_generation_pos
        slide_summit_vec
        unit_simplex_arb
        vertices

===============================================================================
===============================================================================
**Model_execution.py**
Running this script runs an example of the opinion formation simulation with choice of parameters. 
===============================================================================
===============================================================================
*outcome(x0, epsilon, W)*
Computes the final opinions from the natural opinion x0, the communication distance epsilon, and the influence W. 
INPUT:
	x0 (vector, matrix): natural opinions
	epsilon (float): communication distance
	W (vector, matrix): influence vector/matrix
OUTPUT:
	y (vector, matrix): final opinions

===============================================================================
*type_natural_opinion(Num_ppl, Num_party, NAT_TYPE)*
Generates the natural opinions by calling the appropriate function. In the opinion is "Bigaussian", the mean, Delta, ds parameters have to be specified beforehand in the local environment. If Num_party is 2, NAT_TYPE can be 'Bigaussian' or 'Uniform', otherwise it can only be 'Uniform'. 'Uniform' draws the opinions from the appropriate simplex described in the manuscript's SI. 
INPUT:
	Num_ppl (int): number of agents
	Num_party (int): number of parties
	NAT_TYPE (string): type of distribution from which the opinions will be drawn 
OUTPUT: 
	x0 (vector, matrix): natural opinions

===============================================================================
===============================================================================
**Nat_opn_generator.py**
===============================================================================
===============================================================================
*admissible_summits(p)*
Defines the vertices of the admissible opinion space.
INPUT:
	p (int): number of parties
OUTPUT:
	V (matrix): vertices of the admissible opinion space

===============================================================================
*gen_rand_perm_vect(x)*
Generates a random permutation of the vector x.
INPUT:
	x (vector)
OUTPUT:
	p (vector)

===============================================================================
*Nat_opn_2(mu, Delta, sd, num_ppl, R=50)*
Generates gaussian and bigaussian random natural opinion for 2 parties
INPUT:
	mu (float): mean of the whole distribution
	Delta (float): polarization coefficient, i.e., distance between the gaussian peaks
	sd (float): standard deviation of the gaussian(s)
	num_ppl (int): number of opinions to be generated
	R (float): percentage of agents in each peak
OUTPUT:
	x0 (matrix): array of natural opinions (num_ppl,2)

===============================================================================
*rand_opinion(p, V)*
Draws one random opinion uniformly in the admissible space. The of vertices is already given in 'V'. One can get 'V' by using the function 'admissible_summits'. 
INPUT:
	p (int): number of parties
	V (matrix): vertices of the admissible space
OUTPUT:
	x (matrix): random opinion

===============================================================================
*rand_opinion_pos(p, V, k)*
Same as 'rand_opinion' with a fixed party k.
INPUT:
	p (int): number of parties
	V (matrix): vertices of the admissible space
	k (int): index of the party where the opinion is drawn
OUTPUT:
	x (matrix): random opinion

===============================================================================
*Random_generation(num_ppl, p)*
Generates 'num_ppl' random opinions with p parties.
INPUT:
	num_ppl (int): number of agents
	p (int): number of parties
OUTPUT:
	X (matrix): random opinions
===============================================================================
*Random_generation_pos(num_ppl, p, k)*
Generates 'num_ppl' random opinions in party k among p parties.
INPUT:
	num_ppl (int): number of agents
	p (int): number of parties
	k (int): index of the party
OUTPUT:
	X (matrix): random opinions

===============================================================================
*slide_summit_vec(x, a)*
Slides the point x (typically in the simplex) towards the barycenter of the unitary p-simplex. For a = 1, the summit does not move, and for a = 0, the summit reaches the barycenter.
INPUT:
	x (vector): point to be moved towards the barycenter of the simplex
	a (float): factor by which the point is slided
OUTPUT:
	y (vector): moved point

===============================================================================
*unit_simplex_arb(S)*
Draws a random point uniformly from an arbitrary simplex, defined by the columns of S. Follows the idea presented in doi.org/10.13140/RG.2.1.3807.6968.
INPUT:
	S (matrix): list of vertices of the simplex
OUTPUT:
	x (vector): random point in S

===============================================================================
*vertices(p)*
Generates the list of vertices of interest in the unitary simplex with 'p' vertices.
INPUT:
	p (int): number of parties
OUTPUT:
	vertex (matrix): list of the simplex's vertices
===============================================================================


