# 22i - Hypergraph inference
Here is a **rough** documentation for some of the functions defined here. 
The documantation is incomplete and will stay so.

## Hyper-Kuramoto
The function `hyper_k` (line `74` in [`hyper_kuramoto.jl`](https://github.com/r0b1nu5/bucket/blob/master/22i-hyper-inf/hyper_kuramoto.jl)) simulates the dynamics of a set of $N$ coupled Kuramoto oscillators on an hypergraph of order 3, following the dynamics of [Delabays et al. (2024), Eq. (3)](https://arxiv.org/abs/2402.00078), using a 4th-order Runge-Kutta scheme. 
Additive noise can be added, which is space-uncorrelated and can be time-correlated. 
<!--
$$
 \dot{\theta}_i = \omega_i - \sum_{j=1}^N a^{(2)}_{ij}\sin(\theta_i-\theta_j-\phi_2) - \sum_{j,k=1}^Na^{(3)}_{ijk}\sin(2\theta_i-\theta_j-\theta_k-\phi_3)\, ,
$$
-->


    function hyper_k(A2, A3, ω, θ0; τ0=1., ξ0=0., ϕ2=0., ϕ3=0., h=.01, max_iter=10_000, tol=1e-6)

**INPUTS:**

`A2::Matrix{Float64}`: Adjacency matrix of pairwise interactions. Either the actual adjacency $N\times N$ matrix or the list of edges is supported. A list of edges should be given as an array where each row has three columns $[i, j, a^{(2)}_{ij}]$. 

`A3::Union{Array{Float64,3},Matrix{Float64}}`: Adjacency tensor of triadic interactions. Either the actual adjacency $N\times N\times N$ tensor or the list of edges is supported. A list of hyperedges should be given as an array where each row has four columns $[i, j, k, a^{(3)}_{ijk}]$. 

`ω::Vector{Float64}`: Vector of the oscillators' natural frequencies. 

`θ0::Vector{Float64}`: Vector of initial oscillators' angles. 

`τ0::Float64=1.`: Time correlation of the noise.
Default value is `1.0`. 

`ξ0::Float64=0.`: Magnitude of the noise. 
Default value is `0.0` (i.e., no noise). 

`ϕ2::Float64=0.`: Phase frustration in the pairwise interactions. 
Default value is `0.0` (no frustration). 

`ϕ3::Float64=0.`: Phase frustration in the triadic interactions. 
Default value is `0.0` (no frustration). 

`h::Float64=.01`: Time step in the Runge-Kutta numerical scheme. 
Default value is `0.01`. 

`max_iter::Int64=10_000`: Maximum number of iterations.
Simulation stops after `max_iter` steps. 
Default value is `10_000`. 

`tol::Float64=1e-6`: Stopping criterion. 
If the maximal angle update over one time step is smaller than `tol`, the simulation stops. 
Default value is `1e-6`. 

**OUTPUT:**

`θs::Matrix{Float64}`: Trajectories of the oscillators' angles. 
The matrix is of size $N\times T$, where $T$ is the number of time steps. 

`dθs::Matrix{Float64}`: Time derivatives of the oscillators' angles at each time step. 
The matrix is of size $N\times T$, where $T$ is the number of time steps. 

**Remarks:**
- Packages needed are: `DelimitedFiles`, `Distributions`, `LinearAlgebra`, `Statistics`.
- A folder named `temp` is required in the running directory for memory purposes. 
- Noise generation calls the function `cnoise` in [`cnoise.jl`](https://github.com/r0b1nu5/bucket/blob/master/22i-hyper-inf/cnoise.jl). 
Line `89` in [`hyper_kuramoto.jl`](https://github.com/r0b1nu5/bucket/blob/master/22i-hyper-inf/hyper_kuramoto.jl) can be commented out for noiseless simulations. 


## Hyper-Lorenz
Trajectories of coupled Lorenz oscillators are simulated by running the Python script [`Lorenz_hypergraph.py`](https://github.com/r0b1nu5/bucket/blob/master/22i-hyper-inf/Lorenz_hypergraph.py). 
The dynamics is given in [Delabays et al. (2024), page 4](https://arxiv.org/abs/2402.00078). 
<!--
The dynamics is given by
$$
 \dot{x}_i = \sigma(y_i-x_i) + \sum_{j=1}^N a_{ij}^{(2)}(x_j-x_i) + \sum_{j,k=1}^N a_{ijk}^{(3)}\left(x_kx_j^2-x_i^3\right)\\
 \dot{y}_i = x_i(\rho - z_i) - y_i \\
 \dot{z}_i = x_iy_i - \beta z_i\, ,
$$
with $\sigma = 10$, $\rho = 28$, and $\beta = 8/3$. 
-->

Simulation parameters can be tuned at lines `42-51` of [`Lorenz_hypergraph.py`](https://github.com/r0b1nu5/bucket/blob/master/22i-hyper-inf/Lorenz_hypergraph.py). 

The interaction hypergraph is randomly generated using the [XGI library](https://xgi.readthedocs.io/en/stable/). 
The seed of the random hypergraph generator can be changed at line `54` of [`Lorenz_hypergraph.py`](https://github.com/r0b1nu5/bucket/blob/master/22i-hyper-inf/Lorenz_hypergraph.py)


## Hypergraph inference
A (hopefully) clean version of the inference algorithm should be available [here](https://github.com/TaylorBasedHypergraphInference/THIS). 

## Test scripts
### Hyper-Kuramoto
The latest tests of THIS for Kuramoto oscillators were run in [`test_hyper_inf_new.jl`](https://github.com/r0b1nu5/bucket/blob/master/22i-hyper-inf/test_hyper_inf_new.jl).
There is a lot of rubbish in this script, but I want to keep it as it serves as a backup of past tests. 
Furthermore, in order to avoid overloading GitHub, data are kept local. 
I give here a rough description of what each section does. 

**Lines:**

`11-14`: Initialize the problem parameters (`n`: number of nodes, `T`: number of time steps in the simulation, `iters`: number of time steps to consider in the various inferences). 

`17-20`: Initialize the adjacency tensors. 

`23-66`: Choose the type of random hypergraph to consider (uncomment desired). 

`68-92`: Generate or load the random hypergraph. When the `ntw` tag has a `-py` suffix, it means that the hypergraph was generated with one of the Python scripts [`gen_rand_hyperg.py`](https://github.com/r0b1nu5/bucket/blob/master/22i-hyper-inf/gen_rand_hyperg.py) or [`gen_rand_simplicial.py`](https://github.com/r0b1nu5/bucket/blob/master/22i-hyper-inf/gen_rand_simplicial.py).

`172-174`: Initialize the inference parameters.

For length of time series in `iters`:

`180`: Run THIS.

`181-182`: Retrieve the adjecency tensors from the inference. 

`220-232`: Compute and plot the ROC curves. 

### Hyper-Lorenz
The latest tests of THIS for Lorenz oscillators were run in [`test_hyperinf_lorenz_bis.jl`](https://github.com/r0b1nu5/bucket/blob/master/22i-hyper-inf/test_hyperinf_lorenz_bis.jl). 
I give here a rough description of what each section does. 

**Lines:** 

`5-14`: Initialize the script. 

`34-56`: Loading the time series `Xs` and the time derivatives `Ys`. 
We take `m` chunks of time series, each of length `dT` time steps and `T` time steps appart from each other. 
The derivatives are estimated by finite differences. 

`59-71`: Running THIS twice, once the raw version (line `66`) and once with pre-filtering (line `69`). 

`76-124`: Extracting the inferred adjacency tensors. 

`131-175`: Loading the ground truth adjacency tensors. 

`177-198`: Computing and plotting the ROC curves. 

