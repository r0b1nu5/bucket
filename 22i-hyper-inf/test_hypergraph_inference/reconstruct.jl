using PyPlot, DelimitedFiles, ROC

include("./basis_expansion.jl")

function reconstruct(X, DX, NODE, connectivity, BASIS="polynomial_diff", ORDER=6, MODEL="kuramoto1")
    # reconstruct(MODEL, NODE, BASIS, ORDER) returns a ranked list of the inferred incoming connections.

    # Parameters
    # ------------------
    # MODEL: Dynamic model employed. This is only used to specify whether the
    #        time series come from 1D systems like kuramoto1 or 3D systems like
    #        roessler. Thus, it is not used during the actual reconstruction.
    # NODE:  Unit upon the reconstruction takes place.
    # BASIS: Type of basis employed. Currently, polynomial, polynomial_diff,
    #        power_series, fourier, fourier_diff and RBF are supported. For
    #        more detailed information, please see 'Functions/basis_expansion.m'
    #        and Table I in the main manuscript.
    # ORDER: Number of basis in the expansion.
    #
    # Input type
    # ------------------
    # MODEL: string
    # NODE:  integer
    # BASIS: string
    # ORDER: integer
    #
    # Output
    # ------------------
    # list: Sequence of inferred interactions in the order such were detected.
    # cost: Fitting cost for all inferred interactions in the order such were
    #       detected.
    # FPR:  False positives rate for the reconstruction.
    # TPR:  True positives rate for the reconstruction.
    # AUC:  Quality of reconstruction measured in AUC scores.
    #
    # Example
    # ------------------
    # reconstruct("michaelis_menten", 10, "polynomial", 6) reconstructs the
    # connectivity of unit 10 using polynomials up to power 6 as basis functions.
    #
    # Accompanying material to "Model-free inference of direct interactions 
    # from nonlinear collective dynamics".
    #
    # Author: Jose Casadiego
    # Date:   May 2017

    # Stopping criterium: decrease it to recover longer list of possible links
    th = 0.0001

    models = ["kuramoto1", "kuramoto2", "michaelis_menten", "roessler"]
    bases = ["polynomial", "polynomial_diff", "fourier", "fourier_diff", "power_series", "RBF"]

    if MODEL ∉ models
        println("ERROR: MODEL must be a valid string: kuramoto1, kuramoto2, michaelis_menten, roessler")
    elseif BASIS ∉ bases
        println("ERROR: BASIS must be a valid string: polynomial, polynomial_diff, fouries, fourier_diff, power_series, RBF")
    else
        println("Initiating reconstruction...")
        println("Reading data...")
#        data = readdlm("Data/data.dat")
#        data = data'
#        connectivity = readdlm("Data/connectivity.dat")
#        ts_param = readdlm("Data/ts_param.dat")

        N,T = size(X)
	M = T
#=
        x = data

        # Estimating time derivatives and constructing input matrices
        println("Estimating time derivatives and constructing input matrices...")
        Xtemp = zeros(N,0)
        DX = zeros(N,0)
        for s = 1:S
            Ytemp = zeros(N, M)
            DY = zeros(N, M)
            for n = 1:N
                for m = 1:M
                    Ytemp[n, m] = (x[n,m + (s-1)*(M+1)] + x[n,m+1+(s-1)*(M+1)])*0.5
                    DY[n, m] = (-x[n,m + (s-1)*(M+1)] + x[n,m+1+(s-1)*(M+1)])*1/1.0
                end
            end
            Xtemp = [Xtemp Ytemp]
            DX = [DX DY]
        end
=#

#######################################################################
# STILL TODO...
        if MODEL == "roessler"
            # Construction of connectivity matrix for Roessler oscillators
            # including y and z variables
            Ns = ceil(Int, N / 3)
            connectivity2 = zeros(Ns, N)

            for i = 1:Ns
                for j = 1:Ns
                    connectivity2[i, 3 * (j - 1) + 1] = connectivity[i, j]
                    if i == j
                        connectivity2[i, 3 * (j - 1) + 2] = 1
                        connectivity2[i, 3 * (j - 1) + 3] = 1
                    end
                end
            end

            X = Xtemp

            # Beginning of reconstruction algorithm
            println("Performing ARNI...")
            Y = basis_expansion(X, ORDER, BASIS, NODE)
            nolist = 1:N
            list = Int[]
            cost = Float64[]
            b = true
            vec = zeros(Float64, N)
            while !isempty(nolist) && b
                # Composition of inferred subspaces
                Z = zeros(Float64, 0, M)
                for n in list
                    Z = vcat(Z, Y[:, :, n])
                end

                # Projection on remaining composite spaces
                P = zeros(Float64, length(nolist), 2)
                cost_err = zeros(Float64, length(nolist))
                for (idx, n) in enumerate(nolist)
                    # Composition of a possible space
                    R = vcat(Z, Y[:, :, n])
                    # Error of projection on possible composite space
                    P[idx, 1] = std(DX[3 * (NODE - 1) + 1, :] - DX[3 * (NODE - 1) + 1, :] * pinv(R) * R)
                    P[idx, 2] = n
                    # Fitting cost of possible composite space
                    cost_err[idx] = 1 / M * norm(DX[3 * (NODE - 1) + 1, :] - DX[3 * (NODE - 1) + 1, :] * pinv(R) * R)
                    R = []
                end

                if std(P[:, 1]) < th
                    b = false
                    break
                else
                    # Selection of composite space which minimizes projection error
                    MIN, block = findmin(P[:, 1])
                    push!(list, P[block, 2])
                    nolist = filter(x -> x != P[block, 2], nolist)
                    vec[P[block, 2]] = MIN
                    push!(cost, cost_err[block])
                end
            end
            # End of reconstruction algorithm

            adjacency = connectivity2
            adjacency[adjacency .!= 0] = 1

            # Evaluation of results via AUC score
            FPR, TPR, _, AUC = perfcurve(adjacency[NODE, :], abs.(vec), 1)
            println("Reconstruction has finished!")
            println("Quality of reconstruction:")
            println(AUC)
            return list, cost, FPR, TPR, AUC

#######################################################################

        else
#=
            if MODEL == "kuramoto1" || MODEL == "kuramoto2"
                # Transforming data coming from phase oscillators
                X = mod.(Xtemp, 2 * π)
            else
                X = Xtemp
            end
=#
            # Beginning of reconstruction algorithm
            println("Performing ARNI...")
            Y = basis_expansion(X, ORDER, BASIS, NODE)
            nolist = 1:N
            list = Int64[]
            cost = Float64[]
            b = true
            vec = zeros(N)
            while !isempty(nolist) && b
                # Composition of inferred subspaces
                Z = zeros(0, T)
                for n in list
                    Z = [Z;Y[:, :, n]]
                end

                # Projection on remaining composite spaces
                P = zeros(length(nolist), 2)
                cost_err = zeros(length(nolist))
                for (idx, n) in enumerate(nolist)
                    # Composition of a possible space
                    R = [Z;Y[:, :, n]]
                    # Error of projection on possible composite space
                    P[idx, 1] = std(DX[[NODE,], :] - DX[[NODE,], :] * pinv(R) * R)
                    P[idx, 2] = n
                    # Fitting cost of possible composite space
                    cost_err[idx] = 1 / M * norm(DX[[NODE,], :] - DX[[NODE,], :] * pinv(R) * R)
                    R = Float64[]
                end

                if std(P[:, 1]) < th
                    b = false
                    break
                else
                    # Selection of composite space which minimizes projection error
                    MIN, block = findmin(P[:, 1])
#    @info "SHOW MIN: $MIN"
                    push!(list, Int64(P[block, 2]))
                    nolist = setdiff(nolist,list)
                    vec[Int64(P[block, 2])] = MIN
                    push!(cost, cost_err[block])
                end
            end
            # End of reconstruction algorithm

            adjacency = Float64.(abs.(connectivity) .> 1e-8)


            # Adding degradation rate to true adjacency matrix of
            # Michaelis-Menten systems
            if MODEL == "michaelis_menten"
                adjacency[1:N + 1:N * N] .= 1
            end

            # Evaluation of results via AUC score
            xxx = roc(abs.(vec), adjacency[NODE, :])
            fpr = xxx.FPR
            tpr = xxx.TPR
            auc = ROC.AUC(xxx)
            println("Reconstruction has finished!")
            println("Quality of reconstruction:")
            println(auc)
            return vec, list, cost, fpr, tpr, auc
        end
    end
end

