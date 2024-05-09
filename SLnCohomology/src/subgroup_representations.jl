# Computes a flip-permutation representation of a specific subgroup
# of SL(N,p). Available options for (N,p): (3,3) and (4,2).
function flip_permutation_representation(
    N::Integer,
    p::Integer
)
    gens_H, gens_index_two_H = symmetric_subgroups_gens(N,p)
    H_gens_expr = subgroup_gens_expression(gens_H, p) # Express all elements of H as words in the distinguished generators of H.
    index_two_H = keys(subgroup_gens_expression(gens_index_two_H, p)) # Compute the whole index two subgroup of H.
    
    # Compute permutation representation of H (by defining it on the generators of and stored in perm_rep_gens variable). 
    # Return its degree as well.
    if N == 3 && p == 3
        # Linear rep of order 2 of Q₈ embedded in SL(3,3).
        a, b, a_inv, b_inv = gens_H
        perm_rep_gens = Dict(
            a => Matrix(Permutations.Permutation([[1]])), 
            b => Matrix(Permutations.Permutation([[1]])),
            a_inv => Matrix(Permutations.Permutation([[1]])),
            b_inv => Matrix(Permutations.Permutation([[1]]))
        )
        deg = 1
    elseif N == 4 && p == 2
        # Define the permutation representation on the generators s and t and their inverses
        s, t, s_inv, t_inv = gens_H
        perm_rep_gens = Dict(
            s => Matrix(Permutations.Permutation([[1,11,14,4],[3,6,12,9],[2,13],[5,10],[7],[8],[15]])), 
            t => Matrix(Permutations.Permutation([[4,7,11,10,9,5],[1,3,12],[2,15,13],[6,8],[14]])), 
            s_inv => Matrix(Permutations.Permutation([[4,14,11,1],[9,12,6,3],[2,13],[5,10],[7],[8],[15]])), 
            t_inv => Matrix(Permutations.Permutation([[5,9,10,11,7,4],[12,3,1],[13,15,2],[6,8],[14]]))
        )
        deg = 15
    end

    # Compute the flip-permutation representation given by perm_rep_gens, H_gens_expr and index_two_H.
    subgroup_rep = flip_representation(perm_rep_gens, H_gens_expr, index_two_H)

    return subgroup_rep, deg
end

# Given values of a representation π on generators (assumed to be symmetric) of a subgroup H of slnp
# possessing a subgroup K of index 2, compute the representation for the whole subgroup H by
# multiplying the representation π by -1 for elements not belonging to K.
# Store the result in a dictionary.
function flip_representation(
    π, # values of the representation π on the generators of H
    sgp_gens_expr,
    K
)
    result = Dict()
    for h in keys(sgp_gens_expr) # note that H = keys(sgp_gens_expr)
        val = I
        for s in sgp_gens_expr[h]
            val *= π[s]
        end
        factor = (h in K ? 1 : -1)
        result[h] = factor*val
    end
    return result
end

# Given matrix generators gens_ (assumed to be symmetric in the sense that for all its elements 
# their inverses are also in gens_) of a subgroup H of SL(N,p), express all elements of H as words 
# in gens_ (not necessarily uniquely). Save the output in a dictionary.
# Note that the subgroup H itself and slnp are not arguments of this function.
function subgroup_gens_expression(
    gens_,
    p::Integer
)
    result = Dict(x => [x] for x in gens_)
    I_N = Matrix(UniformScaling(Int8(1)),size(first(gens_))[1],size(first(gens_))[1])
    result[I_N] = []
    while true
        old_result_length = length(result)
        for M in keys(result)
            for s in gens_
                sM, Ms = matrix_mod_p(s*M,p), matrix_mod_p(M*s,p)
                if !(sM in keys(result))
                    result[sM] = vcat([s],result[M])
                end
                if !(Ms in keys(result))
                    result[Ms] = vcat(result[M],[s])
                end
            end
        end
        if length(result) == old_result_length
            break
        end
    end
    return result
end

# Computes the generators of the subgroup H of SL(N,p) and its index two subgroup K.
function symmetric_subgroups_gens(N::Integer, p::Integer)
    if (N,p) == (3,3)
        # H equals Q₈ embedded in SL(3,3) by embedding its embedding into SL(2,3) into SL(3,3).
        a = [0 2 0;
             1 0 0;
             0 0 1]
        b = [1 1 0;
             1 2 0;
             0 0 1]
        a_inv = Int8.(AbstractAlgebra.lift.(inv(AbstractAlgebra.matrix(GF(p),a))))
        b_inv = Int8.(AbstractAlgebra.lift.(inv(AbstractAlgebra.matrix(GF(p),b))))
        gens_H = [a,b,a_inv,b_inv]

        # Compute index 2 subgroup of H. This subgroup, isomorphic to C₄, is generated by the following matrix c.
        c = [1 2 0;
             2 2 0;
             0 0 1]
        c_inv = Int8.(AbstractAlgebra.lift.(inv(AbstractAlgebra.matrix(GF(p),c))))
        gens_index_two_H = [c,c_inv]
    elseif (N,p) == (4,2)
        # The subgroup H of SL(4,2) we consider is generated by s and t below. H has order 576.
        s = [1 0 0 0;
             0 0 0 1;
             1 1 0 1;
             1 0 1 1]
        t = [0 1 1 0;
             0 1 1 1;
             1 1 1 1;
             0 0 1 1]
        s_inv = Int8.(AbstractAlgebra.lift.(inv(AbstractAlgebra.matrix(GF(p),s))))
        t_inv = Int8.(AbstractAlgebra.lift.(inv(AbstractAlgebra.matrix(GF(p),t))))
        gens_H = [s,t,s_inv,t_inv]

        # Compute index 2 subgroup of H. This subgroup is generated by the following matrices
        # (we also include their inverses in the case they are not idempotent).
        a = [1 0 1 0;
             0 1 0 0;
             0 0 1 0;
             0 0 0 1]
        b = [0 1 1 1;
             0 1 0 0;
             1 1 1 1;
             0 0 0 1]
        c = [1 1 0 0;
             1 0 0 0;
             1 1 0 1;
             1 0 1 1]
        d = [0 1 1 1;
             1 0 1 1;
             0 0 1 0;
             0 0 0 1]
        e = [1 0 1 0;
             0 1 1 0;
             0 0 1 0;
             0 0 0 1]
        f = [0 1 0 1;
             0 1 0 0;
             0 0 1 0;
             1 1 0 0]
        g = [1 0 0 0;
             1 0 1 1;
             0 0 1 0;
             1 1 1 0]
        b_inv = Int8.(AbstractAlgebra.lift.(inv(AbstractAlgebra.matrix(GF(p),b))))
        c_inv = Int8.(AbstractAlgebra.lift.(inv(AbstractAlgebra.matrix(GF(p),c))))
        g_inv = Int8.(AbstractAlgebra.lift.(inv(AbstractAlgebra.matrix(GF(p),g))))
        gens_index_two_H = [a,b,c,d,e,f,g,b_inv,c_inv,g_inv]
    end

    return gens_H, gens_index_two_H
end
