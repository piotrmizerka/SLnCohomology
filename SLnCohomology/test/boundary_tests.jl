function test_boundary_data(N)
    cells_sln = deserialize(joinpath(@__DIR__, "../scripts/differentials_computation/precomputed_cells/sl"*string(N)*"_cells.sjl"))

    # computing the chain complex with the new functions
    oriented_cells_sln = SLnCohomology.oriented_cells_dict(cells_sln)
    stabiliser_new = SLnCohomology.stabilisers_dict(oriented_cells_sln) # compute the stabilisers
    boundaries_new  = SLnCohomology.boundaries_dict(oriented_cells_sln) # compute the boundaries

    sln_data_old = deserialize(joinpath(@__DIR__, "../scripts/differentials_computation/precomputed_boundaries/sl"*string(N)*"_bound_stab_old.sjl"))
    stabiliser_old = sln_data_old["stabilisers"]
    boundaries_old = sln_data_old["boundaries"]
    for dimension in keys(stabiliser_new)
        for index in eachindex(stabiliser_new[dimension]) # these are the cells in dimension
            set_stabilisers_new = Set(stabiliser_new[dimension][index])
            set_stabilisers_old = Set(stabiliser_old[dimension][index])
            @test set_stabilisers_new == set_stabilisers_old
        end
    end

    for dimension in keys(boundaries_new)
        for index in eachindex(boundaries_new[dimension]) # these are the cells in dimension
            # first check whether the cell has the same boundaries (without orientations...)
            new_cell = boundaries_new[dimension][index]
            old_cell = boundaries_old[dimension][index]
            
            new_cell_boundaries_list = []
            for boundary in new_cell
                push!(new_cell_boundaries_list,boundary["orbit_standard_cell"])
            end

            old_cell_boundaries_list = []
            for boundary in old_cell
                push!(old_cell_boundaries_list,boundary["orbit_standard_cell"])
            end

            @test Multiset(new_cell_boundaries_list) == Multiset(old_cell_boundaries_list)

            # now check whether the orientations agree
            # the test below makes the above test redundant (but is more complicated)
            for boundary_new in new_cell
                orbit_new = boundary_new["orbit_standard_cell"]
                orbit_coset_new = boundary_new["orbit_coset_with_orientation"]
                # try the two different sign conventions (this is not sufficient to say that everything is the same, but necessary)
                orbit_coset_new_with_sign1 = []
                for (matrix,sign) in orbit_coset_new
                    push!(orbit_coset_new_with_sign1,(matrix,sign))
                end
                orbit_coset_new_with_sign2 = []
                for (matrix,sign) in orbit_coset_new
                    push!(orbit_coset_new_with_sign2,(matrix,(-1)*sign))
                end
                found_same_boundary = 0
                for boundary_old in old_cell
                    # check whether we find that boundary in the old boundary list
                    if boundary_old["orbit_standard_cell"] == orbit_new
                        orbit_coset_old = boundary_old["orbit_coset_with_orientation"]
                        if Set(orbit_coset_new_with_sign1) == Set(orbit_coset_old) || Set(orbit_coset_new_with_sign2) == Set(orbit_coset_old)
                            found_same_boundary +=1
                        end
                    end
                end
                @test found_same_boundary == 1
            end
        end
    end
end

@testset "sl3_boundary_data" begin
    test_boundary_data(3)
end

@testset "sl4_boundary_data" begin
    test_boundary_data(4)
end