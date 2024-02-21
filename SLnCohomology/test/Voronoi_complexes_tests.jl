function test_cell_data(N)
    #import list of perfect forms
    A_3 = [2 -1 0
    -1 2 -1
    0 -1 2]

    D_4 = [2 0 1 0
    0 2 -1 0
    1 -1 2 -1
    0 0 -1 2]

    A_4 = [2 -1 0 0
    -1 2 -1 0
    0 -1 2 -1
    0 0 -1 2]

    D_5 = [2 0 1 0 0
    0 2 -1 0 0
    1 -1 2 -1 0
    0 0 -1 2 -1
    0 0 0 -1 2]

    A_5_plus3 = [6 -3 0 0 0
    -3 6 -3 0 3
    0 -3 6 -3 0
    0 0 -3 6 0
    0 3 0 0 4]

    A_5 = [2 -1 0 0 0
    -1 2 -1 0 0
    0 -1 2 -1 0
    0 0 -1 2 -1
    0 0 0 -1 2]

    if N == 3
        forms = [A_3]
    elseif N == 4
        forms = [D_4, A_4]
    elseif N == 5
        forms = [D_5, A_5_plus3, A_5]
    end

    # compute the Voronoi_cells for n=3,4,5
    cells_SLn = SLnCohomology.Voronoi_cells(N,forms)
    cells_SLn_old = deserialize(joinpath(@__DIR__, "./old_data/sl"*string(N)*"_cells.sjl"))

    @test_skip keys(cells_SLn) == keys(cells_SLn_old) # activate if want to check that we have cells in the same dimensions

    # check whether we have the same orbits of cells
    for dimension in keys(cells_SLn_old)
        @test length(cells_SLn[dimension]) == length(cells_SLn_old[dimension])
        for old_cell in cells_SLn_old[dimension]
            @test SLnCohomology.orbit_in_list(old_cell, cells_SLn[dimension])
        end
    end
end

@testset "sl3_cell_data" begin
    test_cell_data(3)
end

@testset "sl4_cell_data" begin
    test_cell_data(4)
end

@testset "sl5_cell_data" begin
    test_cell_data(5)
end