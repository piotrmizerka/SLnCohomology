function test_cell_data(N)
    # compute the Voronoi_cells
    cells_SLn = SLnCohomology.cells_sln(N)
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