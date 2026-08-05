using Test

@testset "MicroclimateMapper" begin
    @testset "landcover"      begin include("landcover.jl") end
    @testset "terrain"        begin include("terrain.jl") end
    @testset "flow graph"     begin include("flow_graph.jl") end
    @testset "output writing" begin include("output_writing.jl") end
    @testset "micro raster"   begin include("micro_raster.jl") end
    @testset "micro vector"   begin include("micro_vector.jl") end
    @testset "soil pedotransfer" begin include("soil_pedotransfer.jl") end
end
