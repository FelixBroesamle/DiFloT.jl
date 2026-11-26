using Test
using DiFloT
import JuMP

include("helpers.jl")

@testset "DiFloT.jl" begin
    include("core.jl")
    include("highs_backend.jl")
    include("algorithms.jl")
    include("readme_examples.jl")

    flag = get(ENV,"DIFLOT_TEST_ORTOOLS","false")
    flag in ("true","false") || error("DIFLOT_TEST_ORTOOLS must be true or false.")
    if flag == "true"
        include("ortools_backend.jl")
    else
        @info "OR-Tools integration not requested; set DIFLOT_TEST_ORTOOLS=true to require it."
    end
end
