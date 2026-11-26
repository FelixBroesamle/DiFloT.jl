using Test

@testset "README executable Julia examples" begin
    root = normpath(joinpath(@__DIR__, ".."))
    readme = read(joinpath(root, "README.md"), String)

    blocks = collect(eachmatch(
        r"(?ms)^```julia[ \t]*\r?\n(.*?)^```[ \t]*\r?$",
        readme
    ))

    @test !isempty(blocks)

    for (i, block) in enumerate(blocks)
        @testset "Example $i" begin
            mktempdir() do dir
                filepath = joinpath(dir, "snippet.jl")
                write(filepath, block.captures[1])

                # Execute each example in a separate module.
                sandbox = Module(gensym(:ReadmeExample))
                Core.eval(sandbox, :(using Base))

                # Resolve relative paths from the repository root.
                cd(root) do
                    Base.include(sandbox, filepath)
                end

                @test true
            end
        end
    end
end