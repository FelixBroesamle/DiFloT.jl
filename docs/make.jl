using Documenter
using DiFloT

makedocs(
    sitename = "DiFloT.jl",
    modules  = [DiFloT],
    strict = false,
    format   = Documenter.HTML(prettyurls = get(ENV, "CI", nothing) == "true"),
    pages    = [
        "Home" => "index.md",
        "File Format" => "format.md",
        "API Reference" => "api.md"
    ]
)