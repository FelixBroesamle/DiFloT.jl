using Documenter
using DiFloT

makedocs(
    sitename = "DiFloT.jl",
    modules  = [DiFloT],
    format   = Documenter.HTML(
        prettyurls = get(ENV, "CI", nothing) == "true",
        assets = ["assets/favicon.ico"]
    ),
    pages    = [
        "Home" => "index.md",
        "User Guide" => [
            "File Formats" => "format.md",
            "End-to-End Workflows" => "workflows.md"
        ],
        "Developer Guide" => [
            "API Reference" => "api.md"
        ]
    ],
    warnonly = [:missing_docs]
)

deploydocs(
    repo = "github.com/FelixBroesamle/DiFloT.jl.git",
    devbranch = "main",
)