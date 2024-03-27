push!(LOAD_PATH,joinpath("..","src"))

using Documenter, DocumenterCitations, IncompressibleAdjoint

bib = CitationBibliography(joinpath(@__DIR__, "src", "docs.bib"); style=:numeric)
pages=[ "Introduction" => "index.md",    
        "API information" => "api_info.md",
        "References" => "references.md",
    ]

makedocs(;
    sitename = "IncompressibleAdjoint.jl",
    modules = [IncompressibleAdjoint],
    format = Documenter.HTML(),
    pages = pages,
    plugins=[bib],
)

deploydocs(
    repo = "github.com/carlodev/IncompressibleAdjoint.jl",
    push_preview = true,
)
