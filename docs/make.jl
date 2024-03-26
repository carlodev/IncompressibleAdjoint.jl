push!(LOAD_PATH,joinpath("..","src"))

using Documenter, DocumenterCitations, IncompressibleAdjoint

bib = CitationBibliography(joinpath(@__DIR__, "src", "docs.bib"); style=:numeric)


makedocs(bib;
    sitename = "IncompressibleAdjoint.jl",
    modules = [IncompressibleAdjoint],
    pages = [
        "Introduction" => "index.md",    
        "API information" => "api_info.md",
        "References" => "references.md",
    ],
)

deploydocs(
    repo = "github.com/carlodev/IncompressibleAdjoint.jl",
    push_preview = true,
)
