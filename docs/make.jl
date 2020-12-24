push!(LOAD_PATH, "..")

using Documenter, DocumenterCitations
using Literate

using LongwaveModePropagator
using LongwaveModePropagator.Fields

bib_filepath = joinpath(dirname(@__FILE__), "longwavemodepropagator.bib")
bib = CitationBibliography(bib_filepath)

#==
Generate examples
==#

const EXAMPLES_DIR = joinpath(@__DIR__, "..", "examples")
const OUTPUT_DIR   = joinpath(@__DIR__, "src", "generated")

examples = [
    # "meshgrid.jl",
    # "meshgrid2.jl",
    "wavefieldintegration.jl"
]


for example in examples
    example_filepath = joinpath(EXAMPLES_DIR, example)
    Literate.markdown(example_filepath, OUTPUT_DIR, documenter=true)
end

#==
Organize page hierarchies
==#

example_pages = [
    "generated/meshgrid.md",
    "generated/meshgrid2.md",
    "generated/wavefieldintegration.md"
]

library_pages = [
    "lib/public.md",
    "lib/internals.md"
]

pages = [
    "Home" => "index.md",
    "Examples" => example_pages,
    "Library" => library_pages,
    "References" => "references.md"
]

#==
Build and deploy docs
==#

format = Documenter.HTML(
    collapselevel = 1,
    prettyurls=get(ENV, "CI", nothing)=="true"
)

makedocs(bib,
    sitename = "LongwaveModePropagator.jl",
    authors = "Forrest Gasdia",
    format = format,
    pages = pages,
    modules = [LongwaveModePropagator]
)
