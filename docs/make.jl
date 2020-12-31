".." in LOAD_PATH || push!(LOAD_PATH, "..")  # from docs/

using Documenter, DocumenterCitations
using Literate

using LongwaveModePropagator
using LongwaveModePropagator.Fields

bib_filepath = joinpath(dirname(@__FILE__), "longwavemodepropagator.bib")
bib = CitationBibliography(bib_filepath)

#==
Generate examples
==#

const ROOT_DIR = abspath(@__DIR__, "..")
const EXAMPLES_DIR = joinpath(ROOT_DIR, "examples")
const BUILD_DIR = joinpath(ROOT_DIR, "docs", "build")
const OUTPUT_DIR = joinpath(ROOT_DIR, "docs", "src", "generated")

# repo_root_url is only assigned with CI, so we have to specify it for local builds
# can't use abspath on Windows (no colons) so we calculate relpath
const REPO_ROOT_URL = relpath(ROOT_DIR, joinpath(BUILD_DIR, "generated"))

examples = [
    "basic.jl",
    "io.jl",
    "meshgrid.jl",
    "meshgrid2.jl",
    "wavefieldintegration.jl",
    "integratedreflection.jl",
    "magneticfield.jl"
]

for example in examples
    example_filepath = joinpath(EXAMPLES_DIR, example)
    if get(ENV, "CI", nothing) == "true"
        Literate.markdown(example_filepath, OUTPUT_DIR, documenter=true)
    else
        # local
        Literate.markdown(example_filepath, OUTPUT_DIR, documenter=true,
                          repo_root_url=REPO_ROOT_URL)
    end
end

#==
Organize page hierarchies
==#

example_pages = [
    "generated/basic.md",
    "generated/io.md",
    "generated/meshgrid.md",
    "generated/meshgrid2.md",
    "generated/wavefieldintegration.md",
    "generated/integratedreflection.md",
    "generated/magneticfield.md"
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
