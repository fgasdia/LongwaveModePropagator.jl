using Documenter
using Literate
using Plots
ENV["GKSwstype"] = "nul"  # to avoid GKS error on docs build ("nul" for headless systems)

using LongwaveModePropagator
using LongwaveModePropagator.Fields

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
    "integratedreflection.jl",
    "wavefieldintegration.jl",
    "magneticfield.jl",
    "interpretinghpbeta.jl",
    "multiplespecies.jl",
    "interpolatingfunctions.jl",
    "ground.jl"
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

example_pages = Any[
    "generated/basic.md",
    "generated/io.md",
    "generated/meshgrid.md",
    "generated/meshgrid2.md",
    "generated/integratedreflection.md",
    "generated/wavefieldintegration.md",
    "generated/magneticfield.md",
    "generated/interpretinghpbeta.md",
    "generated/multiplespecies.md",
    "generated/interpolatingfunctions.md",
    "generated/ground.md"
]

library_pages = Any[
    "lib/public.md",
    "lib/internals.md"
]

pages = Any[
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

makedocs(
    sitename = "LongwaveModePropagator",
    authors = "Forrest Gasdia",
    format = format,
    pages = pages,
    modules = [LongwaveModePropagator]
)

deploydocs(
    repo = "github.com/fgasdia/LongwaveModePropagator.jl.git",
    devbranch = "main",
)
