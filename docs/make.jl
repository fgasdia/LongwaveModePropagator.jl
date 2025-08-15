using Documenter
# using Literate
using Plots
ENV["GKSwstype"] = "nul"  # to avoid GKS error on docs build ("nul" for headless systems)

using LongwaveModePropagator
using LongwaveModePropagator.Fields

# push!(LOAD_PATH,"../src/")

#==
Generate examples
==#

# const ROOT_DIR = abspath(@__DIR__, "..")
# const EXAMPLES_DIR = joinpath(ROOT_DIR, "examples")
# const BUILD_DIR = joinpath(ROOT_DIR, "docs", "build")
# const OUTPUT_DIR = joinpath(ROOT_DIR, "docs", "src", "generated")

# repo_root_url is only assigned with CI, so we have to specify it for local builds
# can't use abspath on Windows (no colons) so we calculate relpath
# const REPO_ROOT_URL = relpath(ROOT_DIR, joinpath(BUILD_DIR, "generated"))

# examples = [
#     "basic.jl",
    # "io.jl",
    # "meshgrid.jl",
    # "meshgrid2.jl",
    # "integratedreflection.jl",
    # "wavefieldintegration.jl",
    # "magneticfield.jl",
    # "interpretinghpbeta.jl",
    # "multiplespecies.jl",
    # "interpolatingfunctions.jl",
    # "ground.jl",
    # "fieldcomponents.jl"
# ]

# for example in examples
#     example_filepath = joinpath(EXAMPLES_DIR, example)
#     Literate.markdown(example_filepath, OUTPUT_DIR, documenter=true)
# end

# Documenter.jl can only link to files within docs/src/
# islink("src/examples") || symlink(abspath("../examples"), "src/examples"; dir_target=true)

#==
Build and deploy docs
==#

makedocs(
    sitename = "LongwaveModePropagator",
    authors = "Forrest Gasdia",
    format = Documenter.HTML(
        collapselevel=1
    ),
    pages = [
        "Home" => "index.md",
        # "Examples" => [
        #     "Introduction and defining scenarios" => "basic.md",
            # "File-based I/O" => "io.md",
            # "Mesh grid for mode finding - Part 1" => "meshgrid.md",
            # "Mesh grid for mode finding - Part 2" => "meshgrid2.md",
            # "Solvers for ionosphere reflection coefficient" => "integratedreflection.md",
            # "Wavefield integration" => "wavefieldintegration.md",
            # "Multiple ionospheric species" => "multiplespecies.md",
            # "Magnetic field direction" => "magneticfield.md",
            # "Density and collision frequency as interpolating functions" => "interpolatingfunctions.md",
            # "Multiple field components" => "fieldcomponents.md"
            # ],
        # "Library" => [
            # "public.md",
            # "internals.md"
        # "References" => "references.md"
    ],
    modules = [LongwaveModePropagator]
)

repo = "github.com/fgasdia/LongwaveModePropagator.jl.git"
withenv("GITHUB_REPOSITORY" => repo) do
    deploydocs(; repo, versions=["stable" => "v^", "dev" => "main"])
end
