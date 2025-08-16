using Documenter
using Plots
ENV["GKSwstype"] = "100"  # to avoid GKS error on docs build

using LongwaveModePropagator
using LongwaveModePropagator.Fields

# push!(LOAD_PATH,"../src/")

#==
Generate examples
==#

# repo_root_url is only assigned with CI, so we have to specify it for local builds
# can't use abspath on Windows (no colons) so we calculate relpath
# const REPO_ROOT_URL = relpath(ROOT_DIR, joinpath(BUILD_DIR, "generated"))

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
            # "Magnetic field direction" => "magneticfield.md",
            # "Interpreting h′ and β" => "interpretinghpbeta.md",
            # "Multiple ionospheric species" => "multiplespecies.md",
            # "Density and collision frequency as interpolating functions" => "interpolatingfunctions.md",
            # "Ground" => "ground.md",
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
