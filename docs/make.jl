using Documenter
using Plots
ENV["GKSwstype"] = "100"  # to avoid GKS error on docs build

using LongwaveModePropagator

makedocs(
    sitename = "LongwaveModePropagator",
    authors = "Forrest Gasdia",
    format = Documenter.HTML(
        collapselevel=1
    ),
    pages = [
        # "Home" => "index.md",
        "Examples" => [
            # "Introduction and defining scenarios" => "basic.md",
            # "File-based I/O" => "io.md",
            # "Mesh grid for mode finding - Part 1" => "meshgrid.md",
            # "Mesh grid for mode finding - Part 2" => "meshgrid2.md",
            # "Solvers for ionosphere reflection coefficient" => "integratedreflection.md",
            # "Wavefield integration" => "wavefieldintegration.md",
            # "Magnetic field direction" => "magneticfield.md",
            # "Interpreting h′ and β" => "interpretinghpbeta.md",
            # "Multiple ionospheric species" => "multiplespecies.md",
            # "Density and collision frequency as interpolating functions" => "interpolatingfunctions.md",
            "Ground" => "ground.md",
            # "Multiple field components" => "fieldcomponents.md"
            ],
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
