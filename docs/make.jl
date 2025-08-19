using Documenter
using LongwaveModePropagator

ENV["GKSwstype"] = "100"  # to avoid GKS error on docs build

makedocs(
    sitename = "LongwaveModePropagator",
    authors = "Forrest Gasdia",
    format = Documenter.HTML(
        size_threshold_ignore = ["internals.md"]
    ),
    pages = [
        "Home" => "index.md",
        "Examples" => [
            "Introduction and defining scenarios" => "basic.md",
            "Interpreting h′ and β" => "interpretinghpbeta.md",
            "Magnetic field direction" => "magneticfield.md",
            "Ground" => "ground.md",
            "Multiple ionospheric species" => "multiplespecies.md",
            "Mesh grid for mode finding - Part 1" => "meshgrid.md",
            "Mesh grid for mode finding - Part 2" => "meshgrid2.md",
            "Solvers for ionosphere reflection coefficient" => "integratedreflection.md",
            "Wavefield integration" => "wavefieldintegration.md",
            "Density and collision frequency as interpolating functions" => "interpolatingfunctions.md",
            "File-based I/O" => "io.md",
            "Multiple field components" => "fieldcomponents.md"
            ],
        "Library" => [
            "Public API" => "public.md",
            "Internals" => "internals.md"
        ],
        "References" => "references.md"
    ],
    modules = [LongwaveModePropagator],
    # pagesonly = true
)

repo = "github.com/fgasdia/LongwaveModePropagator.jl.git"
withenv("GITHUB_REPOSITORY" => repo) do
    deploydocs(; repo, versions=["stable" => "v^", "dev" => "main"])
end
