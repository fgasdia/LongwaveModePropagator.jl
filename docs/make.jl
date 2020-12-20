"../src/" in LOAD_PATH || push!(LOAD_PATH, "../src/")

using Documenter, LongwaveModePropagator

makedocs(
    sitename="LongwaveModePropagator.jl",
    format=Documenter.HTML(
        prettyurls=get(ENV, "CI", nothing)=="true"
    ),
    modules=[LongwaveModePropagator,LongwaveModePropagator.Fields],
    pages=[
        "Home" => "index.md",
        "Manual" => [
            "examples/demo.md"
        ],
        "Library" => [
            "lib/public.md",
            "lib/internals.md"
        ]
    ],
)
