"../src/" in LOAD_PATH || push!(LOAD_PATH, "../src/")

using Documenter, LongwaveModePropagator

makedocs(
    sitename="Longwave\nMode\nPropagator",
    format=Documenter.HTML(
        prettyurls=get(ENV, "CI", nothing)=="true"
    ),
    modules=[LongwaveModePropagator,Fields],
    pages=[
        "Home" => "index.md",
        "Manual" => [
            "man/demo.md"
        ],
        "Library" => [
            "lib/public.md",
            "lib/internals.md"
        ]
    ]
)
