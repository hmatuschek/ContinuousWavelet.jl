using Documenter, Weave, StatsBase
using ContinuousWavelet

# Convert jmd -> md
foreach(
    filename -> weave(
        joinpath(dirname(pathof(ContinuousWavelet)), "..", "docs", "jmd", filename),
        doctype = "github",
        fig_path = joinpath("docs", "fig"),
        fig_ext = ".svg",
        out_path = joinpath("docs", "src"),
    ),
    [
        "wavelets.jmd",
        "transform.jmd",
    ]
)

# assemble docs from md sources
makedocs(
    root = joinpath(dirname(pathof(ContinuousWavelet)), "..", "docs"),
    sitename = "ContinuousWavelet",
    pages = [
        "index.md",
        "wavelets.md",
        "transform.md",
    ],
)

# upload
deploydocs(repo = "github.com/hmatuschek/ContinuousWavelet.jl.git", push_preview = true)
