using Documenter, Weave, StatsBase
using ContinuousWavelet

# Convert jmd -> md
foreach(
    filename -> weave(
        joinpath(dirname(pathof(ContinuousWavelet)), "..", "doc", filename),
        doctype = "github",
        fig_path = joinpath("doc", "fig"),
        fig_ext = ".svg",
        out_path = joinpath("doc"),
    ),
    [
        "wavelets.jmd",
        "transform.jmd",
    ]
)

# assemble docs from md sources
makedocs(
    root = joinpath(dirname(pathof(ContinuousWavelet)), "..", "doc"),
    sitename = "ContinuousWavelet",
    pages = [
        "index.md",
        "wavelets.md",
        "transform.md",
    ],
)

# upload
deploydocs(repo = "github.com/hmatuschek/ContinuousWavelet.jl.git", push_preview = true)
