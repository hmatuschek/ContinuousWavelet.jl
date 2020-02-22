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
        "coherence.jmd",
        "utils.jmd"
    ]
)

# assemble docs from md sources
makedocs(
    root = joinpath(dirname(pathof(ContinuousWavelet)), "..", "docs"),
    sitename = "ContinuousWavelet",
    pages = [
        "index.md",
        "transform.md",
        "wavelets.md",
        "coherence.md",
        "utils.md"
    ],
)

# upload
deploydocs(repo = "github.com/hmatuschek/ContinuousWavelet.jl.git", push_preview = true)
