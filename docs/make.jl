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
        "synthesis.jmd",
        "coherence.jmd",
        "utils.jmd",
        "examples.jmd"
    ]
)

# assemble docs from md sources
makedocs(
    root = joinpath(dirname(pathof(ContinuousWavelet)), "..", "docs"),
    sitename = "ContinuousWavelet",
    pages = [
        "index.md",
        "transform.md",
        "synthesis.md",
        "wavelets.md",
        "coherence.md",
        "utils.md",
        "examples.md"
    ],
)

# upload
deploydocs(repo = "github.com/hmatuschek/ContinuousWavelet.jl.git", push_preview = true)
