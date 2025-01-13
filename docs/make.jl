# Use
#
#     DOCUMENTER_DEBUG=true julia --color=yes make.jl local [nonstrict] [fixdoctests]
#
# for local builds.

using Documenter
using MCBench

makedocs(
    sitename = "MCBench",
    modules = [MCBench],
    format = Documenter.HTML(
        prettyurls = !("local" in ARGS),
        canonical = "https://github.com/tudo-physik-e4/MCBench"
    ),
    pages = [
        "Home" => "index.md",
        "Tutorial" => "tutorial.md",
        "Test cases & Metrics" => "testcases.md",
        "API" => "api.md",
        "LICENSE" => "LICENSE.md",
    ],
    doctest = ("fixdoctests" in ARGS) ? :fix : true,
    linkcheck = ("linkcheck" in ARGS),
    warnonly = !("nonstrict" in ARGS),
)

deploydocs(
    repo = "github.com/tudo-physik-e4/MCBench.git",
    forcepush = true,
    push_preview = true,
    devbranch = "main"
)
