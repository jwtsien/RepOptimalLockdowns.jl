using Documenter, RepOptimalLockdowns

makedocs(modules = [RepOptimalLockdowns],
         sitename = "RepOptimalLockdowns.jl",
         format = Documenter.HTML()
        )
deploydocs(
        repo = "https://github.com/cwzie/RepOptimalLockdowns",
        )