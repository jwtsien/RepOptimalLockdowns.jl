using Test, RepOptimalLockdowns

# it's really hard to test the specific numbers in the results, so here we test whether the code runs well.

## we would like to test whether all the results are exported correctly
@test ispath("/lockdown/figs")
@test ispath("/lockdown/models")
@test ispath("/lockdown/results")
@test ispath("/lockdown/summaryres")

@test isfile("/lockdown/summaryres/summary.csv")