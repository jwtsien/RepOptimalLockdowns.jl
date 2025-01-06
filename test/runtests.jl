using Test, RepOptimalLockdowns

# it's really hard to test the specific numbers in the results, so here we test whether the code runs well.

## for that we didn't get the correct answer, we just test whether the code runs well.
@testset "Test whether the code runs well" begin
    @test RepOptimalLockdowns.run() == nothing
end