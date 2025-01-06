module RepOptimalLockdowns

# export the final function
export run

# the main function
function run()
    include("src/create_files.jl")
    include("src/GSiR_optimal_policy.jl")
    include("src/frontier_GSiR.jl")
    # create the folders
    create_my_directory()
    # get results from the optimal policy function
    tag = "Base"
    USF = GSiR_OptimalPolicy(policy="Uniform", type="SF",  tag=tag)
    UFT = GSiR_OptimalPolicy(policy="Uniform", type="FT", chim=35, tag=tag)
    UEF = GSiR_OptimalPolicy(policy="Uniform", type="EF", l0 = 0, tag=tag)
    STSF = GSiR_OptimalPolicy(policy="SemiTargeted", type="SF", tag=tag)
    STEF = GSiR_OptimalPolicy(policy="SemiTargeted", type="EF",l0= 0.1, tag=tag)
    STFT = GSiR_OptimalPolicy(policy="SemiTargeted", type="FT",  chim=35, tag=tag)
    
    # the other way we tried
    chimseq = [0.01, 10, 40, 60, 80]  # Testing with a sequence of chim values
    frontier_GSiR(chimseq)
end

end # module RepOptimalLockdowns
