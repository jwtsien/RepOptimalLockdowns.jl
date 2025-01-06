module RepOptimalLockdowns

# export the final function
export run

# the main function
function run()
    include("src/CreateFiles.jl")
    include("src/GSiROptimalPolicy.jl")
    # create the folders
    create_my_folders()
    # get results from the optimal policy function
    tag = "Base"
    USF = GSiR_OptimalPolicy(policy="Uniform", type="SF",  tag=tag)
    UFT = GSiR_OptimalPolicy(policy="Uniform", type="FT", chim=35, tag=tag)
    UEF = GSiR_OptimalPolicy(policy="Uniform", type="EF", l0= 0, tag=tag)
    STSF = GSiR_OptimalPolicy(policy="SemiTargeted", type="SF", tag=tag)
    STEF = GSiR_OptimalPolicy(policy="SemiTargeted", type="EF",l0= 0.1, tag=tag)
    STFT = GSiR_OptimalPolicy(policy="SemiTargeted", type="FT",  chim=35, tag=tag)
    
    # the other way we tried
end

end # module RepOptimalLockdowns
