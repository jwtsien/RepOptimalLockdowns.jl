module RepOptimalLockdowns

# export the final function
export run

# the main function
function run()
    include("src/CreateFiles.jl")
    include("src/GSiROptimalPolicy.jl")
    # create the folders
    create_my_folders()
    # get a result
    tag = "Base"
    USF = GSiR_OptimalPolicy(policy="Uniform", type="SF", tag=tag)
end

end # module RepOptimalLockdowns
