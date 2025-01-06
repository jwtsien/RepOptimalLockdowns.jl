using JuMP, Ipopt, DifferentialEquations, LinearAlgebra, CSV, DataFrames
function GSiR_OptimalPolicy(;theta=.75,             # lockdown effectiveness
                            beta= .134*7*.8,       # contagiousness (infection rate)
                            recovery= 18,          # length of resolution in days (ask for separate code if you want to split out exposure and infection periods)
                            vat= 100,              # mean vaccine arrival time if arrival is stochastic (if you want to model stochastic arrival set vat=1 for arrival in about 1yeasr
                            chim=35,               # lagrange multiplier for fatalities in frontier tracking
                            alpha=2,               # matching/contact techonology in SIR
                            rho=1,                 # default contact bw groups
                            nt=79,                 # time horizon in weeks +1 (python cuts off the last point); 
                            eta=.9,                # percent undetected infected
                            ppi =1,                # proportion able to prove immunity after recovery
                            irate=.01,             # interest rate
                            sh_y =.53,             # share of the young
                            sh_m =.26,             # share of the middle
                            sh_o =.21,             # share of the old
                            l0=.5,                 # initial value for lockdown to start OC  
                            bfr =.001,             # base case fatality ratio for young
                            frm = 10, fro =60,     # relative fatality for m/y .. o/y
                            wr=1, wro=0.26,        # ratio of wages m vs y; ratio of wages o vs y      
                            thetao = .75,          # lockdown effectiveness old
                            rhoo= 1,               # default contact rate of o with y and m
                            kappa=1,               # penalty for high infection rates (death rate increasing with infections)
                            Ri=.005, Ii =.005, 
                            Si=.99,                # initial conditions of recovered, infected, susceptible
                            pwfh =.4,              # proportional loss of productivity under lockdown
                            wfh=0,                 # proportion of work that can be done from home w/loss of productivity
                            policy="Uniform",      # policty set to Uniform, SemiTargeted, or Targeted
                            type="SF",             # type of control: SF (safety first), EF (economy first), or FT (frontier tracking)
                            safetyfirst = .001,    # Fatality Constraint in Safety First Policy
                            economyfirst = .1,     # Econ Cost loss Constraint in Economy First Policy
                            # plot=true,             # Plotting True or False
                            tag="Test",            # Supply a tag for files or figures; can be any character string
                            # etatag=".9",           # supply a tag for testing; can be any character string
                            icubound=false,        # Turn off ICU constraint; set True to turn on
                            htbound =.016,         # ICU implied bound on weighted infection rate
                            # ftag="",               # additional file tag to put in front of model id resut.
                            path="lockdown/",      # path to folder where models and results will be stored/ require subfolders /figs, /summaryres
                            summary=true)         # turn on summary of results being written out in a table
    
    # Policy can be "Uniform", "SemiTargeted", "Targeted"
    # Type can be "SF" (safety first) and "EF" (economy first)
    
    ##################### The orignal code used GEKKO and we use JuMP to replace it #########################
    model = Model(Ipopt.Optimizer)   # initialize the model

    #################### Set solver options ####################
    set_optimizer_attribute(model, "max_iter", 1000)   # usually <300 is enough
    set_optimizer_attribute(model, "tol", 1e-6)        # tolerance for convergence

    ######################## Modeling ##########################
    # Define variables and parameters
    time = LinRange(1, nt, nt)         # time space
    p = zeros(nt); p[end] = 1.0     # mark final time point
    final = p                       # parameter for final time point
    t = 1:nt                  # weeks

    # Economic Parameters
    w = (1/52)/(sh_y * 1 + sh_m * wr + sh_o * wro)  # weekly wage, normalized GDP to 1 per year
    beta = beta          # infection rate
    r = 7 * irate / 365  # interest rate
    nu = 1 / (vat * 52)  # intensity of vaccine arrival
    wr = wr              # wage ratio m/y
    wro = wro            # wage ratio o/y
    ksi = 1 - pwfh       # fraction of work that cannot be done from home
    frm = frm            # fatility ratio m/y
    fro = fro            # fatility ratio o/y
    chi = chim           # emotional cost of death
    safetyfirst = safetyfirst  # safety first constraint
    economyfirst = economyfirst  # economy first constraint

    # Social Interaction Parameters
    theta = theta        # lockdown effectiveness
    thetao = thetao      # lockdown effectiveness - old
    rho = rho            # group interaction
    rhoo = rhoo          # group interaction - old with others
    alpha = alpha        # matching parameter

    # Infection Parameters
    eta = eta            # lockdown effectiveness
    bfr = bfr            # base fatality ratio
    gamma = 7 / recovery # infection exit rate
    gammad = gamma * bfr # base fatality rate for young per week
    gammadm = gammad * frm # base fatality rate for middle
    gammado = gammad * fro # base fatality rate for old
    kappa = kappa        # mortality increasing in the infection rate
    htbound = htbound    # ICU constraint

    # State and Control Variables
    sh_y = sh_y          # share of young
    wh = sh_y * wfh      # share of young who can work from home
    @variable(model, s[1:nt] >= 0, start = Si * sh_y, base_name = "s")                 # susceptible young/total initial population
    @variable(model, i[1:nt] >= 0, start = Ii * sh_y, base_name = "i")                 # share of infected young
    @variable(model, 0 <= l <= 0.7 - wfh, start = l0, base_name = "l")                 # lockdown initial value (control)
    @variable(model, d[1:nt] >= 0, start = 0, base_name =  "d")                 # share of death young
    @variable(model, re[1:nt], start = Ri * sh_y, base_name = "re")                # share of recovered young

    # Similarly define `sh_m`, `sh_o`, and other parameters for the middle and old age groups.
    sh_m = sh_m
    whm = sh_m * wfh
    @variable(model, sm[1:nt] >= 0, start = Si * sh_m, base_name = "sm")
    @constraint(model, sm[1] == Si * sh_m)
    @variable(model, im[1:nt] >= 0, start = Ii * sh_m, base_name = "im")
    @constraint(model, im[1] == Ii * sh_m)
    @variable(model, 0 <= lm <= 0.7 - wfh, start = l0, base_name = "lm")
    @variable(model, dm[1:nt] >= 0, start = 0.01, base_name =  "dm")
    @variable(model, rem[1:nt], start = Ri * sh_m, base_name = "rem")

    sh_o = sh_o
    who = sh_o * wfh
    @variable(model, so[1:nt] >= 0, start = Si * sh_o, base_name = "so")
    # we hope that the first value is the start we set
    @constraint(model, so[1] == Si * sh_o)
    @variable(model, io[1:nt] >= 0, start = Ii * sh_o, base_name = "io")
    @constraint(model, io[1] == Ii * sh_o)
    @variable(model, 0 <= lo <= 1, start = 1, base_name = "lo")
    @variable(model, dol[1:nt] >= 0, start = 0.01, base_name =  "do")
    @variable(model, reo[1:nt], start = Ri * sh_o, base_name = "reo")


    # Cumulative Costs
    @variable(model, cf[1:nt] >= 0, start = 0.01, base_name = "cf")  # cost function
    @variable(model, co[1:nt] >= 0, start = 0.01, base_name = "co")  # cost output
    @variable(model, cd[1:nt] >= 0, start = 0.01, base_name = "cd")  # cost death
    @variable(model, ce[1:nt] >= 0, start = 0.01, base_name = "ce")  # cost emotional

    # Define intermediate variables
    ## Economic value of Life
    @expression(model, ev, (w / r) * (1 - 1 / (1 + r)^(32.5 * 52)))      # pv econ value of life
    @expression(model, evm, (wr * w / r) * (1 - 1 / (1 + r)^(10 * 52)))      # pv econ value of life for the middle
    @expression(model, evo, (wr * w / r) * (1 - 1 / (1 + r)^(2.5 * 52)))    # pv econ value of life for the old
    @expression(model, evr, evm / ev)   # pv econ value ration m/y
    @expression(model, evro, evo / ev)  # pv econ value ration o/y

    ## Cumulative Losses
    @expression(model, death_loss, d[end] + dm[end] + dol[end])  # fatalities
    @expression(model, econ_loss, co[end] + cd[end])        # output loss

    # Rate Functions
    @expression(model, df, exp.(-(r + nu) .* t))  # discount factor
    @expression(model, ht, (1 .* i + frm .* im + fro .* io) ./ (1 .* sh_y + frm .* sh_m + fro .* sh_o))            # hospitalization needs (times constant)
    @expression(model, rd, gammad .* (1 .+ kappa .* ht) .* i)  # death rate for young
    @expression(model, rdm, gammadm .* (1 .+ kappa .* ht) .* im)  # death rate for middle
    @expression(model, rdo, gammado .* (1 .+ kappa .* ht) .* io)  # death rate for old
    @expression(model, rre, i .* gamma .- rd)  # recovery rate for young
    @expression(model, rrem, im .* gamma .- rdm)  # recovery rate for middle
    @expression(model, rreo, io .* gamma .- rdo)  # recovery rate for old

    ## Unisolated infected
    @expression(model, ui, eta .* i)  # un-isolated infected for young
    @expression(model, uim, eta .* im)  # un-isolated infected for middle
    @expression(model, uio, eta .* io)  # un-isolated infected for old

    ## Matching/Network Multipliers
    @expression(model, gmm, (((s .+ ui .+ (1 .- ppi) .* re) .* (1 .- theta .* l .- theta .* wh) .+ ppi .* re) .+ rho .* ((sm .+ uim .+ (1 - ppi) .* rem) .* (1 - theta .* lm .- theta .* whm) .+ ppi .* rem) .+ rhoo .* ((so .+ uio .+ (1 .- ppi) .* reo) .* (1 .- thetao .* lo .- thetao .* who) .+ ppi .* reo)).^(alpha - 2))
    @expression(model, gmmm, (rho .* ((s .+ ui .+ (1 .- ppi) .* re) .* (1 - theta .* l - theta .* wh) + ppi .* re) .+ ((sm .+ uim .+ (1 - ppi) .* rem) .* (1 - theta .* lm .- theta .* whm) .+ ppi .* rem) .+ rhoo .* ((so .+ uio .+ (1 - ppi) .* reo) .* (1 .- thetao .* lo - thetao .* who) .+ ppi .* reo)).^(alpha - 2))
    @expression(model, gmmo, (rhoo .* ((s .+ ui .+ (1 .- ppi) .* re) .* (1 .- theta .* l .- theta .* wh) .+ ppi .* re) .+ rhoo .* ((sm .+ uim .+ (1 .- ppi) .* rem) .* (1 .- theta .* lm .- theta .* whm) .+ ppi .* rem) .+ ((so .+ uio .+ (1 .- ppi) .* reo) .* (1 .- thetao .* lo .- thetao .* who) .+ ppi .* reo)).^(alpha - 2))

    ## Change in Susceptibles y, m, o
    @expression(model, rs, -beta .* s .* (1 - theta .* l - theta .* wh) .* (ui .* (1 .- theta .* l .- theta .* wh) .+ rho .* uim .* (1 .- theta .* lm .- theta .* whm) .+ rhoo .* uio .* (1 .- thetao .* lo .- thetao .* who)) .* gmm)
    @expression(model, rsm, -beta .* sm .* (1 .- theta .* lm .- theta .* whm) .* (rho .* ui .* (1 .- theta .* l .- theta .* wh) .+ uim .* (1 .- theta .* lm .- theta .* whm) .+ rhoo .* uio .* (1 .- thetao .* lo .- thetao .* who)) .* gmmm)
    @expression(model, rso, -beta .* so .* (1 .- thetao .* lo .- thetao .* who) .* (rhoo .* ui .* (1 .- theta .* l .- theta .* wh) .+ rhoo .* uim .* (1 .- theta .* lm .- theta .* whm) .+ uio .* (1 .- thetao .* lo .- thetao .* who)) .* gmmo)

    # Statistical Reproduction Rate or Rt (ratio of new infected cases over new resolved cases)
    # Note only that the population one has structural meaning; others have statistical meaning
    @expression(model, rt, - rs ./ (i .* gamma))       # RR y
    @expression(model, rtm, - rsm ./ (im .* gamma))    # RR m
    @expression(model, rto, - rso ./ (io .* gamma))    # RR o
    @expression(model, rtpop, - (rs + rsm + rso) ./ (i .* gamma + im .* gamma + io .* gamma))               # RR population

    # Define the SIR Equations
    @constraint(model, [t in 1:nt-1], (s[t+1] .- s[t]) .== rs)              # change in susceptible young
    @constraint(model, [t in 1:nt-1], (i[t+1] .- i[t]) .== -rs .- gamma * i[t])  # change in infected young
    @constraint(model, [t in 1:nt-1], (d[t+1] .- d[t]) .== rd)              # change in death young
    @constraint(model, [t in 1:nt-1], (re[t+1] .- re[t]) .== rre)             # change in recovered young
    # then we do the same to middle and old
    @constraint(model, [t in 1:nt-1], (sm[t+1] .- sm[t]) .== rsm)
    @constraint(model, [t in 1:nt-1], (im[t+1] .- im[t]) .== -rsm .- gamma * im[t])
    @constraint(model, [t in 1:nt-1], (dm[t+1] .- dm[t]) .== rdm)
    @constraint(model, [t in 1:nt-1], (rem[t+1].- rem[t]) .== rrem)
    @constraint(model, [t in 1:nt-1], (so[t+1] .- so[t]) .== rso)
    @constraint(model, [t in 1:nt-1], (io[t+1] .- io[t]) .== -rso .- gamma * io[t])
    @constraint(model, [t in 1:nt-1], (dol[t+1] .- dol[t]) .== rdo)
    @constraint(model, [t in 1:nt-1], (reo[t+1] .- reo[t]) .== rreo)

    # Equation for Cost Flows
    @constraint(model, [t in 1:nt-1], (co[t+1] - co[t]) .== ksi .* df .* (w .* (l .* s .+ l .* (1 .- ppi) .* re .+ (1 .- eta .* (1 .- l)) .* i) .+ w .* wr .* (lm .* sm .+ lm .* (1 .- ppi) .* rem .+ (1 .- eta .* (1 .- lm)) .* im) .+ w .* wro .* (lo .* so .+ lo .* (1 .- ppi) .* reo .+ (1 .- eta .* (1 .- lo)) .* io)))  # cost function -- flow
    @constraint(model, [t in 1:nt-1], (cd[t+1] - cd[t]) .== df .* (ev .* rd .+ ev .* evr .* rdm .+ ev .* evro .* rdo))  # death cost -- flow
    @constraint(model, [t in 1:nt-1], (ce[t+1] - ce[t]) .== df .* chi .* (rd .+ rdm .+ rdo))  # emotional cost -- flow
    @constraint(model, [t in 1:nt-1], (cf[t+1] - cf[t]) .== (co[t+1]-co[t])+(cd[t+1]-cd[t])+(ce[t+1]-ce[t]))  # financial cost -- flow

    # Implied ICU constraint
    if icubound == true
        @constraint(model, ht .<= htbound)
    end

    # Uniform and Semi-Targeted Policy Constraints
    if policy == "Uniform"
        @constraint(model, l == lm)  # lockdown y = lockdown m
        @constraint(model, l == lo)  # lockdown o = lockdown y
    elseif policy == "SemiTargeted"
        @constraint(model, l == lm)  # lockdown y = lockdown m
    end

    # Safety-first and Economy-first constraints
    if type == "SF"
        @constraint(model, death_loss .* final .<= safetyfirst)
        @objective(model, Min, econ_loss)
        typefull = "Safety Focused Optimal"
    elseif type == "EF"
        @constraint(model, econ_loss .* final .<= economyfirst)
        @objective(model, Min, death_loss)
        typefull = "Economy Focused Optimal"
    elseif type == "FT" # Frontier Tracking
        @objective(model, Min, cf[end])
        typefull = "Fixed χ Optimal"
    end

    # Solve the optimization problem
    optimize!(model)

    ################### GETTING THE RESULTS ####################
    # Access the results
    # optimal_values = value.(model[:s])  # for example, to get the optimal values of the susceptible young
    
    # Create result container (assuming `result_container` is a defined structure in your code)
    result = Dict()
    result["tag"] = tag
    result["policy"] = policy
    result["theta"] = theta[end]
    result["chi"] = chim
    result["rho"] = rho[end]
    result["alpha"] = alpha[end]
    result["eta"] = eta[end]
    result["nt"] = nt
    result["econ_loss"] = value.(econ_loss)[end]
    result["death_loss"] = value.(death_loss)[end]
    result["output_loss"] = value.(co)[end]
    result["futoutput_loss"] = value.(cd)[end]
    result["fr"] = value.(d)[end] / value.(sh_y)[end]
    result["frm"] = value.(dm)[end] / value.(sh_m)[end]
    result["fro"] = value.(dol)[end] / value.(sh_o)[end]
    result["d"] = value.(d)[end]
    result["dtraj"] = death_loss
    result["dm"] = value.(dm)[end]
    result["do"] = value.(dol)[end]
    result["rt"] = rt[2:nt]
    result["rtm"] = rtm[2:nt]
    result["rto"] = rto[2:nt]
    result["lost_l"] = sum(value.(l)) / nt  # lost proportion of y labor
    result["lost_lm"] = sum(value.(lm)) / nt  # lost proportion of m labor
    result["lost_lo"] = sum(value.(lo)) / nt  # lost proportion of o labor

    sumtable = [
        result["tag"], result["policy"], result["theta"], result["chi"], result["eta"], result["rho"], 
        result["alpha"], result["nt"], result["econ_loss"], result["death_loss"], result["output_loss"], result["futoutput_loss"], 
        result["fr"], result["frm"], result["fro"], result["d"], result["dm"], result["do"], 
        result["lost_l"], result["lost_lm"], result["lost_lo"]
    ]

    columns = [
        "tag", "Policy", "theta", "chi", "eta", "rho",
        "alpha", "horison (weeks)", "Economic Loss", "Fatalities",
        "Output Loss", "Future Output Loss", "FR_y", "FR_m", "FR_o",
        "d_y", "d_m", "d_o",
        "sum(ly)/T", "sum(lm)/T", "sum(lo)/T"
    ]

    summaryres_dir = path * "summaryres/"
    # save the result into a .csv file
    if summary == true
        fname = summaryres_dir * "summary.csv"
        if !isfile(fname)
            CSV.write(fname, DataFrame(permutedims(sumtable), columns))
        else
            CSV.write(fname, DataFrame(permutedims(sumtable), :auto), append = true)
        end
    end
    return result
end