using DifferentialEquations, JuMP, Ipopt, Plots

# Initial Conditions
Si = 0.99
Ii = 0.005
Ri = 0.005

# Parameters
const beta = 0.134 * 7 * 0.8  # Infection rate
const recovery = 18.0  # Recovery period (days)
const vat = 100.0  # Vaccine arrival time (weeks)
const contact_alpha = 2.0  # Contact technology
const rho = 1.0  # Interaction between groups
const rhoo = 1.0  # Interaction between middle-aged and old groups
const nt = 78  # Time horizon (weeks)
const l0 = 0.7  # Initial lockdown
const irate = 0.01  # Interest rate
const frm = 10.0  # Fatality ratio (middle to young)
const fro = 60.0  # Fatality ratio (old to young)
const wr = 1.0  # Wage ratio (middle/young)
const wro = 0.26  # Wage ratio (old/young)
const sh_y = 0.53  # Share of young
const sh_m = 0.26  # Share of middle-aged
const sh_o = 0.21  # Share of old
const kappa = 1.0  # Penalty for infection rates
const gamma = 7.0 / recovery  # Exit rate per week
const htbound = 0.016  # ICU capacity bound

# Calculate Annual GDP
const gdp = (sh_y * wr + sh_m * wr + sh_o * wro)
# Constants for Population
const total_pop = Si * (sh_y + sh_m + sh_o)


# Ensure Subfolder Exists
save_path = "lockdown/figs"
mkpath(save_path)  # Creates the folder if it doesn't exist

function sir!(du, u, p, t)
    s = u[1:3]  # Susceptible states for each group
    i = u[4:6]  # Infected states for each group
    r = u[7:9]  # Recovered states for each group

    beta, theta_y, theta_m, theta_o, gamma, rho, rhoo, frm, fro = p

    du[1] = -beta * s[1] * (i[1] * (1 - theta_y) + rho * i[2] * (1 - theta_m) + rhoo * i[3] * (1 - theta_o))
    du[2] = -beta * s[2] * (i[1] * (1 - theta_y) + rho * i[2] * (1 - theta_m) + rhoo * i[3] * (1 - theta_o))
    du[3] = -beta * s[3] * (i[1] * (1 - theta_y) + rho * i[2] * (1 - theta_m) + rhoo * i[3] * (1 - theta_o))

    du[4] = beta * s[1] * (i[1] * (1 - theta_y) + rho * i[2] * (1 - theta_m) + rhoo * i[3] * (1 - theta_o)) - gamma * i[1]
    du[5] = beta * s[2] * (i[1] * (1 - theta_y) + rho * i[2] * (1 - theta_m) + rhoo * i[3] * (1 - theta_o)) - gamma * i[2]
    du[6] = beta * s[3] * (i[1] * (1 - theta_y) + rho * i[2] * (1 - theta_m) + rhoo * i[3] * (1 - theta_o)) - gamma * i[3]

    du[7] = gamma * i[1]
    du[8] = gamma * i[2]
    du[9] = gamma * i[3]
end

#objective function
function objective(theta_y, theta_m, theta_o, s, i, r, nt, chim)
    econ_loss = 0.0
    death_loss = 0.0
    for t in 1:nt
        econ_loss += ((1 - theta_y[t]) * (s[1][t] + r[1][t]) +
                      (1 - theta_m[t]) * (s[2][t] + r[2][t]) * wr +
                      (1 - theta_o[t]) * (s[3][t] + r[3][t]) * wro) / gdp
        death_loss += (i[1][t] * gamma * chim +
                       i[2][t] * gamma * frm * chim +
                       i[3][t] * gamma * fro * chim) / gdp
    end
    return econ_loss + death_loss
end


# Simulate Function (Global Scope)
function simulate(theta_y, theta_m, theta_o)
    s = [[Si * sh_y], [Si * sh_m], [Si * sh_o]]
    i = [[Ii * sh_y], [Ii * sh_m], [Ii * sh_o]]
    r = [[Ri * sh_y], [Ri * sh_m], [Ri * sh_o]]

    for t in 1:nt
        u0 = [s[1][end], s[2][end], s[3][end], i[1][end], i[2][end], i[3][end], r[1][end], r[2][end], r[3][end]]
        p = [beta, theta_y[t], theta_m[t], theta_o[t], gamma, rho, rhoo, frm, fro]
        prob = ODEProblem(sir!, u0, (0.0, 1.0), p)
        sol = solve(prob, Tsit5(), saveat=1.0)
        for k in 1:3
            push!(s[k], sol[k, end])
            push!(i[k], sol[k+3, end])
            push!(r[k], sol[k+6, end])
        end
    end
    return s, i, r
end


# Optimize Lockdown Policy
function optimize_policy(chim::Float64, policy::String="Uniform")
    m = Model(Ipopt.Optimizer)

    @variable(m, 0 <= theta_y[1:nt] <= 0.7)
    @variable(m, 0 <= theta_m[1:nt] <= 0.7)
    @variable(m, 0 <= theta_o[1:nt] <= 1.0)

    if policy == "Uniform"
        @constraint(m, theta_y .== theta_m)
        @constraint(m, theta_m .== theta_o)
    elseif policy == "SemiTargeted"
        @constraint(m, theta_y .== theta_m)
    end

    @NLobjective(m, Min,
        sum((1 - theta_y[t]) * (Si * sh_y + Ri * sh_y) +
            (1 - theta_m[t]) * (Si * sh_m + Ri * sh_m) * wr +
            (1 - theta_o[t]) * (Si * sh_o + Ri * sh_o) * wro +
            chim * (Ii * sh_y * gamma + Ii * sh_m * gamma * frm + Ii * sh_o * gamma * fro)
            for t in 1:nt)
    )

    optimize!(m)
    return (value.(theta_y), value.(theta_m), value.(theta_o))
end

# Frontier Function
function frontier_GSiR(chimseq::Vector{Float64})
    econ_losses_uniform = Float64[]
    death_losses_uniform = Float64[]
    econ_losses_semi = Float64[]
    death_losses_semi = Float64[]
    econ_losses_targeted = Float64[]
    death_losses_targeted = Float64[]

    for chim in chimseq
        # Uniform Policy
        theta_y, theta_m, theta_o = optimize_policy(chim, "Uniform")
        s, i, r = simulate(theta_y, theta_m, theta_o)
        loss = objective(theta_y, theta_m, theta_o, s, i, r, nt, chim)
        fatalities = (sum(i[1]) * gamma + sum(i[2]) * gamma * frm + sum(i[3]) * gamma * fro) / total_pop
        push!(econ_losses_uniform, loss)
        push!(death_losses_uniform, fatalities)

        # Semi-Targeted Policy
        theta_y, theta_m, theta_o = optimize_policy(chim, "SemiTargeted")
        s, i, r = simulate(theta_y, theta_m, theta_o)
        loss = objective(theta_y, theta_m, theta_o, s, i, r, nt, chim)
        fatalities = (sum(i[1]) * gamma + sum(i[2]) * gamma * frm + sum(i[3]) * gamma * fro) / total_pop
        push!(econ_losses_semi, loss)
        push!(death_losses_semi, fatalities)

        # Targeted Policy
        theta_y, theta_m, theta_o = optimize_policy(chim, "Targeted")
        s, i, r = simulate(theta_y, theta_m, theta_o)
        loss = objective(theta_y, theta_m, theta_o, s, i, r, nt, chim)
        fatalities = (sum(i[1]) * gamma + sum(i[2]) * gamma * frm + sum(i[3]) * gamma * fro) / total_pop
        push!(econ_losses_targeted, loss)
        push!(death_losses_targeted, fatalities)
    end

    # Plot Fatalities vs Economic Losses
    key_indices = [1, 3, 5]  # Markers at selected points

    plot(death_losses_uniform, econ_losses_uniform,
         label="Uniform policy", lw=2, line=:dashdot,
         color=:red)

    plot!(death_losses_semi, econ_losses_semi,
          label="Semi-targeted policy", lw=2, line=:dash,
          color=:green)

    plot!(death_losses_targeted, econ_losses_targeted,
          label="Targeted policy", lw=2, line=:dot,
          color=:blue)

    xlabel!("Fatalities as fraction of adult population")
    ylabel!("PDV of economic loss in fractions of GDP")

    plot!(legend=:topright)

    savefig(joinpath(save_path, "Figure3_Replication.pdf"))
end


###problem: wrong scale - indicating that the calculation of share is wrong, and no labeling of different policy objectives on the frontier 


# Keeping the existing three curves, I also want the graph to mark: 1. the trade-off achieved by the optimal policy that is uniform and under type == EF; 2. the trade-off achieved by uniform policy and under type