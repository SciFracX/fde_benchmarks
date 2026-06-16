"""
    Bio-Ethanol Fermentation: Parameter Fitting and Model Comparison
    
    This script performs a comprehensive analysis of bio-ethanol production through
    microbial fermentation by comparing integer-order ODE and fractional-order ODE models.
    
    Key objectives:
    1. Load experimental fermentation data (biomass, substrate, product concentrations)
    2. Fit both classical integer-order and fractional-order differential equation models
    3. Compare model predictions with experimental measurements
    4. Perform train-test split validation to assess generalization
    5. Apply information-theoretic model selection criteria (AIC/BIC)
    
    Models:
    - Integer-order ODE: Classical Lotka-Volterra predator-prey dynamics
    - Integer-order ODE: Advanced Moser-Luong model with growth kinetics
    - Fractional-order ODE: Lotka-Volterra with fractional derivatives for memory effects
    
    The analysis demonstrates whether fractional-order models better capture
    the memory-dependent and non-local temporal dynamics of fermentation processes.
"""

# ========================================================================
# SECTION 0: Import Required Libraries
# ========================================================================

using CSV          # CSV file reading and writing
using Optim        # Optimization algorithms (BFGS, Fminbox)
using DataFrames   # Tabular data structures
using OrdinaryDiffEq  # ODE solver suite
using FractionalDiffEq  # Fractional-order ODE solvers
using StatsBase    # Statistical functions (RMSD, quantiles)
using CairoMakie   # High-quality plotting and visualization
using Printf       # Formatted printing

# ========================================================================
# Fermentation System Parameters and ODE Model Definition
# ========================================================================
# These parameters define a Lotka-Volterra model applied to fermentation:
#   - Predator (y) represents substrate concentration
#   - Prey (x) represents biomass concentration
# The parameters have been empirically determined from literature:

kc = 0.0041        # Transmission coefficient from infected individuals
                   # (In fermentation context: growth rate coefficient)
km = 2.3875e-14    # Relative transmissibility of hospitalized patients
                   # (In fermentation: very small mortality/degradation rate)
ks = 0.0585        # Transmission coefficient due to super-spreaders
                   # (In fermentation: substrate consumption rate)
kp = 0.0156        # Rate at which exposed become infectious
                   # (In fermentation: product formation rate)

function lotka_volterra(du, u, p, t)
    """
        lotka_volterra(du, u, p, t)
    
    Classical Lotka-Volterra predator-prey dynamics applied to fermentation.
    State variables (u):
        u[1]: Biomass concentration (g/L) - "prey"
        u[2]: Substrate concentration (g/L) - "predator"
        u[3]: Product concentration (g/L) - cumulative product
    
    Parameters (p):
        p[1] = kc: Growth interaction coefficient
        p[2] = km: Death rate of biomass
        p[3] = ks: Consumption rate of substrate
        p[4] = kp: Product formation rate
    
    Dynamics:
        d(biomass)/dt       = kc * biomass * substrate - km * biomass
        d(substrate)/dt     = -ks * biomass * substrate
        d(product)/dt       = kp * biomass * substrate
    """
    kc, km, ks, kp = p
    du[1] = kc * u[1] * u[2] - km * u[1]  # Biomass: growth (interaction) - death
    du[2] = -ks * u[1] * u[2]              # Substrate: consumption (both organisms)
    du[3] = kp * u[1] * u[2]               # Product: formation from substrate
end
# ========================================================================
# Experimental Data Loading and Preparation
# ========================================================================
# Load fermentation kinetics data from CSV files
# Each file contains time-series measurements of concentrations

biomass = DataFrame(CSV.File("/Users/quqingyu/SciFracX/paper/applications/bio_ethanol/biomass.csv"))
    # Column 1: Time (hours)
    # Column 2: Biomass concentration (g/L)

product = DataFrame(CSV.File("/Users/quqingyu/SciFracX/paper/applications/bio_ethanol/product.csv"))
    # Column 1: Time (hours)
    # Column 2: Product concentration (g/L)

substrate = DataFrame(CSV.File("/Users/quqingyu/SciFracX/paper/applications/bio_ethanol/substrate.csv"))
    # Column 1: Time (hours)
    # Column 2: Substrate concentration (g/L)

# Assemble measurement data matrix: rows = variables, columns = time points
# Order: [biomass; substrate; product]
data = Matrix(hcat(biomass[!, 2], substrate[!, 2], product[!, 2])')
    # data[1, :] = biomass at measurement times
    # data[2, :] = substrate at measurement times
    # data[3, :] = product at measurement times

# Measurement time points (hours)
t = [0, 2, 4, 6, 8, 10, 12, 16, 28, 44, 48]

# Pack model parameters
par = [kc, km, ks, kp]  # Parameter vector for ODE model

# Initial conditions for fermentation state
X0 = [0.5, 90.0, 2.5]   # [biomass (g/L), substrate (g/L), product (g/L)]

# Temporal integration domain
tspan = (0, 50)  # Simulate from 0 to 50 hours
# ========================================================================
# Loss Functions for Model Fitting
# ========================================================================
# Three loss functions correspond to three different models:
# 1. Integer-order ODE (Lotka-Volterra)
# 2. Fractional-order ODE (Lotka-Volterra with memory)
# 3. Advanced integer-order model (Moser-Luong kinetics)

function loss_1(b)
    """
        loss_1(b) → Float64
    
    Loss function for integer-order ODE model.
    Solves the classical Lotka-Volterra system and compares predictions
    with experimental measurements.
    
    Arguments:
        b::Vector: Parameter vector [kc, km, ks, kp] to optimize
    
    Returns:
        Normalized root-mean-square error (NRMSE) between model and data
    """
    par = copy(b)  # Copy to avoid modifying original
    prob = ODEProblem(lotka_volterra, X0, tspan, par)
    sol = solve(prob, Tsit5(), saveat=t)  # Tsit5: 5th-order Tsitouras integrator
    appX = reduce(hcat, sol.u)  # Stack solution vectors into matrix
    rmsd(data, appX; normalize=true)  # Normalized RMSE
end

function loss_2(b)
    """
        loss_2(b) → Float64
    
    Loss function for fractional-order ODE model.
    Solves the fractional Lotka-Volterra system where derivatives are
    replaced with fractional derivatives of orders α₁, α₂, α₃.
    
    Arguments:
        b::Vector: [α₁, α₂, α₃, kc, km, ks, kp]
                   First 3 elements: fractional derivative orders
                   Last 4 elements: model parameters
    
    Returns:
        Normalized root-mean-square error between fractional model and data
    
    Note:
        The fractional model adds memory/non-local effects to fermentation dynamics,
        potentially capturing substrate aging, biofilm formation, or metabolic memory.
    """
    par = b[4:end]    # Extract parameters [kc, km, ks, kp]
    order = b[1:3]    # Extract fractional orders [α₁, α₂, α₃]
    prob = FODEProblem(lotka_volterra, order, X0, tspan, par)
    osol = solve(prob, PITrap(), dt = 1)  # Predictor-Integrator Trap-rule method
    # Extract solution at measurement time indices
    app = [osol[1], osol[3], osol[5], osol[7], osol[9], osol[11], osol[13], osol[17], osol[29], osol[45], osol[49]]
    appX = reduce(hcat, app)
    rmsd(data, appX; normalize=true)
end

# Advanced fermentation kinetics model incorporating growth inhibition
# and product inhibition effects

function moser_luong_lotka_volterra(du, u, p, t)
    """
        moser_luong_lotka_volterra(du, u, p, t)
    
    Advanced fermentation kinetics model combining Monod growth kinetics,
    product inhibition, and substrate uptake dynamics.
    
    State variables (u):
        u[1]: Biomass concentration (g/L)
        u[2]: Substrate concentration (g/L)
        u[3]: Product concentration (g/L)
    
    Parameters (p):
        p[1] = um: Maximum specific growth rate (1/h)
        p[2] = ks: Substrate saturation constant (g/L)
        p[3] = kp: Product inhibition constant (g/L)
        p[4] = ms: Maintenance coefficient for substrate
        p[5] = yxs: Biomass yield on substrate (g biomass / g substrate)
        p[6] = yps: Product yield on substrate (g product / g substrate)
        p[7] = alpha: Basal product formation rate
        p[8] = beta: Growth-associated product formation rate
    
    Model Features:
        - Monod kinetics (n=1) for substrate limitation
        - Product inhibition (m=9, supresses growth at high product)
        - Growth-phase dependence via logistic term: 1/(1+exp(-t+1))
        - Maintenance coefficient for basal metabolic needs
    """
    um, ks, kp, ms, yxs, yps, alpha, beta = p
    n = 1    # Monod order (first-order substrate limitation)
    m = 9    # Product inhibition order (steep inhibition)
    
    # Specific growth rate with substrate limitation (Monod) and product inhibition
    mu = um * u[2]^n / (ks + u[2]^n) * (1 - (u[3] / kp)^m) / (1 + exp(-t + 1))
    
    # Biomass accumulation: growth minus maintenance
    du[1] = mu * u[1]
    
    # Substrate consumption: for growth and maintenance
    du[2] = -(1 / yxs + 1 / yps) * mu * u[1] - ms * u[1]
    
    # Product formation: growth-associated and non-growth associated components
    du[3] = alpha * yps * mu * u[1] + beta * u[1]
end

function loss_3(b)
    """
        loss_3(b) → Float64
    
    Loss function for advanced Moser-Luong kinetics model.
    Fit this model to compare with simpler Lotka-Volterra approaches.
    
    Arguments:
        b::Vector: Parameter vector for Moser-Luong model
    
    Returns:
        Normalized root-mean-square error
    """
    par = copy(b)
    prob = ODEProblem(moser_luong_lotka_volterra, X0, tspan, par)
    sol = solve(prob, Tsit5(), saveat=t)
    appX = reduce(hcat, sol.u)
    rmsd(data, appX; normalize=true)
end

# ========================================================================
# Parameter Optimization for Model 1 (Integer-Order ODE)
# ========================================================================
# Use box-constrained optimization (Fminbox) with BFGS descent algorithm
# to fit the classical Lotka-Volterra model to experimental data.

# Define parameter search bounds
p_lo_1 = [0.003, 0.0, 0.05, 0.015]      # Lower bounds for [kc, km, ks, kp]
p_up_1 = [0.005, 0.0005, 0.06, 0.016]   # Upper bounds for [kc, km, ks, kp]
p_vec_1 = [0.004, 0.0001, 0.055, 0.0155]  # Initial guess for optimization

# Perform optimization using L-BFGS with box constraints
# BFGS: Broyden–Fletcher–Goldfarb–Shanno quasi-Newton algorithm
# Excellent for smooth, moderately-sized nonlinear problems
Res1 = optimize(
    loss_1,                    # Objective function to minimize
    p_lo_1,                    # Lower bounds
    p_up_1,                    # Upper bounds
    p_vec_1,                   # Initial point
    Fminbox(BFGS()),          # Box constraints with BFGS algorithm
    Optim.Options(
        outer_iterations = 10,   # Number of outer box-constrained iterations
        iterations = 10000,      # Maximum inner BFGS iterations
        show_trace = true,       # Print iteration details
        show_every = 1           # Show trace every iteration
    )
)

# Extract optimized parameters
p1 = vcat(Optim.minimizer(Res1))  # Convert to vector format



# ========================================================================
# Parameter Optimization for Model 2 (Fractional-Order ODE)
# ========================================================================
# Fit a fractional-order variant of the Lotka-Volterra model.
# Parameters include both the fractional derivative orders and
# the original fermentation kinetics parameters.

# Define parameter search bounds
# Parameters: [α₁, α₂, α₃, kc, km, ks, kp] (7 total)
p_lo_2 = [0.7, 0.8, 0.9, 0.003, 0.0, 0.05, 0.015]   # Lower bounds
p_up_2 = [0.8, 0.9, 1.0, 0.005, 0.0005, 0.06, 0.016] # Upper bounds
p_vec_2 = [0.75, 0.85, 0.95, 0.004, 0.0001, 0.055, 0.0155]  # Initial guess

# Optimize fractional-order model
Res2 = optimize(
    loss_2,
    p_lo_2,
    p_up_2,
    p_vec_2,
    Fminbox(BFGS()),
    Optim.Options(
        outer_iterations = 10,
        iterations = 10000,
        show_trace = true,
        show_every = 1
    )
)

p2 = vcat(Optim.minimizer(Res2))





# ========================================================================
# Parameter Optimization for Model 3 (Moser-Luong Kinetics)
# ========================================================================
# Fit the advanced fermentation kinetics model incorporating
# Monod growth saturation and product inhibition.

# Define parameter search bounds
# Parameters: [um, ks, kp, ms, yxs, yps, alpha, beta] (8 total)
p_lo_3 = [0.3, 15, 120, 0.0, 0.0, 0.4, 20, 0.0]      # Lower bounds
p_up_3 = [0.4, 25, 140, 0.1, 0.1, 0.5, 25, 0.05]     # Upper bounds
p_vec_3 = [0.35, 19, 125, 0.01, 0.01, 0.45, 24, 0.03]  # Initial guess

# Optimize Moser-Luong model
Res3 = optimize(
    loss_3,
    p_lo_3,
    p_up_3,
    p_vec_3,
    Fminbox(BFGS()),
    Optim.Options(
        outer_iterations = 10,
        iterations = 10000,
        show_trace = true,
        show_every = 1
    )
)

p3 = vcat(Optim.minimizer(Res3))


# ========================================================================
# Model Solutions and Performance Comparison
# ========================================================================
# Solve both ODE and FODE models using optimized parameters
# and compute prediction errors.

using Plots  # Additional plotting library for extended visualization

# Solve optimized integer-order ODE model using Tsit5
aprob1 = ODEProblem(lotka_volterra, X0, tspan, p1)
ode_sol = solve(aprob1, Tsit5())
    # ode_sol.t: solution time points
    # ode_sol.u: solution values at each time point

# Solve optimized fractional-order ODE model using PITrap
aprob2 = FODEProblem(lotka_volterra, p2[1:3], X0, tspan, p2[4:end])
fde_sol = solve(aprob2, PITrap(), dt = 0.01)
    # fde_sol[i]: solution value at time t=i

# Evaluate ODE predictions at measurement times
sol1_err = solve(aprob1, Tsit5(), saveat=t)

# Evaluate FODE predictions at measurement times
osol = solve(aprob2, PITrap(), dt = 1)  # dt=1 for exact time points
sol2_err = [osol[1], osol[3], osol[5], osol[7], osol[9], osol[11], osol[13], osol[17], osol[29], osol[45], osol[49]]

# Calculate normalized RMSE for each model (using measurement times only)
err1 = rmsd(data, reduce(hcat, sol1_err.u); normalize=true)  # ODE error
err2 = rmsd(data, reduce(hcat, sol2_err); normalize=true)    # FODE error

# ========================================================================
# Visualization of Model Fits and Comparison
# ========================================================================
# Create a comprehensive four-panel comparison figure showing:
#   - Biomass, substrate, and product concentrations
#   - Experimental measurements (scatter points)
#   - Model predictions (continuous vs dashed lines for ODE vs FODE)
#   - Performance comparison via RMSE bar chart (bottom-right)

fig = Figure(size = (700, 490))
axb = Axis(fig[1, 1], xlabel = "Time (h)", ylabel = "Concentration of biomass (g/L)")

axs = Axis(fig[1, 2], xlabel = "Time (h)", ylabel = "Concentration of substrate (g/L)")

axp = Axis(fig[2, 1], xlabel = "Time (h)", ylabel = "Concentration of product (g/L)")

ax_bar = Axis(fig[2, 2], xticks = (1:2, ["Integer order model", "Fractional order model"]), title = "Normalized RMSE")
ax_bar.xlabelsize=3
scb = CairoMakie.scatter!(axb, biomass[!, 1], biomass[!, 2], color = :blue, label = "Experimental Biomass")
scp = CairoMakie.scatter!(axp, product[!, 1], product[!, 2], color = :green, label = "Experimental Product")
scs = CairoMakie.scatter!(axs, substrate[!, 1], substrate[!, 2], color = :red, label = "Experimental Substrate")
fde_lineb=lines!(axb, fde_sol.t, fde_sol[1, :], color = :blue, label = "Biomass (model)")
fde_lines=lines!(axs, fde_sol.t, fde_sol[2, :], color = :red, label = "Substrate (model)")
fde_linep=lines!(axp, fde_sol.t, fde_sol[3, :], color = :green, label = "Product (model)")
ode_lineb=lines!(axb, ode_sol.t, ode_sol[1, :], linestyle = :dash, color = :blue, label = "Biomass (model)")
ode_lines=lines!(axs, ode_sol.t, ode_sol[2, :], linestyle = :dash, color = :red, label = "Substrate (model)")
ode_linep=lines!(axp, ode_sol.t, ode_sol[3, :], linestyle = :dash, color = :green, label = "Product (model)")
barplot!(ax_bar, [1, 2], [err1, err2])

axislegend(axb, [fde_lineb, ode_lineb], ["Fractional order", "Integer order"], position = :rb, labelsize=10, rowgap = -5, patchsize = (12, 22))
axislegend(axs, [fde_lines, ode_lines], ["Fractional order", "Integer order"], position = :rt, labelsize=10, rowgap = -5, patchsize = (12, 22))
axislegend(axp, [fde_linep, ode_linep], ["Fractional order", "Integer order"], position = :rb, labelsize=10, rowgap = -5, patchsize = (12, 22))
fig
save("bio_ethanol.svg", fig)


# ========================================================================
# Train-Test Split Validation (Time-Based 80/20 Split)
# ========================================================================
# Assess model generalization by training on early time points
# and testing on later time points. This avoids data leakage in temporal data.
#
# Motivation: In fermentation monitoring, early predictions are valuable for
# process control. Models must generalize to unseen future times.

# Determine train/test split boundary
n_points = length(t)  # Total measurement points: 11
n_train = max(1, floor(Int, 0.8 * n_points))  # 80% for training: 9 points
train_idx = 1:n_train
test_idx = (n_train + 1):n_points

# Partition measurement times
t_train = t[train_idx]  # Times: [0, 2, 4, 6, 8, 10, 12, 16, 28]
t_test = t[test_idx]    # Times: [44, 48]

# Partition experimental data
data_train = data[:, train_idx]  # Training data
data_test = data[:, test_idx]    # Test data

# Helper function to convert time to solution index (for PITrap with dt=1)
time_to_index(tt) = Int(round(tt)) + 1
train_solution_idx = time_to_index.(t_train)  # Convert time points to indices
test_solution_idx = time_to_index.(t_test)

# Loss functions for training data only (not used for final comparison)
function loss_ode_train(b)
    """
    ODE loss evaluated only on training set.
    Used for parameter fitting to avoid test set contamination.
    """
    par_local = copy(b)
    prob_local = ODEProblem(lotka_volterra, X0, tspan, par_local)
    sol_local = solve(prob_local, Tsit5(), saveat=t_train)
    pred_local = reduce(hcat, sol_local.u)
    rmsd(data_train, pred_local; normalize=true)
end

function loss_fode_train(b)
    """
    FODE loss evaluated only on training set.
    Used for parameter fitting to avoid test set contamination.
    """
    par_local = b[4:end]
    order_local = b[1:3]
    prob_local = FODEProblem(lotka_volterra, order_local, X0, tspan, par_local)
    sol_local = solve(prob_local, PITrap(), dt = 1)
    pred_local = [sol_local[i] for i in train_solution_idx]
    pred_mat = reduce(hcat, pred_local)
    rmsd(data_train, pred_mat; normalize=true)
end

# Refit both models using TRAINING data only
Res1_split = optimize(
    loss_ode_train,
    p_lo_1,
    p_up_1,
    p_vec_1,
    Fminbox(BFGS()),
    Optim.Options(outer_iterations = 10, iterations = 10000, show_trace = true, show_every = 1)
)
p1_split = vcat(Optim.minimizer(Res1_split))

Res2_split = optimize(
    loss_fode_train,
    p_lo_2,
    p_up_2,
    p_vec_2,
    Fminbox(BFGS()),
    Optim.Options(outer_iterations = 10, iterations = 10000, show_trace = true, show_every = 1)
)
p2_split = vcat(Optim.minimizer(Res2_split))

# Evaluate both models on FULL time range (for visualization)
ode_split_prob = ODEProblem(lotka_volterra, X0, tspan, p1_split)
ode_split_sol = solve(ode_split_prob, Tsit5(), saveat=t)

# Evaluate on training set specifically
ode_split_train = solve(ode_split_prob, Tsit5(), saveat=t_train)
ode_train_rmse = rmsd(data_train, reduce(hcat, ode_split_train.u); normalize=true)

# Evaluate on test set specifically
ode_split_test = solve(ode_split_prob, Tsit5(), saveat=t_test)
ode_test_rmse = rmsd(data_test, reduce(hcat, ode_split_test.u); normalize=true)

# Repeat for FODE model
fode_split_prob = FODEProblem(lotka_volterra, p2_split[1:3], X0, tspan, p2_split[4:end])
fode_split_sol = solve(fode_split_prob, PITrap(), dt = 1)

fode_split_train = reduce(hcat, [fode_split_sol[i] for i in train_solution_idx])
fode_train_rmse = rmsd(data_train, fode_split_train; normalize=true)

fode_split_test = reduce(hcat, [fode_split_sol[i] for i in test_solution_idx])
fode_test_rmse = rmsd(data_test, fode_split_test; normalize=true)

# ========================================================================
# Visualization of Train-Test Split Results
# ========================================================================
# Create comparison figure showing model predictions on full trajectory
# along with separate train/test error metrics.

fig_split = Figure(size = (700, 490))

# Create four subplots
axb_split = Axis(fig_split[1, 1], xlabel = "Time (h)", ylabel = "Concentration of biomass (g/L)")
axs_split = Axis(fig_split[1, 2], xlabel = "Time (h)", ylabel = "Concentration of substrate (g/L)")
axp_split = Axis(fig_split[2, 1], xlabel = "Time (h)", ylabel = "Concentration of product (g/L)")
ax_split_bar = Axis(
    fig_split[2, 2],
    xlabel = "Dataset",
    xticks = (1:2, ["Train", "Test"]),
    title = "Normalized RMSE"
)
ax_split_bar.xlabelsize = 2

# Add experimental data points
scb_split = CairoMakie.scatter!(axb_split, biomass[!, 1], biomass[!, 2], color = :blue, label = "Experimental Biomass")
scp_split = CairoMakie.scatter!(axp_split, product[!, 1], product[!, 2], color = :green, label = "Experimental Product")
scs_split = CairoMakie.scatter!(axs_split, substrate[!, 1], substrate[!, 2], color = :red, label = "Experimental Substrate")

# Plot model predictions on full time domain
node_lineb_split = lines!(axb_split, ode_split_sol.t, ode_split_sol[1, :], linestyle = :dash, color = :blue, label = "ODE (test)")
ode_lines_split = lines!(axs_split, ode_split_sol.t, ode_split_sol[2, :], linestyle = :dash, color = :red, label = "ODE (test)")
ode_linep_split = lines!(axp_split, ode_split_sol.t, ode_split_sol[3, :], linestyle = :dash, color = :green, label = "ODE (test)")

fode_lineb_split = lines!(axb_split, fode_split_sol.t, fode_split_sol[1, :], color = :blue, label = "FODE (test)")
fode_lines_split = lines!(axs_split, fode_split_sol.t, fode_split_sol[2, :], color = :red, label = "FODE (test)")
fode_linep_split = lines!(axp_split, fode_split_sol.t, fode_split_sol[3, :], color = :green, label = "FODE (test)")

# Plot separate train/test RMSE as grouped bar chart
# Each model gets two bars: one for training error, one for test error
barplot!(
    ax_split_bar,
    [1, 2],
    [ode_train_rmse, ode_test_rmse],
    label = "Integer-order Model",
    dodge = 1,
    n_dodge = 2,
    label_size = 10,
    width = 0.35
)

barplot!(
    ax_split_bar,
    [1, 2],
    [fode_train_rmse, fode_test_rmse],
    label = "Fractional-order Model",
    dodge = 2,
    n_dodge = 2,
    label_size = 10,
    width = 0.35
)

# Add legends distinguishing model types
axislegend(
    axb_split,
    [fode_lineb_split, ode_lineb_split],
    ["Fractional order", "Integer order"],
    position = :rb,
    labelsize = 10,
    rowgap = -5,
    patchsize = (12, 22)
)

axislegend(
    axs_split,
    [fode_lines_split, ode_lines_split],
    ["Fractional order", "Integer order"],
    position = :rt,
    labelsize = 10,
    rowgap = -5,
    patchsize = (12, 22)
)

axislegend(
    axp_split,
    [fode_linep_split, ode_linep_split],
    ["Fractional order", "Integer order"],
    position = :rb,
    labelsize = 10,
    rowgap = -5,
    patchsize = (12, 22)
)

axislegend(ax_split_bar, position = :lt)

fig_split
save("bio_ethanol_train_test_split.svg", fig_split)

# Print train-test validation results
println("\n===== Train-test split validation (time-based 80/20) =====")
println("Train time points: ", collect(t_train))
println("Test time points: ", collect(t_test))
println("ODE train RMSE: ", ode_train_rmse)
println("ODE test RMSE: ", ode_test_rmse)
println("FODE train RMSE: ", fode_train_rmse)
println("FODE test RMSE: ", fode_test_rmse)

# ========================================================================
# AIC/BIC Model Selection Criteria (Information Theory)
# ========================================================================
# Akaike Information Criterion (AIC) and Bayesian Information Criterion (BIC)
# provide quantitative model comparison accounting for both fit quality and
# model complexity. These are derived under maximum likelihood framework.
#
# Assumptions:
#   - Gaussian residuals with constant variance
#   - Log-likelihood derived from residual sum of squares (RSS)
#
# Formulas:
#   AIC  = n * ln(RSS/n) + 2*k
#   BIC  = n * ln(RSS/n) + k * ln(n)
#
# where:
#   k = number of free model parameters (penalizes complexity)
#   n = total number of scalar observations (variables × time points)
#   RSS = residual sum of squares (model-data misfit)
#
# Interpretation:
#   - Lower AIC/BIC indicates better model
#   - ΔAIC > 10 or ΔBIC > 10 indicates strong evidence for better model
#   - AIC penalizes complexity less than BIC (AIC prefers complex models)

# ========================================================================
# Compute predictions on full dataset using optimized parameters
# ========================================================================

sol1_full = solve(ODEProblem(lotka_volterra, X0, tspan, p1), Tsit5(), saveat=t)
pred_ode = reduce(hcat, sol1_full.u)  # Predictions: 3 variables × 11 times

osol_full = solve(FODEProblem(lotka_volterra, p2[1:3], X0, tspan, p2[4:end]), PITrap(), dt = 1)
fode_idxs = time_to_index.(t)  # Convert times to indices
pred_fode = reduce(hcat, [osol_full[i] for i in fode_idxs])  # 3 variables × 11 times

# Compute residuals (model errors)
resid_ode = data .- pred_ode
resid_fode = data .- pred_fode

# Compute residual sum of squares (RSS)
rss_ode = sum(resid_ode .^ 2)
rss_fode = sum(resid_fode .^ 2)

# Count total observations
n_obs = length(data)  # 3 variables × 11 time points = 33 scalar observations

# Count free parameters (model complexity)
k_ode = 4    # ODE model: [kc, km, ks, kp]
k_fode = 7   # FODE model: [α₁, α₂, α₃] + [kc, km, ks, kp]

# ========================================================================
# Compute information criteria
# ========================================================================

aic_ode = n_obs * log(rss_ode / n_obs) + 2 * k_ode
bic_ode = n_obs * log(rss_ode / n_obs) + k_ode * log(n_obs)

aic_fode = n_obs * log(rss_fode / n_obs) + 2 * k_fode
bic_fode = n_obs * log(rss_fode / n_obs) + k_fode * log(n_obs)

# ========================================================================
# Print comprehensive results
# ========================================================================

println("\n===== AIC / BIC comparison (full dataset) =====")
println("Number of observations: ", n_obs)
@printf("ODE  (k=%d) : AIC = %.4f  |  BIC = %.4f\n", k_ode, aic_ode, bic_ode)
@printf("FODE (k=%d) : AIC = %.4f  |  BIC = %.4f\n", k_fode, aic_fode, bic_fode)
println()

# Display model preference and effect sizes
if aic_fode < aic_ode
    println("AIC favors: FODE  (ΔAIC = ", round(aic_ode - aic_fode, digits=4), ")")
else
    println("AIC favors: ODE   (ΔAIC = ", round(aic_fode - aic_ode, digits=4), ")")
end

if bic_fode < bic_ode
    println("BIC favors: FODE  (ΔBIC = ", round(bic_ode - bic_fode, digits=4), ")")
else
    println("BIC favors: ODE   (ΔBIC = ", round(bic_fode - bic_ode, digits=4), ")")
end
