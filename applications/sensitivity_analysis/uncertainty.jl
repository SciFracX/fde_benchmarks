"""
    Uncertainty and Sensitivity Analysis for Fractional-Order Oscillator Systems
    
    This script performs a comprehensive uncertainty quantification and sensitivity analysis
    on a fractional-order damped harmonic oscillator system subject to periodic forcing.
    
    Key Components:
    1. Parameter Space Definition: Defines the ranges for six system parameters including
       the fractional differentiation order α, which is the key parameter of interest.
    
    2. Uncertainty Propagation: Uses quasi-Monte Carlo sampling methods to propagate
       parametric uncertainties through the nonlinear fractional differential equations.
    
    3. Multi-Order Analysis: Illustrates how variations in the fractional derivative order
       affect the fundamental dynamical behavior of the system.
    
    4. Statistical Analysis: Computes quantile-based uncertainty bands, mean/median
       trajectories, and provides visual representation of uncertainty intervals.
"""

using QuasiMonteCarlo
using FractionalDiffEq
using Statistics
using CairoMakie

# ========================================================================
# Basic Setup and Parameter Definitions
# ========================================================================
# This section defines the parameter space for the fractional-order
# oscillator system. The system model is:
#     m*d²x/dt² + c*dᵅx/dtᵅ + k*x = F₀*sin(ωt)
# where:
#   - m: mass (inertial coefficient)
#   - c: damping coefficient
#   - α: fractional derivation order (key parameter, typically α ∈ [1, 2])
#   - k: stiffness coefficient
#   - F₀: forcing amplitude
#   - ω: excitation frequency

param_names = [
    "m",   # Mass (kg) - controls system inertia
    "c",   # Damping coefficient - controls energy dissipation
    "α",   # Fractional order - defines memory/hereditary behavior
    "k",   # Stiffness coefficient - controls restoring force
    "F₀",  # Forcing amplitude - magnitude of external excitation
    "ω"    # Forcing frequency - frequency of external excitation
]

# Lower bounds for each parameter (normalized around typical values)
lb = [
    0.8,   # Lower bound for mass
    0.5,   # Lower bound for damping coefficient
    1.1,   # Lower bound for fractional order (slightly above integer derivative)
    8.0,   # Lower bound for stiffness
    0.8,   # Lower bound for forcing amplitude
    1.0    # Lower bound for forcing frequency
]

# Upper bounds for each parameter
ub = [
    1.2,   # Upper bound for mass
    8.0,   # Upper bound for damping coefficient (allows wide variation)
    1.9,   # Upper bound for fractional order (approaching integer order 2)
    12.0,  # Upper bound for stiffness
    1.2,   # Upper bound for forcing amplitude
    1.5    # Upper bound for forcing frequency
]

# Temporal discretization for solution output
# Output times: from t=0 to t=20 with uniform spacing of 0.02 time units
# This resolution (500 time points) captures the transient and steady-state behavior
saveat = collect(0.0:0.02:20.0)
N_t = length(saveat)  # Total number of time points for analysis

# ========================================================================
# Solution Interpolation and Resampling
# ========================================================================
# The fractional differential equation solver produces solutions at
# internally determined time points. This function performs linear
# interpolation to evaluate the solution at specified output times (saveat).
# This is necessary because the solver's output grid may differ from the
# desired analysis grid.

function linear_resample(t, u, saveat)
    """
        linear_resample(t, u, saveat)
    
    Performs linear interpolation on solution data.
    
    Arguments:
        t::Vector: Original time points from the solver
        u::Vector: Solution values at the original time points
        saveat::Vector: Desired output time points for resampling
    
    Returns:
        y::Vector: Interpolated solution values at requested time points
    
    Algorithm:
        Uses a forward-searching algorithm with linear interpolation between
        adjacent data points. For each requested time point tj in saveat,
        we find the bracketing original time points t[k] and t[k+1],
        then compute interpolation weight θ and blend the surrounding values.
    """
    y = similar(saveat, Float64)

    k = 1  # Index tracking current position in original time grid
    for (j, tj) in pairs(saveat)
        # Advance k to find the correct bracket for current time point tj
        while k < length(t) - 1 && t[k + 1] < tj
            k += 1
        end

        # Handle boundary conditions and compute interpolation
        if tj <= t[1]
            # Before first point: use first value
            y[j] = u[1]
        elseif tj >= t[end]
            # After last point: use last value
            y[j] = u[end]
        else
            # Interior point: perform linear interpolation
            # θ ∈ [0,1] is the normalized position between t[k] and t[k+1]
            θ = (tj - t[k]) / (t[k + 1] - t[k])
            # Linear blend: y = (1-θ)*u[k] + θ*u[k+1]
            y[j] = (1 - θ) * u[k] + θ * u[k + 1]
        end
    end

    return y
end

# ========================================================================
# Fractional-Order Oscillator System Solver
# ========================================================================
# This function solves the multi-term fractional damped oscillator equation:
#     m*d²x/dt² + c*dᵅx/dtᵅ + k*x = F₀*sin(ωt)
# The system is converted to state-space form for numerical integration.

function simulate_fractional_oscillator(p, saveat)
    """
        simulate_fractional_oscillator(p, saveat)
    
    Solves a fractional-order forced oscillator system.
    
    Arguments:
        p::Vector: Parameter vector [m, c, α, k, F₀, ω]
                   - m: mass (inertial parameter)
                   - c: damping coefficient
                   - α: fractional derivative order
                   - k: stiffness
                   - F₀: forcing amplitude
                   - ω: forcing frequency
        saveat::Vector: Time points at which to return solution values
    
    Returns:
        x::Vector: Displacement response at requested time points
    
    Implementation Details:
        - Uses multi-term fractional ODE formulation
        - Employs MTPITrap solver (Multi-Term Predictor-Integrator Trap-rule method)
        - dt=0.01 provides adequate temporal resolution for accuracy
        - Solution is resampled to match requested time grid
    """
    # Unpack parameter vector
    m, c, α, k, F0, ω = p

    # Temporal integration domain: from start to steady-state
    tspan = (0.0, 20.0)
    
    # Initial conditions: system starts at rest
    x0 = 0.0   # Initial displacement
    v0 = 0.0   # Initial velocity

    # Define the multi-term fractional ODE problem
    # The equation is rewritten as: m*d²x/dt² + c*dᵅx/dtᵅ + k*x = F₀*sin(ωt)
    # In multi-term format: [m, c, k] are coefficients (order of precedence)
    #                       [2.0, α, 0.0] are corresponding derivative orders
    #                       These form: m*Dᵗ² + c*Dᵗᵅ + k*Dᵗ⁰ where Dᵗ denotes fractional differentiation
    prob = MultiTermsFODEProblem(
        [m, c, k],                    # Coefficient vector
        [2.0, α, 0.0],                # Derivative order vector (α ∈ [1,2])
        (u, p, t) -> F0 * sin(ω * t), # Forcing function: periodic sine wave
        [x0, v0],                      # Initial conditions [displacement, velocity]
        tspan                          # Time span for integration
    )

    # Solve using Multi-Term Predictor-Integrator Trapezoidal rule
    # This is a high-order numerical method suitable for fractional equations
    sol = solve(prob, MTPITrap(), dt = 0.01)

    # Extract raw solution data from solver
    t_raw = sol.t      # Original time grid from solver
    u_raw = Array(sol.u)  # Solution array

    # Extract displacement component (first state variable)
    x_raw = if eltype(u_raw) <: Number
        # If solution is scalar, use directly
        u_raw
    else
        # If solution is multi-dimensional (state vector), extract first component
        [ui[1] for ui in sol.u]
    end

    # Resample solution to target output grid (saveat)
    return linear_resample(t_raw, x_raw, saveat)
end

# ========================================================================
# Parametric Sensitivity Analysis - Varying Fractional Order α
# ========================================================================
# This section investigates how the fractional derivative order α affects
# the system's frequency response and damping characteristics. By varying α
# while keeping other parameters fixed, we can isolate the influence of
# memory effects and non-local temporal effects on system dynamics.
#
# The fractional order α controls the nature of the damping:
#   - α = 1: Standard Newtonian viscous damping (~velocity)
#   - α ∈ (1,2): Anomalous/memory-based damping (fractional behavior)
#   - α = 2: Pure inertial response (no damping)

p_base = [
    1.0,   # m:  mass = 1.0 kg
    3.0,   # c:  damping coefficient = 3.0
    1.5,   # α:  baseline fractional order = 1.5 (mid-range)
    10.0,  # k:  stiffness = 10.0 N/m
    1.0,   # F₀: forcing amplitude = 1.0 N
    1.2    # ω:  forcing frequency = 1.2 rad/s
]

# Test values for fractional order spanning the entire range of interest
# These values demonstrate the transition from anomalous damping to classical behavior
α_values = [1.1, 1.3, 1.5, 1.7, 1.9]  # Varies from strongly fractional to nearly classical

# Create visualization figure
fig = Figure(size = (1500, 520), fontsize = 18)

# First subplot: α sensitivity analysis
ax1 = Axis(
    fig[1, 1],
    xlabel = "Time",
    ylabel = "Displacement",
    title = "Effect of fractional order α on system response",
    xgridvisible = true,
    ygridvisible = true,
)
ylims!(ax1, -0.5, 0.5)

# Simulate and plot trajectories for each α value
for α in α_values
    p = copy(p_base)        # Copy base parameters
    p[3] = α                # Override α with current test value
    x = simulate_fractional_oscillator(p, saveat)  # Solve the system
    lines!(ax1, saveat, x, linewidth = 2, label = L"\alpha = %$α")
end

# Configure legend to show current α values
axislegend(ax1, position = :rt, framevisible = true)

# ========================================================================
# Uncertainty Quantification via Quasi-Monte Carlo Sampling
# ========================================================================
# This section performs a comprehensive uncertainty analysis by:
# 1. Sampling the 6-dimensional parameter space using quasi-random sequences
# 2. Propagating parametric uncertainty through the nonlinear system
# 3. Computing statistical moments (mean, quantiles) of the response
# 4. Visualizing the impact of combined parameter variations
#
# Quasi-Monte Carlo (Sobol sequences) provides:
#   - Better space-filling properties than random sampling
#   - Faster convergence for small sample sizes
#   - Low-discrepancy point distributions

N = 1000  # Number of parameter samples for uncertainty analysis

# Generate quasi-random samples in unit hypercube [0,1]^6 using Sobol sequence
# Sobol sequence is a low-discrepancy sequence optimal for multidimensional sampling
samples_unit = QuasiMonteCarlo.sample(N, 6, SobolSample())

# Scale unit hypercube samples to actual parameter ranges [lb, ub]
# For each parameter j and sample i: param[j,i] = lb[j] + samples_unit[j,i] * (ub[j] - lb[j])
samples = [
    lb[j] + samples_unit[j, i] * (ub[j] - lb[j])
    for j in 1:6, i in 1:N
]

# Initialize storage for solution trajectories
# trajectories[t, i] = displacement at time t for i-th parameter sample
trajectories = zeros(N_t, N)

# Forward propagate parametric uncertainty through the fractional oscillator
for i in 1:N
    p = samples[:, i]  # Extract i-th parameter sample
    trajectories[:, i] = simulate_fractional_oscillator(p, saveat)  # Solve system
end

# Compute statistical summaries across all sample trajectories
# At each time point, we have N values from N parameter realizations

# Mean displacement trajectory: E[x(t)] over all parameter samples
mean_traj = vec(mean(trajectories, dims = 2))

# Median displacement trajectory: 50th percentile of response ensemble
median_traj = [median(trajectories[j, :]) for j in 1:N_t]

# 90% Confidence interval: captures central 90% of ensemble predictions
# lower_05: 5th percentile (lower bound, 5% of samples below)
lower_05 = [quantile(trajectories[j, :], 0.05) for j in 1:N_t]
# upper_95: 95th percentile (upper bound, 5% of samples above)
upper_95 = [quantile(trajectories[j, :], 0.95) for j in 1:N_t]

# Interquartile range (IQR): captures central 50% of ensemble predictions
# lower_25: 25th percentile (lower quartile)
lower_25 = [quantile(trajectories[j, :], 0.25) for j in 1:N_t]
# upper_75: 75th percentile (upper quartile)
upper_75 = [quantile(trajectories[j, :], 0.75) for j in 1:N_t]
# Create second figure for uncertainty visualization
fig2 = Figure(size = (780, 460), fontsize = 18)

# Second subplot: Uncertainty quantification visualization
ax2 = Axis(
    fig[1, 2],
    xlabel = "Time",
    ylabel = "Displacement",
    title = "Propagation of parameter uncertainty in the fractional oscillator",
    titlesize = 20,
    xlabelsize = 18,
    ylabelsize = 18,
    xticklabelsize = 15,
    yticklabelsize = 15,
    xgridvisible = true,
    ygridvisible = true,
)

# Plot outer uncertainty band (90% confidence interval)
# This shaded region represents where 90% of all possible responses lie
# The 5% outside this band represents extreme parameter combinations
band!(
    ax2,
    saveat,
    lower_05,      # 5th percentile (lower boundary)
    upper_95,      # 95th percentile (upper boundary)
    color = (:steelblue, 0.28),
    label = "90% uncertainty interval"  # Captures main uncertainty region
)

# Plot inner uncertainty band (50% interquartile range)
# This shaded region represents the "typical" behavior
# It shows where the central half of parameter combinations lead
band!(
    ax2,
    saveat,
    lower_25,      # 25th percentile (lower quartile)
    upper_75,      # 75th percentile (upper quartile)
    color = (:orange, 0.35),
    label = "50% interquartile interval"  # Central tendency region
)

# Plot ensemble mean trajectory
# This is the expected value of the response across all parameter samples
lines!(
    ax2,
    saveat,
    mean_traj,
    color = :black,
    linewidth = 3,
    label = "Mean response",  # E[x(t)]
)

# Plot ensemble median trajectory  
# For skewed distributions, median may differ from mean
# Median is robust to outliers and extreme parameter values
lines!(
    ax2,
    saveat,
    median_traj,
    color = :firebrick,
    linestyle = :dash,
    linewidth = 2.5,
    label = "Median response",  # 50th percentile
)

# Add horizontal reference line at zero displacement
# This helps visualize oscillation behavior around equilibrium
hlines!(
    ax2,
    [0.0],
    color = (:gray, 0.5),
    linestyle = :dot,
    linewidth = 1.5
)

# Set y-axis limits to properly display uncertainty bands
ylims!(ax2, -0.6, 0.6)

# Configure and position the legend
axislegend(
    ax2,
    position = :rt,           # Position: right-top
    framevisible = true,      # Show legend border
    backgroundcolor = (:white, 0.85),  # Semi-transparent white background
    labelsize = 15
)

# Export combined visualization to file formats for publication/presentation
# The composite figure contains both sensitivity analysis and uncertainty quantification

# Save as PDF
save("fractional_uncertainty.pdf", fig)

# Save as PNG
save("fractional_uncertainty.png", fig)

# Display the figure in the REPL for immediate viewing
fig