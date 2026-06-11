using QuasiMonteCarlo
using FractionalDiffEq
using Statistics
using CairoMakie

# -----------------------------
# Basic setup
# -----------------------------

param_names = ["m", "c", "α", "k", "F₀", "ω"]

lb = [0.8, 0.5, 1.1, 8.0, 0.8, 1.0]
ub = [1.2, 8.0, 1.9, 12.0, 1.2, 1.5]

saveat = collect(0.0:0.02:20.0)
N_t = length(saveat)

# -----------------------------
# Resampling
# -----------------------------

function linear_resample(t, u, saveat)
    y = similar(saveat, Float64)

    k = 1
    for (j, tj) in pairs(saveat)
        while k < length(t) - 1 && t[k + 1] < tj
            k += 1
        end

        if tj <= t[1]
            y[j] = u[1]
        elseif tj >= t[end]
            y[j] = u[end]
        else
            θ = (tj - t[k]) / (t[k + 1] - t[k])
            y[j] = (1 - θ) * u[k] + θ * u[k + 1]
        end
    end

    return y
end

# -----------------------------
# Fractional oscillator solver
# -----------------------------

function simulate_fractional_oscillator(p, saveat)
    m, c, α, k, F0, ω = p

    tspan = (0.0, 20.0)
    x0 = 0.0
    v0 = 0.0

    prob = MultiTermsFODEProblem(
        [m, c, k],
        [2.0, α, 0.0],
        (u, p, t) -> F0 * sin(ω * t),
        [x0, v0],
        tspan
    )

    sol = solve(prob, MTPITrap(), dt = 0.01)

    t_raw = sol.t
    u_raw = Array(sol.u)

    x_raw = if eltype(u_raw) <: Number
        u_raw
    else
        [ui[1] for ui in sol.u]
    end

    return linear_resample(t_raw, x_raw, saveat)
end

# ============================================================
# Part 1: vary α only
# ============================================================

p_base = [
    1.0,   # m
    3.0,   # c
    1.5,   # α
    10.0,  # k
    1.0,   # F₀
    1.2    # ω
]

α_values = [1.1, 1.3, 1.5, 1.7, 1.9]

fig = Figure(size = (1500, 520), fontsize = 18)

ax1 = Axis(
    fig[1, 1],
    xlabel = "Time",
    ylabel = "Displacement",
    title = "Effect of fractional order on system response"
)
ylims!(ax1, -0.5, 0.5)
for α in α_values
    p = copy(p_base)
    p[3] = α
    x = simulate_fractional_oscillator(p, saveat)
    lines!(ax1, saveat, x, linewidth = 2, label = L"\alpha = %$α")
end

axislegend(ax1, position = :rt)

# ============================================================
# Part 2: uncertainty propagation with all parameters
# ============================================================

N = 1000
samples_unit = QuasiMonteCarlo.sample(N, 6, SobolSample())

samples = [
    lb[j] + samples_unit[j, i] * (ub[j] - lb[j])
    for j in 1:6, i in 1:N
]

trajectories = zeros(N_t, N)

for i in 1:N
    p = samples[:, i]
    trajectories[:, i] = simulate_fractional_oscillator(p, saveat)
end

mean_traj = vec(mean(trajectories, dims = 2))
median_traj = [median(trajectories[j, :]) for j in 1:N_t]

lower_05 = [quantile(trajectories[j, :], 0.05) for j in 1:N_t]
upper_95 = [quantile(trajectories[j, :], 0.95) for j in 1:N_t]

lower_25 = [quantile(trajectories[j, :], 0.25) for j in 1:N_t]
upper_75 = [quantile(trajectories[j, :], 0.75) for j in 1:N_t]
fig2 = Figure(size = (780, 460), fontsize = 18)

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

# Outer uncertainty interval

band!(
    ax2,
    saveat,
    lower_05,
    upper_95,
    color = (:steelblue, 0.28),
    label = "90% uncertainty interval"
)

# Inner uncertainty interval

band!(
    ax2,
    saveat,
    lower_25,
    upper_75,
    color = (:orange, 0.35),
    label = "50% interquartile interval"
)

# Mean and median

lines!(
    ax2,
    saveat,
    mean_traj,
    color = :black,
    linewidth = 3,
    label = "Mean response"
)

lines!(
    ax2,
    saveat,
    median_traj,
    color = :firebrick,
    linestyle = :dash,
    linewidth = 2.5,
    label = "Median response"
)

# Zero reference line

hlines!(
    ax2,
    [0.0],
    color = (:gray, 0.5),
    linestyle = :dot,
    linewidth = 1.5
)

# set y limits to better show the uncertainty bands
ylims!(ax2, -0.6, 0.6)
axislegend(
    ax2,
    position = :rt,
    framevisible = true,
    backgroundcolor = (:white, 0.85),
    labelsize = 15
)

save("fractional_uncertainty.pdf", fig)

save("fractional_uncertainty.png", fig)

fig