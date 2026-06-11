# Dataset subset
using CSV
using Optim
using DataFrames
using OrdinaryDiffEq
using FractionalDiffEq
using StatsBase
using CairoMakie
using Printf

# parameters
kc=0.0041 # Transmission coeﬃcient from infected individuals
km=2.3875*10^(-14) # Relative transmissibility of hospitalized patients
ks=0.0585 # Transmission coeﬃcient due to super-spreaders
kp=0.0156 # Rate at which exposed become infectious
function lotka_volterra(du, u, p, t)
    kc, km, ks, kp = p
    du[1] = kc * u[1] * u[2] - km * u[1]
    du[2] = -ks * u[1] * u[2]
    du[3] = kp * u[1] * u[2]
end
biomass = DataFrame(CSV.File("/Users/quqingyu/SciFracX/paper/applications/bio_ethanol/biomass.csv"))
product = DataFrame(CSV.File("/Users/quqingyu/SciFracX/paper/applications/bio_ethanol/product.csv"))
substrate = DataFrame(CSV.File("/Users/quqingyu/SciFracX/paper/applications/bio_ethanol/substrate.csv"))

X0=[0.5, 90.0, 2.5] # initial values
tspan=(0, 50) # time span [initial time, final time]
data = Matrix(hcat(biomass[!, 2], substrate[!, 2], product[!, 2])')
t = [0,2,4,6,8,10,12,16,28,44,48]
par=[kc, km, ks, kp] # parameters
## optimazation of β for integer order model

function loss_1(b)# loss function
	par=copy(b)
    prob = ODEProblem(lotka_volterra, X0, tspan, par)
	#_, x = FDEsolver(SIR, tspan, X0, ones(8), par, h = .1)
    sol = solve(prob, Tsit5(), saveat=t)
    #appX=vec(sum(x[1:10:end,[3,4,6]], dims=2))
    appX = reduce(hcat, sol.u)
    rmsd(data, appX; normalize=:true) # Normalized root-mean-square error
end

function loss_2(b)# loss function
	par=b[4:end]
    order = b[1:3]
    prob = FODEProblem(lotka_volterra, order, X0, tspan, par)
	#_, x = FDEsolver(SIR, tspan, X0, ones(8), par, h = .1)
    osol = solve(prob, PITrap(), dt = 1)
    #appX=vec(sum(x[1:10:end,[3,4,6]], dims=2))
    app = [osol[1], osol[3], osol[5], osol[7], osol[9], osol[11], osol[13], osol[17], osol[29], osol[45], osol[49]]
    appX = reduce(hcat, app)
    rmsd(data, appX; normalize=:true) # Normalized root-mean-square error
end

function moser_luong_lotka_volterra(du, u, p, t)
    um, ks, kp, ms, yxs, yps, alpha, beta = p
    n=1;m=9
    mu = um*u[2]^n/(ks+u[2]^n)*(1-(u[3]/kp)^m)/(1+exp(-t+1))
    du[1] = mu*u[1]
    du[2] = -(1/yxs + 1/yps) * mu * u[1] - ms * u[1]
    du[3] = alpha*yps*mu*u[1] + beta*u[1]
end
function loss_3(b)# loss function
	par=copy(b)
    prob = ODEProblem(moser_luong_lotka_volterra, X0, tspan, par)
	#_, x = FDEsolver(SIR, tspan, X0, ones(8), par, h = .1)
    sol = solve(prob, Tsit5(), saveat=t)
    #appX=vec(sum(x[1:10:end,[3,4,6]], dims=2))
    appX = reduce(hcat, sol.u)
    rmsd(data, appX; normalize=:true) # Normalized root-mean-square error
end

p_lo_1=[0.003, 0.0, 0.05, 0.015] #lower bound for β
p_up_1=[0.005, 0.0005, 0.06, 0.016] # upper bound for β
p_vec_1=[0.004, 0.0001, 0.055, 0.0155] #  initial guess for β
Res1=optimize(loss_1,p_lo_1,p_up_1,p_vec_1,Fminbox(BFGS()),# Broyden–Fletcher–Goldfarb–Shanno algorithm
# Result=optimize(loss_1,p_lo_1,p_up_1,p_vec_1,SAMIN(rt=.99), # Simulated Annealing algorithm (sometimes it has better perfomance than (L-)BFGS)
			Optim.Options(outer_iterations = 10,
						  iterations=10000,
						  show_trace=true,
						  show_every=1))
p1=vcat(Optim.minimizer(Res1))



p_lo_2=[0.7, 0.8, 0.9, 0.003, 0.0, 0.05, 0.015] #lower bound for β
p_up_2=[0.8, 0.9, 1.0, 0.005, 0.0005, 0.06, 0.016] # upper bound for β
p_vec_2=[0.75, 0.85, 0.95, 0.004, 0.0001, 0.055, 0.0155] #  initial guess for β
Res2=optimize(loss_2,p_lo_2,p_up_2,p_vec_2,Fminbox(BFGS()),# Broyden–Fletcher–Goldfarb–Shanno algorithm
# Result=optimize(loss_1,p_lo_1,p_up_1,p_vec_1,SAMIN(rt=.99), # Simulated Annealing algorithm (sometimes it has better perfomance than (L-)BFGS)
			Optim.Options(outer_iterations = 10,
						  iterations=10000,
						  show_trace=true,
						  show_every=1))
p2=vcat(Optim.minimizer(Res2))





p_lo_3=[0.3, 15, 120, 0.0, 0.0, 0.4, 20, 0.0] #lower bound for β
p_up_3=[0.4, 25, 140, 0.1, 0.1, 0.5, 25, 0.05] # upper bound for β
p_vec_3=[0.35, 19, 125, 0.01, 0.01, 0.45, 24, 0.03] #  initial guess for β
Res3=optimize(loss_3,p_lo_3,p_up_3,p_vec_3,Fminbox(BFGS()),# Broyden–Fletcher–Goldfarb–Shanno algorithm
# Result=optimize(loss_1,p_lo_1,p_up_1,p_vec_1,SAMIN(rt=.99), # Simulated Annealing algorithm (sometimes it has better perfomance than (L-)BFGS)
			Optim.Options(outer_iterations = 10,
						  iterations=10000,
						  show_trace=true,
						  show_every=1))
p3=vcat(Optim.minimizer(Res3))


using Plots
aprob1 = ODEProblem(lotka_volterra, X0, tspan, p1)
ode_sol = solve(aprob1, Tsit5())

aprob2 = FODEProblem(lotka_volterra, p2[1:3], X0, tspan, p2[4:end])
fde_sol = solve(aprob2, PITrap(), dt = 0.01)

sol1_err = solve(aprob1, Tsit5(), saveat=t)
osol = solve(aprob2, PITrap(), dt = 1)

sol2_err = [osol[1], osol[3], osol[5], osol[7], osol[9], osol[11], osol[13], osol[17], osol[29], osol[45], osol[49]]

err1 = rmsd(data, reduce(hcat, sol1_err.u); normalize=true)
err2 = rmsd(data, reduce(hcat, sol2_err); normalize=true)


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


# =========================
# Train-test split experiment
# =========================
# Split by time order: first 80% for training, last 20% for testing.
n_points = length(t)
n_train = max(1, floor(Int, 0.8 * n_points))
train_idx = 1:n_train
test_idx = (n_train + 1):n_points

t_train = t[train_idx]
t_test = t[test_idx]

data_train = data[:, train_idx]
data_test = data[:, test_idx]

# Convert saved times (0,2,4,...) to indices of PITrap(dt=1) solution (1-based).
time_to_index(tt) = Int(round(tt)) + 1
train_solution_idx = time_to_index.(t_train)
test_solution_idx = time_to_index.(t_test)

# ODE model: fit only on training points.
function loss_ode_train(b)
    par_local = copy(b)
    prob_local = ODEProblem(lotka_volterra, X0, tspan, par_local)
    sol_local = solve(prob_local, Tsit5(), saveat=t_train)
    pred_local = reduce(hcat, sol_local.u)
    rmsd(data_train, pred_local; normalize=true)
end

# FODE model: fit only on training points.
function loss_fode_train(b)
    par_local = b[4:end]
    order_local = b[1:3]
    prob_local = FODEProblem(lotka_volterra, order_local, X0, tspan, par_local)
    sol_local = solve(prob_local, PITrap(), dt = 1)
    pred_local = [sol_local[i] for i in train_solution_idx]
    pred_mat = reduce(hcat, pred_local)
    rmsd(data_train, pred_mat; normalize=true)
end

# Refit models using training data only.
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

# Evaluate both models on train and test sets.
ode_split_prob = ODEProblem(lotka_volterra, X0, tspan, p1_split)
ode_split_sol = solve(ode_split_prob, Tsit5(), saveat=t)
ode_split_train = solve(ode_split_prob, Tsit5(), saveat=t_train)
ode_split_test = solve(ode_split_prob, Tsit5(), saveat=t_test)
ode_train_rmse = rmsd(data_train, reduce(hcat, ode_split_train.u); normalize=true)
ode_test_rmse = rmsd(data_test, reduce(hcat, ode_split_test.u); normalize=true)

fode_split_prob = FODEProblem(lotka_volterra, p2_split[1:3], X0, tspan, p2_split[4:end])
fode_split_sol = solve(fode_split_prob, PITrap(), dt = 1)
fode_split_train = reduce(hcat, [fode_split_sol[i] for i in train_solution_idx])
fode_split_test = reduce(hcat, [fode_split_sol[i] for i in test_solution_idx])
fode_train_rmse = rmsd(data_train, fode_split_train; normalize=true)
fode_test_rmse = rmsd(data_test, fode_split_test; normalize=true)

# plot the solution of ODE and FODE with the fitted parameters
fig_split = Figure(size = (700, 490))
axb_split = Axis(fig_split[1, 1], xlabel = "Time (h)", ylabel = "Concentration of biomass (g/L)")
axs_split = Axis(fig_split[1, 2], xlabel = "Time (h)", ylabel = "Concentration of substrate (g/L)")
axp_split = Axis(fig_split[2, 1], xlabel = "Time (h)", ylabel = "Concentration of product (g/L)")
ax_split_bar = Axis(fig_split[2, 2], xlabel = "Dataset", xticks = (1:2, ["Train", "Test"]), title = "Normalized RMSE")
ax_split_bar.xlabelsize=2
scb_split = CairoMakie.scatter!(axb_split, biomass[!, 1], biomass[!, 2], color = :blue, label = "Experimental Biomass")
scp_split = CairoMakie.scatter!(axp_split, product[!, 1], product[!, 2], color = :green, label = "Experimental Product")
scs_split = CairoMakie.scatter!(axs_split, substrate[!, 1], substrate[!, 2], color = :red, label = "Experimental Substrate")
ode_lineb_split = lines!(axb_split, ode_split_sol.t, ode_split_sol[1, :], linestyle = :dash, color = :blue, label = "ODE (test)")
ode_lines_split = lines!(axs_split, ode_split_sol.t, ode_split_sol[2, :], linestyle = :dash, color = :red, label = "ODE (test)")
ode_linep_split = lines!(axp_split, ode_split_sol.t, ode_split_sol[3, :], linestyle = :dash, color = :green, label = "ODE (test)")
fode_lineb_split = lines!(axb_split, fode_split_sol.t, fode_split_sol[1, :], color = :blue, label = "FODE (test)")
fode_lines_split = lines!(axs_split, fode_split_sol.t, fode_split_sol[2, :], color = :red, label = "FODE (test)")
fode_linep_split = lines!(axp_split, fode_split_sol.t, fode_split_sol[3, :], color = :green, label = "FODE (test)")
barplot!(ax_split_bar, [1, 2], [ode_train_rmse, ode_test_rmse], label = "Integer-order Model", dodge = 1, n_dodge = 2, label_size=10, width = 0.35)
barplot!(ax_split_bar, [1, 2], [fode_train_rmse, fode_test_rmse], label = "Fractional-order Model", dodge = 2, n_dodge = 2, label_size=10, width = 0.35)
axislegend(axb_split, [fode_lineb_split, ode_lineb_split], ["Fractional order", "Integer order"], position = :rb, labelsize=10, rowgap = -5, patchsize = (12, 22))
axislegend(axs_split, [fode_lines_split, ode_lines_split], ["Fractional order", "Integer order"], position = :rt, labelsize=10, rowgap = -5, patchsize = (12, 22))
axislegend(axp_split, [fode_linep_split, ode_linep_split], ["Fractional order", "Integer order"], position = :rb, labelsize=10, rowgap = -5, patchsize = (12, 22))
axislegend(
    ax_split_bar,
    position = :lt,
)
fig_split
save("bio_ethanol_train_test_split.svg", fig_split)

println("\n===== Train-test split validation (time-based 80/20) =====")
println("Train time points: ", collect(t_train))
println("Test time points: ", collect(t_test))
println("ODE train RMSE: ", ode_train_rmse)
println("ODE test RMSE: ", ode_test_rmse)
println("FODE train RMSE: ", fode_train_rmse)
println("FODE test RMSE: ", fode_test_rmse)












# =========================
# AIC / BIC model comparison (full dataset)
# =========================
# Both criteria assume Gaussian residuals; log-likelihood is derived from RSS.
#
#   AIC  = n * ln(RSS/n) + 2k
#   BIC  = n * ln(RSS/n) + k * ln(n)
#
# k   = number of free parameters
# n   = total number of scalar observations (n_vars × n_timepoints)
# RSS = sum of squared residuals over all variables and time points

# --- Predictions on the full dataset ---
sol1_full = solve(ODEProblem(lotka_volterra, X0, tspan, p1), Tsit5(), saveat=t)
pred_ode = reduce(hcat, sol1_full.u)              # 3 × 11

osol_full = solve(FODEProblem(lotka_volterra, p2[1:3], X0, tspan, p2[4:end]), PITrap(), dt = 1)
fode_idxs = time_to_index.(t)
pred_fode = reduce(hcat, [osol_full[i] for i in fode_idxs])  # 3 × 11

# --- Residuals and RSS ---
resid_ode  = data .- pred_ode
resid_fode = data .- pred_fode

rss_ode  = sum(resid_ode .^ 2)
rss_fode = sum(resid_fode .^ 2)

# --- Observation count ---
n_obs = length(data)   # 3 variables × 11 time points = 33

# --- Parameter counts ---
# ODE:  4 model parameters (kc, km, ks, kp)
# FODE: 4 model parameters + 3 fractional orders = 7
k_ode  = 4
k_fode = 7

# --- AIC and BIC ---
aic_ode  = n_obs * log(rss_ode  / n_obs) + 2 * k_ode
bic_ode  = n_obs * log(rss_ode  / n_obs) + k_ode  * log(n_obs)

aic_fode = n_obs * log(rss_fode / n_obs) + 2 * k_fode
bic_fode = n_obs * log(rss_fode / n_obs) + k_fode * log(n_obs)

println("\n===== AIC / BIC comparison (full dataset) =====")
println("Number of observations: ", n_obs)
@printf("ODE  (k=%d) : AIC = %.4f  |  BIC = %.4f\n", k_ode,  aic_ode,  bic_ode)
@printf("FODE (k=%d) : AIC = %.4f  |  BIC = %.4f\n", k_fode, aic_fode, bic_fode)
println()
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











