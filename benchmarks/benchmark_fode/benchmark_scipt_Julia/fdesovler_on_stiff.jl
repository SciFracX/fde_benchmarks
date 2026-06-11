using FdeSolver

beta = [0.5, 0.5, 0.5]
β = 0.5
#σ₁ = β; σ₂ = 2β; σ₃ = 1+β; σ₄ = 5β; σ₅ = 2; σ₆ = 2+β
σ = [β, 2β, 1+β, 5β, 2, 2+β]
a1 = 0.5; a2 = 0.8; a3 = 1; a4 = 1; a5 = 1; a6 = 1

A = [-10000 0 1;
     -0.05 -0.08 -0.2;
     1 0 -1]
B = [-0.6 0 0.2;
     -0.1 -0.2 0;
     0 -0.5 -0.8]
u0 = [1, 1, 1]
Γ(k) = gamma(σ[k]+1)/gamma(σ[k]-β+1)
g(t) = [a1*Γ(1)*t^(σ[1]-β) + a2*Γ(2)*t^(σ[2]-β); a3*Γ(3)*t^(σ[3]-β) + a4*Γ(4)*t^(σ[4]-β); a5*Γ(5)*t^(σ[5]-β) + a6*Γ(6)*t^(σ[6]-β)] - (A + B) * [a1*t^σ[1] + a2*t^σ[2] + u0[1]; a3*t^σ[3] + a4*t^σ[4] + u0[2]; a5*t^σ[5] + a6*t^σ[6] + u0[3]]

# Inputs
I0 = 0.001;             # intial value of infected
tSpan = [0, 100];       # [intial time, final time]
y0 = [1 , 1, 1];   # initial values [S0,I0,R0]
α = [0.5, 0.5, 0.5];          # order of derivatives
h = 0.1;                # step size of computation (default = 0.01)
par = [0.4, 0.04];      # parameters [β, recovery rate]

## ODE model
function F(t, y, par)

    # parameters
    β = par[1]    # infection rate
    γ = par[2]    # recovery rate

    S = y[1]   # Susceptible
    I = y[2]   # Infectious
    R = y[3]   # Recovered

    # System equation
    dSdt = - β .* S .* I
    dIdt = β .* S .* I .- γ .* I
    dRdt = γ .* I

    return [dSdt, dIdt, dRdt]

end

## Jacobian of ODE system
function JacobF(t, y, par)

    # parameters
    β = par[1]     # infection rate
    γ = par[2]     # recovery rate

    S = y[1]    # Susceptible
    I = y[2]    # Infectious
    R = y[3]    # Recovered

    # System equation
    J11 = - β * I
    J12 = - β * S
    J13 =  0
    J21 =  β * I
    J22 =  β * S - γ
    J23 =  0
    J31 =  0
    J32 =  γ
    J33 =  0

    J = [J11 J12 J13
         J21 J22 J23
         J31 J32 J33]

    return J

end

## Solution
t, Yapp = FDEsolver(F, tSpan, y0, α, par, JF = JacobF, h = h);

# Plot
plot(t, Yapp, linewidth = 5, title = "Numerical solution of SIR model",
     xaxis = "Time (t)", yaxis = "SIR populations", label = ["Susceptible" "Infectious" "Recovered"]);