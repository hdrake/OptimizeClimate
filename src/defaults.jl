# Model domain
present_year = 2020.
dt = 1. # 1-year timestep, can make longer to speed up the model
t = Array(present_year:dt:2200);

# Physics
CO₂_init = 415.
δT_init = 1.1
ECS = 3.0; # "Best-guess equilibrium climate sensitivity"
ocean_fraction = 0.71
H = 4000. * ocean_fraction; # effective depth of deep ocean [m]
ρ = 1000.; # density of liquid water [kg m^-3]
Cp = 4180.0; # specific heat capacity of liquid  water [J kg^-1 K^-1]
Cd = Cp * ρ * H # upper ocean heat capacity
τd = 200. # deep ocean relaxation time scale [years] (similar to Gregory 2000, Held 2009)
γ = Cd / τd

# Economics
GWP = 100. # global world product (trillion $ / year)
β = 0.02*GWP/(3.0)^2 # damages (trillion USD / year / celsius^2)
utility_discount_rate = 0.0 # ρ (relative low value from Stern review)

# Control technology cost scales, as fraction of GWP (cost scale is for full deployment, α=1.)
mitigate_cost = 0.01*GWP;
remove_cost = 0.02*GWP;
geoeng_cost = 0.05*GWP;
adapt_cost = 0.03*GWP;

"""
    Economics()

Create data structure for economic input parameters for `ClimateModel` struct with default values.

Default parameters are:
- `β`= 0.222 × 10^12 USD / (°C)^2
- `utility_discount_rate` = 0.0 (compare with Stern review median value of 1.4% and ~3% Nordhaus values)
- `mitigate_cost` = 1. × 10^12 USD
- `remove_cost` = 2. × 10^12 USD
- `geoeng_cost` = 5. × 10^12 USD
- `adapt_cost` = 3. × 10^12 USD
- `[control]_init` = 0.
- `baseline_emissions` = baseline_emissions(t::Array{Float64,1}, 10., 2080., 40.)

The default baseline emissions scenario corresponds to flat emissions of 10 ppm / year
from 2020 to 2080 and linearly decreasing from 10 ppm / year in 2080 to 0 ppm / year in 2120.

See also: [`ClimateModel`](@ref), [`baseline_emissions`](@ref)
"""
Economics(t) = Economics(
    β, utility_discount_rate,
    mitigate_cost, remove_cost, geoeng_cost, adapt_cost,
    0., 0., 0., 0., # Assumed initial condition of zero control deployments in 2020
    baseline_emissions(t)
)
Economics() = Economics(t)


Physics(ECS) = Physics(ECS::Float64, CO₂_init, δT_init, Cd, γ)
Physics() = Physics(ECS, CO₂_init, δT_init, Cd, γ)

ClimateModel(name::String) = ClimateModel(
    name,
    t,
    dt,
    present_year,
    Economics(),
    Physics(),
    init_zero_controls(t)
)

ClimateModel(;t::Array{Float64,1}, dt::Float64) = ClimateModel(
    "default",
    t,
    dt,
    present_year,
    Economics(t),
    Physics(),
    init_zero_controls(t)
)

ClimateModel(;ECS::Float64) = ClimateModel(
    "default",
    t,
    dt,
    present_year,
    Economics(t),
    Physics(ECS),
    init_zero_controls(t)
)


ClimateModel() = ClimateModel("default")