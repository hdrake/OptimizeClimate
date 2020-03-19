# Model domain
present_year = 2020.
dt = 5. # years
t = Array(present_year:dt:2200);
sec_per_year = (365. * 24. * 60. * 60.)

## Physics

# Two-layer EBM (Gregory 2000) parameters from Geoffroy 2013
B = 1.13 * sec_per_year; # Feedback parameter [J yr^-1 m^-2 K^-1]
Cd = 106 * sec_per_year; # Deep ocean heat capacity [J m^-2 K^-1]
κ = 0.73 * sec_per_year; # Heat exchange coefficient [J yr^-1 m^2 K^-1]
δT_init = 1.1 # Berkeley Earth Surface Temperature (Rohde 2013)

# Carbon model
CO₂_init = 460.
r = 0.5 # fraction of emissions remaining after biosphere and ocean uptake (Solomon 2009)

# Economics
GWP0 = 100. # global world product (trillion $ / year)
GWP = GWP0 * exp.((t .- t[1]) / 50.)

β = 0.02/(3.0)^2 # damages (%GWP / celsius^2)
utility_discount_rate = 0.0

# Control technology cost scales, as fraction of GWP (cost scale is for full deployment, α=1.)
mitigate_cost = 0.01*GWP0; # [10^12$ yr^-1] 
remove_cost = 0.05*GWP0; # [10^12$ yr^-1] # Estimate cost from Fuss 2018 (see synthesis Figure 14)
geoeng_cost = 2 * β * (6.0^2); # [%GWP]
adapt_cost = 0.03*GWP0; # [10^12$ yr^-1]

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
    GWP, β, utility_discount_rate,
    mitigate_cost, remove_cost, geoeng_cost, adapt_cost,
    1. /6., 0., 0., nothing, # Initial condition on control deployments at t[1]
    baseline_emissions(t)
)

Economics0(t) = Economics(
    GWP, β, utility_discount_rate,
    mitigate_cost, remove_cost, geoeng_cost, adapt_cost,
    0., 0., 0., nothing, # Initial condition on control deployments at t[1]
    baseline_emissions(t)
)

Economics() = Economics(t)
Economics0() = Economics0(t)

Physics() = Physics(CO₂_init, δT_init, B, Cd, κ, r)

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

ClimateModel() = ClimateModel("default")