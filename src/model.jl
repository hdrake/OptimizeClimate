
"""
    Controls(reduce, remove, geoeng, adapt)

Create data structure for climate control arrays Arrays along time axis,
with bounded values ∈ [0,1].

See also: [`ClimateModel`](@ref)
"""
struct Controls
    reduce::Array{Float64,1}
    remove::Array{Float64,1}
    geoeng::Array{Float64,1}
    adapt::Array{Float64,1}
end

"""
    Economics(
        β, utility_discount_rate,
        reduce_cost, remove_cost, geoeng_cost, adapt_cost,
        reduce_init, remove_init, geoeng_init, adapt_init,
        baseline_emissions
    )

Create data structure for economic input parameters for `ClimateModel` struct,
including a baseline emissions scenario.

### Arguments
- `β::Float64`: climate damage parameter [10^12 USD / (°C)^2].
- `utility_discount_rate::Float64`: typically denoted ρ in economic references.
- `[control]_cost::Float64`: scaling cost of full control deployment [10^12 USD].
- `[control]_init::Float64`: fixed initial condition for control deployment [10^12 USD].
- `baseline_emissions::Array{Float64,1}`: prescribed baseline CO₂ equivalent emissions [ppm / yr].

See also: [`ClimateModel`](@ref), [`baseline_emissions`](@ref).

"""
struct Economics
    β::Float64
    utility_discount_rate::Float64
    
    reduce_cost::Float64
    remove_cost::Float64
    geoeng_cost::Float64
    adapt_cost::Float64
    
    reduce_init::Float64
    remove_init::Float64
    geoeng_init::Float64
    adapt_init::Float64
    
    baseline_emissions::Array{Float64,1}
end

"""
    Economics()

Create data structure for economic input parameters for `ClimateModel` struct with default values.

Default parameters are:
- `β`=1. × 10^12 USD / (°C)^2
- `utility_discount_rate` = 0.014 (roughly Stern review median value of 1.4%)
- `reduce_cost` = 5. × 10^12 USD
- `remove_cost` = 5. × 10^12 USD
- `geoeng_cost` = 25. × 10^12 USD
- `adapt_cost` = 15. × 10^12 USD
- `[control]_init` = 0. USD
- `baseline_emissions` = baseline_emissions(t::Array{Float64,1}, 5., 2060., 40.)

The default baseline emissions scenario corresponds to flat emissions of 5 ppm / year
from 2020 to 2060 and linearly decreasing from 5 ppm / year in 2060 to 0 ppm / year in 2100.

See also: [`ClimateModel`](@ref), [`baseline_emissions`](@ref)
"""
Economics() = Economics(
    1., 0.014,
    0.05*100., 0.05*100., 0.25*100., 0.15*100.,
    0., 0., 0., 0.,
    baseline_emissions(Array(2020:1.:2100))
)

"Return a non-dimensional Array of size(t) which represents a linear increase from 0. to 1."
nondim_linear(t::Array) = (t .- t[1])/(t[end] - t[1]);

"""
    init_linear_controls(t)

Create physical (but arbitrary) initial guess of climate controls in which all controls increase
linearly with time.

See also: [`Controls`](@ref)
"""
function init_linear_controls(t::Array{Float64,1})
    c = Controls(
        nondim_linear(t)/3.,
        nondim_linear(t)/2.,
        nondim_linear(t)/8.,
        nondim_linear(t)/10.
    )
    return c
end

"""
    ClimateModel(name, ECS, domain, controls, economics, present_year, CO₂_init, δT_pre, ϵ)

Create instance of an extremely idealized integrated-assessment
climate model in which the climate response (`CO₂` and `temperature`) is a function
of physical input parameters (`CO₂_init`, `δT_pre`, `ECS`, `ϵ`), economic input parameters
(`economics`), and climate control policies (`controls`) over some time frame (`domain`).

See also: [`Controls`](@ref), [`Economics`](@ref), [`CO₂`](@ref), [`δT`](@ref),
[`optimize!`](@ref)
"""
struct ClimateModel
    name::String
    ECS::Float64
    domain::Array{Float64,1}
    controls::Controls
    economics::Economics
    present_year::Float64
    CO₂_init::Float64
    δT_pre::Float64
    
    ϵ::Float64
    
    function ClimateModel(name, ECS, domain, controls, economics, present_year, CO₂_init, δT_pre)
        
        ϵ = ECS/log(2.); # Transient Warming Parameter
        
        return new(
            name, ECS, domain, controls, economics, present_year, CO₂_init, δT_pre,
            ϵ)
    end
end

ClimateModel(name, ECS, domain, controls, economics, present_year) = ClimateModel(
    name, ECS, domain, controls, economics, present_year, 415., 1.1
)