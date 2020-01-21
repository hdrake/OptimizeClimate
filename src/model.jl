
struct Physics
    ECS::Float64
    CO₂_init::Float64
    δT_init::Float64
    Cd::Float64
    γ::Float64

    B::Float64
    τs::Float64
    function Physics(ECS, CO₂_init, δT_init, Cd, γ)
        B = (3.48 * (60. * 60. * 24. * 365.25)) / ECS; # Transient Warming Parameter [K (W m^-2 s yr^-1)^-1] = [K (J yr^-1)^-1 m^2]
        τs = (Cd/B) * (B+γ)/γ
        return new(ECS, CO₂_init, δT_init, Cd, γ, B, τs)
    end
end

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
    extra_CO₂::Array{Float64,1}
end

function Economics(β, utility_discount_rate, reduce_cost, remove_cost, geoeng_cost, adapt_cost, reduce_init, remove_init, geoeng_init, adapt_init, baseline_emissions)
    return Economics(
        β::Float64,
        utility_discount_rate::Float64,
        reduce_cost::Float64,
        remove_cost::Float64,
        geoeng_cost::Float64,
        adapt_cost::Float64,
        reduce_init::Float64,
        remove_init::Float64,
        geoeng_init::Float64,
        adapt_init::Float64,
        baseline_emissions::Array{Float64,1},
        zeros(size(baseline_emissions))
    )
end

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
    init_zero_controls(t)

Create initial guess of zero climate controls.

See also: [`Controls`](@ref)
"""
function init_zero_controls(t::Array{Float64,1})
    c = Controls(
        nondim_linear(t)*0.,
        nondim_linear(t)*0.,
        nondim_linear(t)*0.,
        nondim_linear(t)*0.
    )
    return c
end


"""
    ClimateModel(name, domain, dt, present_year, economics, physics, controls)

Create instance of an extremely idealized integrated-assessment
climate model in which the climate response (`CO₂` and `temperature`) is a function
of physical input parameters (`CO₂_init`, `δT_init`, `ECS`, `ϵ`), economic input parameters
(`economics`), and climate control policies (`controls`) over some time frame (`domain`).

See also: [`Controls`](@ref), [`Economics`](@ref), [`CO₂`](@ref), [`δT`](@ref),
[`optimize!`](@ref)
"""
struct ClimateModel
    name::String
    domain::Array{Float64,1}
    dt::Float64
    present_year::Float64
    economics::Economics
    physics::Physics
    controls::Controls
end