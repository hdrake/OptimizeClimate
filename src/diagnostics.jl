
f_low(α::Array) = (α ./ (1. .+ α)).^2 # shape of individual cost functions
f_med(α::Array) = α.^2 # shape of individual cost functions
f_high(α::Array) = (α ./ (1. .- α)).^2 # shape of individual cost functions

function baseline_emissions(t::Array{Float64,1}, q0::Float64, t0::Float64, Δt::Float64)
    q = zeros(size(t))
    q[t .<= t0] = q0 .* ones(size(t[t .<= t0])); # emissions scenario
    q[(t .> t0) .& (t .<= (t0+Δt))] .= q0 * (Δt .- (t[(t .> t0) .& (t .<= (t0 + Δt))] .- t0))/Δt
    q[t .> (t0 + Δt)] .= 0.
    return q
end

baseline_emissions(t::Array{Float64,1}) = baseline_emissions(t::Array{Float64,1}, 5., 2060., 40.)

effective_emissions(model::ClimateModel) = (
    model.economics.baseline_emissions .* (1. .- model.controls.reduce) .-
    model.economics.baseline_emissions[1] .* model.controls.remove
)

function CO₂_baseline(model::ClimateModel)
    
    CO₂_baseline = zeros(size(model.domain))
    CO₂_baseline[model.domain .<= model.present_year] = CO₂(model)[model.domain .<= model.present_year]
    CO₂_baseline[model.domain .> model.present_year] = (
        CO₂_baseline[model.domain .== model.present_year] .+
        cumsum(
            model.economics.baseline_emissions[model.domain .> model.present_year] .*
            model.dt
        )
    )
    
    return CO₂_baseline
end

CO₂(model::ClimateModel) = (
    model.physics.CO₂_init .+
    cumsum(model.economics.baseline_emissions .* (1. .- model.controls.reduce) .*
        model.dt) .-
    cumsum(model.economics.baseline_emissions[1] .* model.controls.remove .*
        model.dt)
);

FCO₂_baseline(model::ClimateModel) = (
    (5.35 .* log.( CO₂_baseline(model)./ CO₂_baseline(model)[1]))
    * (60. * 60. * 24. * 365.25) # (W m^-2 s yr^-1)
)

FCO₂(model::ClimateModel) = (
    (5.35 .* log.( CO₂(model)./ CO₂(model)[1]))
    * (60. * 60. * 24. * 365.25)
)

δT_baseline(model::ClimateModel) = (
    model.physics.δT_init .+
    (FCO₂_baseline(model) .+ model.physics.γ * 
        (
            (model.physics.τs * model.physics.B)^(-1) .*
            exp.( - model.domain / model.physics.τs ) .*
            cumsum(
                exp.( model.domain / model.physics.τs ) .*
                FCO₂_baseline(model)
                .* model.dt
            )
        )
    ) .* (model.physics.B + model.physics.γ)^-1
)

δT_no_geoeng(model::ClimateModel) = (
    model.physics.δT_init .+
    (FCO₂(model) .+ model.physics.γ * 
        (
            (model.physics.τs * model.physics.B)^(-1) .*
            exp.( - (model.domain .- model.domain[1]) / model.physics.τs ) .*
            cumsum(
                exp.( (model.domain .- model.domain[1]) / model.physics.τs ) .*
                FCO₂(model)
                .* model.dt
            )
        )
    ) .* (model.physics.B + model.physics.γ)^-1
)

δT(model::ClimateModel) = ((
        model.physics.δT_init .+
        (FCO₂(model) .+ model.physics.γ * 
            (
                (model.physics.τs * model.physics.B)^(-1) .*
                exp.( - (model.domain .- model.domain[1]) / model.physics.τs ) .*
                cumsum(
                    exp.( (model.domain .- model.domain[1]) / model.physics.τs ) .*
                    FCO₂(model)
                    .* model.dt
                )
            )
        ) .* (model.physics.B + model.physics.γ)^-1
    ) .* (1. .- model.controls.geoeng)
)

function discounting(model::ClimateModel)
    discount = (1. .+ model.economics.utility_discount_rate) .^ (-(t .- model.present_year))
    discount[t .< model.present_year] .= 0.
    
    return discount
end
    
damage_cost_baseline(model::ClimateModel) = (
    model.economics.β .* δT_baseline(model).^2
)

discounted_damage_cost_baseline(model::ClimateModel) = (
    model.economics.β .* δT_baseline(model).^2 .*
    discounting(model)
)

damage_cost(model::ClimateModel) = (
    (1. .- model.controls.adapt) .*
    model.economics.β .* δT(model).^2
)

discounted_damage_cost(model::ClimateModel) = (
    (1. .- model.controls.adapt) .*
    model.economics.β .* δT(model).^2 .*
    discounting(model)
)

control_cost(model::ClimateModel) = (
    model.economics.reduce_cost .* f_med(model.controls.reduce) .+
    model.economics.remove_cost .* f_med(model.controls.remove) .+
    model.economics.geoeng_cost .* f_med(model.controls.geoeng) .+
    model.economics.adapt_cost .* f_med(model.controls.adapt)
)

discounted_control_cost(model::ClimateModel) = (
    (
        model.economics.reduce_cost .* f_med(model.controls.reduce) .+
        model.economics.remove_cost .* f_med(model.controls.remove) .+
        model.economics.geoeng_cost .* f_med(model.controls.geoeng) .+
        model.economics.adapt_cost .* f_med(model.controls.adapt)
    ) .* discounting(model)
)

discounted_total_damage_cost(model::ClimateModel) = (
    sum(discounted_damage_cost(model) .* model.dt)
)

net_cost(model::ClimateModel) = (
    damage_cost(model) .+ control_cost(model)
)

discounted_net_cost(model::ClimateModel) = (
    (damage_cost(model) .+ control_cost(model)) .* discounting(model)
)

total_cost(model::ClimateModel) = (
    sum(net_cost(model) .* model.dt)
)

discounted_total_cost(model::ClimateModel) = (
    sum(net_cost(model) .* discounting(model)  .* model.dt)
)

function extra_ton(model::ClimateModel, year::Float64)
    
    econ = model.economics
    
    year_idx = argmin(abs.(model.domain .- year))
    
    extra_CO₂ = zeros(size(model.domain))
    extra_CO₂[year_idx:end] .= 1. /(1.25e10)
    
    new_economics = Economics(
        econ.β, econ.utility_discount_rate,
        econ.reduce_cost, econ.remove_cost,
        econ.geoeng_cost, econ.adapt_cost,
        0., 0., 0., 0.,
        econ.baseline_emissions,
        extra_CO₂
    );
    
    return ClimateModel(
        model.name, model.domain, model.dt, model.present_year,
        new_economics, model.physics, model.controls
    )
end

extra_ton(model::ClimateModel) = extra_ton(model::ClimateModel, model.domain[1])

SCC(model::ClimateModel, year::Float64) = round((
     discounted_total_cost(extra_ton(model, year)) -
     discounted_total_cost(model)
)*1.e12, digits=2)

SCC(model::ClimateModel) = SCC(model::ClimateModel, model.domain[1])

GWP = 100. # global world product (trillion $ / year)

β = 0.02*GWP/(3.0)^2 # damages (trillion USD / year / celsius^2)
utility_discount_rate = 0.0 # ρ (relative low value from Stern review)

# Control technology cost scales, as fraction of GWP (cost scale is for full deployment, α=1.)
reduce_cost = 0.01*GWP;
remove_cost = 0.02*GWP;
geoeng_cost = 0.05*GWP;
adapt_cost = 0.03*GWP;

"""
    Economics()

Create data structure for economic input parameters for `ClimateModel` struct with default values.

Default parameters are:
- `β`= 0.222 × 10^12 USD / (°C)^2
- `utility_discount_rate` = 0.025 (between Stern review median value of 1.4% and Nordhaus values)
- `reduce_cost` = 1. × 10^12 USD
- `remove_cost` = 2. × 10^12 USD
- `geoeng_cost` = 5. × 10^12 USD
- `adapt_cost` = 3. × 10^12 USD
- `[control]_init` = 0.
- `baseline_emissions` = baseline_emissions(t::Array{Float64,1}, 5., 2080., 40.)

The default baseline emissions scenario corresponds to flat emissions of 5 ppm / year
from 2020 to 2080 and linearly decreasing from 5 ppm / year in 2080 to 0 ppm / year in 2120.

See also: [`ClimateModel`](@ref), [`baseline_emissions`](@ref)
"""
Economics(t) = Economics(
    GWP, β,
    reduce_cost, remove_cost, geoeng_cost, adapt_cost,
    0., 0., 0., 0., # Assumed initial condition of zero control deployments in 2020
    baseline_emissions(t, 5., 2080., 40.)
)