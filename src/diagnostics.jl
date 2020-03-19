
f(α::Array; p=2.) = α.^p # shape of individual cost functions

# Following RCP8.5 CO2e concentrations 
# Raw data at https://www.iiasa.ac.at/web-apps/tnt/RcpDb/dsd?Action=htmlpage&page=compare
#
# See below link for 2020 initial condition:
# https://www.eea.europa.eu/data-and-maps/indicators/atmospheric-greenhouse-gas-concentrations-6/assessment-1
function baseline_emissions(t::Array{Float64,1}, q0::Float64, t1::Float64, t2::Float64)
    t0 = t[1]
    Δt0 = t1 - t0
    Δt1 = t2 - t1
    q = zeros(size(t))
    increase_idx = (t .<= t1)
    decrease_idx = ((t .> t1) .& (t .<= t2))
    q[increase_idx] .= q0 * (1 .+ (t[increase_idx] .- t0)/Δt0)
    q[decrease_idx] .= 2. * q0 * (t2 .- t[decrease_idx])/Δt1
    q[t .> t2] .= 0.
    return q
end
baseline_emissions(t::Array{Float64,1}) = baseline_emissions(t::Array{Float64,1}, 15., 2100., 2140.)

effective_baseline_emissions(model::ClimateModel) = (
    model.physics.r * model.economics.baseline_emissions
)

controlled_emissions(model::ClimateModel) = (
    (1. .- model.controls.mitigate) .* model.economics.baseline_emissions
)

effective_emissions(model::ClimateModel) = (
    model.physics.r * model.economics.baseline_emissions .* (1. .- model.controls.mitigate) .-
    model.economics.baseline_emissions[1] .* model.controls.remove
)

function CO₂_baseline(model::ClimateModel)
    
    CO₂_baseline = zeros(size(model.domain))
    CO₂_baseline[model.domain .<= model.present_year] = CO₂(model)[model.domain .<= model.present_year]
    CO₂_baseline[model.domain .> model.present_year] = (
        CO₂_baseline[model.domain .== model.present_year] .+
        model.physics.r * cumsum(
            model.economics.baseline_emissions[model.domain .> model.present_year] .*
            model.dt
        )
    )
    
    return CO₂_baseline
end

CO₂(model::ClimateModel) = (
    model.physics.CO₂_init .+ (
        model.physics.r * cumsum(model.economics.baseline_emissions .* (1. .- model.controls.mitigate) .*
            model.dt) .-
        cumsum(model.economics.baseline_emissions[1] .* model.controls.remove .*
            model.dt)
    )
);

FCO₂_baseline(model::ClimateModel) = (
    5.35 .* log.( (CO₂_baseline(model) .+ model.physics.r * model.economics.extra_CO₂)./ CO₂_baseline(model)[1])
    * (60. * 60. * 24. * 365.25) # (W m^-2 s yr^-1)
)

FCO₂(model::ClimateModel) = (
    (5.35 .* log.( (CO₂(model) .+ model.physics.r * model.economics.extra_CO₂)./ CO₂(model)[1]) -
        8.5*model.controls.geoeng)
    * (60. * 60. * 24. * 365.25)
)

FCO₂_no_geoeng(model::ClimateModel) = (
    5.35 .* log.( (CO₂(model) .+ model.physics.r * model.economics.extra_CO₂)./ CO₂(model)[1])
    * (60. * 60. * 24. * 365.25)
)

δT_baseline(model::ClimateModel) = (
    model.physics.δT_init .+
    (FCO₂_baseline(model) .+ model.physics.κ * 
        (
            (model.physics.τd * model.physics.B)^(-1) .*
            exp.( - model.domain / model.physics.τd ) .*
            cumsum(
                exp.( model.domain / model.physics.τd ) .*
                FCO₂_baseline(model)
                .* model.dt
            )
        )
    ) .* (model.physics.B + model.physics.κ)^-1
)

δT_no_geoeng(model::ClimateModel) = (
    model.physics.δT_init .+
    (FCO₂_no_geoeng(model) .+ model.physics.κ * 
        (
            (model.physics.τd * model.physics.B)^(-1) .*
            exp.( - (model.domain .- model.domain[1]) / model.physics.τd ) .*
            cumsum(
                exp.( (model.domain .- model.domain[1]) / model.physics.τd ) .*
                FCO₂_no_geoeng(model)
                .* model.dt
            )
        )
    ) .* (model.physics.B + model.physics.κ)^-1
)

δT(model::ClimateModel) = ((
        model.physics.δT_init .+
        (FCO₂(model) .+ model.physics.κ * 
            (
                (model.physics.τd * model.physics.B)^(-1) .*
                exp.( - (model.domain .- model.domain[1]) / model.physics.τd ) .*
                cumsum(
                    exp.( (model.domain .- model.domain[1]) / model.physics.τd ) .*
                    FCO₂(model)
                    .* model.dt
                )
            )
        ) .* (model.physics.B + model.physics.κ)^-1
    )
)

function discounting(model::ClimateModel)
    discount = (1. .+ model.economics.utility_discount_rate) .^ (-(t .- model.present_year))
    discount[t .< model.present_year] .= 0.
    
    return discount
end
    
damage_cost_baseline(model::ClimateModel) = (
    model.economics.β .* model.economics.GWP .* δT_baseline(model).^2
)

discounted_damage_cost_baseline(model::ClimateModel) = (
    damage_cost_baseline(model) .* discounting(model)
)

damage_cost(model::ClimateModel) = (
    (1. .- model.controls.adapt) .*
    model.economics.β .* model.economics.GWP .* δT(model).^2
)

discounted_damage_cost(model::ClimateModel) = (
    damage_cost(model) .* discounting(model)
)

control_cost(model::ClimateModel) = (
    model.economics.mitigate_cost .* f(model.controls.mitigate) .+
    model.economics.remove_cost .* f(model.controls.remove) .+
    model.economics.geoeng_cost .* model.economics.GWP.*
    f(model.controls.geoeng) .+
    model.economics.adapt_cost .* f(model.controls.adapt)
)

discounted_control_cost(model::ClimateModel) = (
    control_cost(model) .* discounting(model)
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

mAtmos = 5.e18 # kg
tCO2_to_ppm(tCO2) = tCO2 / (mAtmos/1.e3) * 1.e6 
GtCO2_to_ppm(GtCO2) = GtCO2 * (1.e9) / (mAtmos/1.e3) * 1.e6

ppm_to_tCO2(ppm) = ppm / 1.e6 * (mAtmos/1.e3)
ppm_to_GtCO2(ppm) = ppm / 1.e6 / 1.e9 * (mAtmos/1.e3)

function extra_ton(model::ClimateModel, year::Float64)
    
    econ = model.economics
    
    year_idx = argmin(abs.(model.domain .- year))
    
    extra_CO₂ = zeros(size(model.domain))
    extra_CO₂[year_idx:end] .= tCO2_to_ppm(1.)
    
    new_economics = Economics(
        econ.β, econ.utility_discount_rate,
        econ.mitigate_cost, econ.remove_cost,
        econ.geoeng_cost, econ.adapt_cost,
        econ.mitigate_init, econ.remove_init, econ.geoeng_init, econ.adapt_init,
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