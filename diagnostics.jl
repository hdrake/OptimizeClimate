
f_low(α::Array) = (α ./ (1. .+ α)).^2 # shape of individual cost functions
f_med(α::Array) = α.^2 # shape of individual cost functions
f_high(α::Array) = (α ./ (1. .- α)).^2 # shape of individual cost functions

function baseline_emissions(t::Array{Float64,1})
    q = zeros(size(t))
    q[(t.<=2060)] = 5. .* ones(size(t[(t.<=2060)])); # emissions scenario
    q[(t.>2060) .& (t.<=2100)] .= 5. * (40. .- (t[(t.>2060) .& (t.<=2100)] .- 2060.))/40.
    q[(t.>2100)] .= 0.
    return q
end

CO₂_baseline(model::ClimateModel) = model.CO₂_init .+ cumsum(baseline_emissions(model.domain))

CO₂(model::ClimateModel) = (
    CO₂_baseline(model) .-
    cumsum(baseline_emissions(model.domain) .* model.controls.reduce) .-
    cumsum(baseline_emissions(model.domain)[1] .* model.controls.remove)
);

δT_baseline(model::ClimateModel) = (
        model.δT_pre .+
        model.ϵ .* log.( CO₂_baseline(model)./ model.CO₂_init )
)

δT_no_geoeng(model::ClimateModel) = (
        model.δT_pre .+
        model.ϵ .* log.( CO₂(model)./ model.CO₂_init )
)

δT(model::ClimateModel) = (
    (
        model.δT_pre .+
        model.ϵ .* log.( CO₂(model)./ model.CO₂_init )
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

net_cost(model::ClimateModel) = (
    damage_cost(model) .+ control_cost(model)
)

discounted_net_cost(model::ClimateModel) = (
    (damage_cost(model) .+ control_cost(model)) .* discounting(model)
)

total_cost(model::ClimateModel) = (
    sum(net_cost(model))
)

discounted_total_cost(model::ClimateModel) = (
    sum(net_cost(model) .* discounting(model))
)

total_cost_constrained(model::ClimateModel) = (
    total_cost(model) +
    200. * (
        sum(diff(model.controls.reduce).^2) +
        sum(diff(model.controls.remove).^2) +
        sum(diff(model.controls.geoeng).^2) +
        sum(diff(model.controls.adapt).^2)
    ) +
    500. * (
        (model.controls.reduce[1] - model.economics.reduce_init).^2 .+
        (model.controls.remove[1] - model.economics.remove_init).^2 .+
        (model.controls.geoeng[1] - model.economics.geoeng_init).^2 .+
        (model.controls.adapt[1] - model.economics.adapt_init).^2
    )
)

discounted_total_cost_constrained(model::ClimateModel) = (
    discounted_total_cost(model) +
    200. * (
        sum(diff(model.controls.reduce).^2) +
        sum(diff(model.controls.remove).^2) +
        sum(diff(model.controls.geoeng).^2) +
        sum(diff(model.controls.adapt).^2)
    ) +
    500. * (
        (model.controls.reduce[1] - model.economics.reduce_init).^2 .+
        (model.controls.remove[1] - model.economics.remove_init).^2 .+
        (model.controls.geoeng[1] - model.economics.geoeng_init).^2 .+
        (model.controls.adapt[1] - model.economics.adapt_init).^2
    )
)

function perturbed_model(model::ClimateModel, controlname::Symbol, time_idx::Int, Δcontrol::Float64)
    perturbed_controls = Controls(
        deepcopy(model.controls.reduce),
        deepcopy(model.controls.remove),
        deepcopy(model.controls.geoeng),
        deepcopy(model.controls.adapt)
    )
    getfield(perturbed_controls, controlname)[time_idx] = getfield(perturbed_controls, controlname)[time_idx] + Δcontrol
    return ClimateModel(
        model.name, model.ECS, model.domain, perturbed_controls,
        model.economics, model.present_year, model.CO₂_init, model.δT_pre
    )
end