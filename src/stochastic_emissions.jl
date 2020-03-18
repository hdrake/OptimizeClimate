function add_emissions_bump(model::ClimateModel, Δt::Float64, Δq::Float64)

    present_year = deepcopy(model.present_year) + Δt
    present_idx = deepcopy(argmin(abs.(model.domain .- present_year)))
    name = string(Int64(round(present_year)));
    
    controls = Controls(
        deepcopy(model.controls.mitigate),
        deepcopy(model.controls.remove),
        deepcopy(model.controls.geoeng),
        deepcopy(model.controls.adapt)
    )
    
    future = (model.domain .>= model.present_year)
    near_future = future .& (model.domain .<= model.present_year + Δt)
    near_future1 = near_future .& (model.domain .< model.present_year + Δt/2)
    near_future2 = near_future .& (model.domain .>= model.present_year + Δt/2)
    
    new_emissions = deepcopy(model.economics.baseline_emissions)
    new_emissions[near_future1] .+= (
        Δq * (model.domain .- model.present_year) / (Δt/2.)
    )[near_future1]
    new_emissions[near_future2] .+= (
        Δq * (1. .- (model.domain .- (model.present_year .+ Δt/2.)) / (Δt/2.))
    )[near_future2]
    
    econ = model.economics
    economics = Economics(
        econ.β, econ.utility_discount_rate,
        econ.mitigate_cost, econ.remove_cost, econ.geoeng_cost, econ.adapt_cost,
        econ.mitigate_init, econ.remove_init, econ.geoeng_init, econ.adapt_init,
        new_emissions
    )
    model = ClimateModel(
        name, model.domain, model.dt, present_year, economics, model.physics, controls,
    );
    return model
end