function optimize_controls!(
        model::ClimateModel; maxslope = 1. /30.,
        obj_option="net_cost", temp_goal = 2., budget=10., expenditure = 0.5
    )
    model_optimizer = Model(with_optimizer(Ipopt.Optimizer, print_level=0))

    f_med_JuMP(α) = α^2
    register(model_optimizer, :f_med_JuMP, 1, f_med_JuMP, autodiff=true)

    function log_JuMP(x)
        if x == 0
            return -1000.0
        else
            return log(x)
        end
    end

    register(model_optimizer, :log_JuMP, 1, log_JuMP, autodiff=true)

    function discounting_JuMP(t)
        if t < model.present_year
            return 0.
        else
            return (
                (1. .+ model.economics.utility_discount_rate) .^
                (-(t .- model.present_year))
            )
        end
    end

    register(model_optimizer, :discounting_JuMP, 1, discounting_JuMP, autodiff=true)

    q = model.economics.baseline_emissions
    N = length(model.domain)

    # constraints on control variables
    @variables(model_optimizer, begin
            0. <= ϕ[1:N] <= 1.  # negative emissions
            0. <= φ[1:N] <= 1.  # emissions reductions
            0. <= χ[1:N] <= 1.  # geoengineering
            0. <= λ[1:N] <= 1.  # adapt
    end)

    ϕ₀ = model.economics.remove_init
    φ₀ = model.economics.reduce_init
    λ₀ = model.economics.geoeng_init
    χ₀ = model.economics.adapt_init

    fix(ϕ[1], ϕ₀; force = true)
    fix(φ[1], φ₀; force = true)
    fix(λ[1], λ₀; force = true)
    fix(χ[1], χ₀; force = true)
    
    domain_idx = (model.domain .>= model.present_year)
    
    for idx in 2:length(model.domain[.~domain_idx])
        fix(ϕ[idx], model.controls.remove[idx]; force = true)
        fix(φ[idx], model.controls.reduce[idx]; force = true)
        fix(λ[idx], model.controls.geoeng[idx]; force = true)
        fix(χ[idx], model.controls.adapt[idx]; force = true)
    end

    # add integral function as a new variable defined by first order finite differences
    @variable(model_optimizer, cumsum_qφϕ[1:N]);
    for i in 1:N-1
        @constraint(
            model_optimizer, cumsum_qφϕ[i+1] - cumsum_qφϕ[i] ==
            (model.dt * (q[i+1] * (1. - φ[i+1]) - q[1] * ϕ[i+1]))
        )
    end
    @constraint(
        model_optimizer, cumsum_qφϕ[1] == 
        (model.dt * (q[1] * (1. - φ[1]) - q[1] * ϕ[1]))
    );

    # Add constraint of rate of changes
    @variables(model_optimizer, begin
            -maxslope <= dϕdt[1:N-1] <= maxslope
            -maxslope <= dφdt[1:N-1] <= maxslope
            -maxslope <= dλdt[1:N-1] <= maxslope
            -maxslope <= dχdt[1:N-1] <= maxslope
    end);

    for i in 1:N-1
        @constraint(model_optimizer, dϕdt[i] == (ϕ[i+1] - ϕ[i]) / model.dt)
        @constraint(model_optimizer, dφdt[i] == (φ[i+1] - φ[i]) / model.dt)
        @constraint(model_optimizer, dλdt[i] == (λ[i+1] - λ[i]) / model.dt)
        @constraint(model_optimizer, dχdt[i] == (χ[i+1] - χ[i]) / model.dt)
    end

    if obj_option == "net_cost"
        # objective function to minimize
        @NLobjective(model_optimizer, Min, 
            sum(
                (
                    (1 - χ[i]) * model.economics.β *
                    (model.δT_init + model.ϵ * log_JuMP(
                        (model.CO₂_init + cumsum_qφϕ[i]) /
                        (model.CO₂_init + cumsum_qφϕ[1])
                    ))^2 *
                    (1 - λ[i])^2 +
                    model.economics.remove_cost * f_med_JuMP(ϕ[i]) +
                    model.economics.reduce_cost * f_med_JuMP(φ[i]) +
                    model.economics.geoeng_cost * f_med_JuMP(λ[i]) +
                    model.economics.adapt_cost * f_med_JuMP(χ[i])
                ) *
                discounting_JuMP(model.domain[i]) *
                model.dt
            for i=1:N)
        )
        
    elseif obj_option == "temp"
        @NLobjective(model_optimizer, Min,
            sum(
                (
                    model.economics.remove_cost * f_med_JuMP(ϕ[i]) +
                    model.economics.reduce_cost * f_med_JuMP(φ[i]) +
                    model.economics.geoeng_cost * f_med_JuMP(λ[i]) +
                    model.economics.adapt_cost * f_med_JuMP(χ[i])
                ) *
                discounting_JuMP(model.domain[i]) *
                model.dt
            for i=1:N)
        )
        
        temp_goal_idx = argmin(abs.(model.domain .- 2100.))
        
        for i in 1:N
            @NLconstraint(model_optimizer,
                (
                    (model.δT_init + model.ϵ * log_JuMP(
                        (model.CO₂_init + cumsum_qφϕ[i]) /
                        (model.CO₂_init + cumsum_qφϕ[1])
                    )) * (1 - λ[i])
                ) <= temp_goal
            )
        end

    elseif obj_option == "budget"
        @NLobjective(model_optimizer, Min,
            sum(
                (
                    (1 - χ[i]) * model.economics.β *
                    (model.δT_init + model.ϵ * log_JuMP(
                        (model.CO₂_init + cumsum_qφϕ[i]) /
                        (model.CO₂_init + cumsum_qφϕ[1])
                    ))^2 *
                    (1 - λ[i])^2
                ) *
                discounting_JuMP(model.domain[i]) *
                model.dt
            for i=1:N)
        )
        
        @NLconstraint(model_optimizer,
            sum(
                (
                    model.economics.remove_cost * f_med_JuMP(ϕ[i]) +
                    model.economics.reduce_cost * f_med_JuMP(φ[i]) +
                    model.economics.geoeng_cost * f_med_JuMP(λ[i]) +
                    model.economics.adapt_cost * f_med_JuMP(χ[i])
                ) *
                discounting_JuMP(model.domain[i]) *
                model.dt
            for i=1:N) <= budget
        )
        
    elseif obj_option == "expenditure"
        @NLobjective(model_optimizer, Min,
            sum(
                (
                    (1 - χ[i]) * model.economics.β *
                    (model.δT_init + model.ϵ * log_JuMP(
                        (model.CO₂_init + cumsum_qφϕ[i]) /
                        (model.CO₂_init + cumsum_qφϕ[1])
                    ))^2 *
                    (1 - λ[i])^2
                ) *
                discounting_JuMP(model.domain[i]) *
                model.dt
            for i=1:N)
        )
        
        for i in 1:N
            @NLconstraint(model_optimizer,
                (
                    model.economics.remove_cost * f_med_JuMP(ϕ[i]) +
                    model.economics.reduce_cost * f_med_JuMP(φ[i]) +
                    model.economics.geoeng_cost * f_med_JuMP(λ[i]) +
                    model.economics.adapt_cost * f_med_JuMP(χ[i])
                ) <= expenditure
            )
        end
    end
    
    optimize!(model_optimizer)
    print("Found optimal solution for model name: ", model.name)
    
    getfield(model.controls, :remove)[domain_idx] = value.(ϕ)[domain_idx]
    getfield(model.controls, :reduce)[domain_idx] = value.(φ)[domain_idx]
    getfield(model.controls, :geoeng)[domain_idx] = value.(λ)[domain_idx]
    getfield(model.controls, :adapt)[domain_idx] = value.(χ)[domain_idx]
    
end


function step_forward(model::ClimateModel, Δt::Float64, q0::Float64, t0::Float64, Δt0::Float64)

    present_year = deepcopy(model.present_year) + Δt
    present_idx = deepcopy(argmin(abs.(model.domain .- present_year)))
    name = string(Int64(round(present_year)));
    
    controls = Controls(
        deepcopy(model.controls.reduce),
        deepcopy(model.controls.remove),
        deepcopy(model.controls.geoeng),
        deepcopy(model.controls.adapt)
    )
    
    new_emissions = zeros(size(model.domain))
    new_emissions[model.domain .< model.present_year] = deepcopy(model.economics.baseline_emissions)[model.domain .< model.present_year]
    new_emissions[model.domain .>= model.present_year] = deepcopy(baseline_emissions(model.domain, q0, t0, Δt0))[model.domain .>= model.present_year]
    
    economics = Economics(
        β, utility_discount_rate,
        reduce_cost, remove_cost, geoeng_cost, adapt_cost,
        0., 0., 0., 0.,
        new_emissions
    )
    model = ClimateModel(
        name, model.ECS, model.domain, model.dt, controls, economics, present_year,
    );
    return model
end

step_forward(model::ClimateModel, Δt::Float64, q0::Float64) = step_forward(model::ClimateModel, Δt::Float64, q0::Float64, 2080., 40.)
step_forward(model::ClimateModel, Δt::Float64) = step_forward(model::ClimateModel, Δt::Float64, 5., 2080., 40.)