function optimize_controls!(
        model::ClimateModel;
        obj_option = "temp", temp_goal = 2., budget=10., expenditure = 0.5,
        max_deployment = Dict("mitigate"=>1., "remove"=>1., "geoeng"=>1., "adapt"=>1.),
        maxslope = Dict("mitigate"=>1. /20., "remove"=>1. /40., "geoeng"=>1. /20., "adapt"=>1. /20.),
        temp_final = nothing,
        start_deployment = Dict(
            "mitigate"=>model.domain[1],
            "remove"=>model.domain[1]+20,
            "geoeng"=> model.domain[1]+40,
            "adapt"=>model.domain[1]
        ),
        cost_exponent = 2
    )
    
    if temp_final == nothing
        temp_final = temp_goal
    elseif temp_final >= temp_goal
        temp_final = temp_goal
    end
    
    for (key, item) in start_deployment
        if item == nothing
            start_deployment[key] = model.present_year
        end
    end
    
    model_optimizer = Model(optimizer_with_attributes(Ipopt.Optimizer))#, "print_level" => 0))

    f_JuMP(α) = α^cost_exponent
    register(model_optimizer, :f_JuMP, 1, f_JuMP, autodiff=true)

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
            0. <= M[1:N] <= max_deployment["mitigate"]  # emissions reductions
            0. <= R[1:N] <= max_deployment["remove"]  # negative emissions
            0. <= G[1:N] <= max_deployment["geoeng"]  # geoengineering
            0. <= A[1:N] <= max_deployment["adapt"]  # adapt
    end)

    control_vars = Dict(
        "mitigate" => M,
        "remove" => R,
        "geoeng" => G,
        "adapt" => A
    )
    controls = Dict(
        "mitigate" => model.controls.mitigate,
        "remove" => model.controls.remove,
        "geoeng" => model.controls.geoeng,
        "adapt" => model.controls.adapt
    )
    
    domain_idx = (model.domain .>= model.present_year)
    
    M₀ = model.economics.mitigate_init
    R₀ = model.economics.remove_init
    G₀ = model.economics.geoeng_init
    A₀ = model.economics.adapt_init
    
    control_inits = Dict(
        "mitigate" => M₀,
        "remove" => R₀,
        "geoeng" => G₀,
        "adapt" => A₀
    )
    
    for (key, control) in control_vars
        fix(control_vars[key][1], control_inits[key]; force = true)

        for idx in 2:N
            if idx <= length(model.domain[.~domain_idx])
                fix(control_vars[key][idx], controls[key][idx]; force = true)
            else
                if model.domain[idx] < start_deployment[key]
                    fix(control_vars[key][idx], control_inits[key]; force = true)
                end
            end
        end
    end

    # add integral function as a new variable defined by first order finite differences
    @variable(model_optimizer, cumsum_qMR[1:N]);
    for i in 1:N-1
        @constraint(
            model_optimizer, cumsum_qMR[i+1] - cumsum_qMR[i] ==
            (model.dt * (q[i+1] * (1. - M[i+1]) - q[1] * R[i+1]))
        )
    end
    @constraint(
        model_optimizer, cumsum_qMR[1] == 
        (model.dt * (q[1] * (1. - M[1]) - q[1] * R[1]))
    );
    
    # add temperature kernel as new variable defined by first order finite difference
    @variable(model_optimizer, cumsum_KFdt[1:N]);
    for i in 1:N-1
        @NLconstraint(
            model_optimizer, cumsum_KFdt[i+1] - cumsum_KFdt[i] ==
            (
                model.dt *
                exp( - model.domain[i+1] / model.physics.τs ) *
                5.35 * log_JuMP(
                    (model.physics.CO₂_init + cumsum_qMR[i+1]) /
                    (model.physics.CO₂_init + cumsum_qMR[1])
                ) * (60. * 60. * 24. * 365.25)
            )
        )
    end
    @constraint(
        model_optimizer, cumsum_KFdt[1] == 0.
    );

    # Add constraint of rate of changes
    if typeof(maxslope) == Float64
        @variables(model_optimizer, begin
                -maxslope <= dMdt[1:N-1] <= maxslope
                -maxslope <= dRdt[1:N-1] <= maxslope
                -maxslope <= dGdt[1:N-1] <= maxslope
                -maxslope <= dAdt[1:N-1] <= maxslope
        end);
    elseif typeof(maxslope) == Dict{String,Float64}
        @variables(model_optimizer, begin
                -maxslope["mitigate"] <= dMdt[1:N-1] <= maxslope["mitigate"]
                -maxslope["remove"] <= dRdt[1:N-1] <= maxslope["remove"]
                -maxslope["geoeng"] <= dGdt[1:N-1] <= maxslope["geoeng"]
                -maxslope["adapt"] <= dAdt[1:N-1] <= maxslope["adapt"]
        end);
    end

    for i in 1:N-1
        @constraint(model_optimizer, dMdt[i] == (M[i+1] - M[i]) / model.dt)
        @constraint(model_optimizer, dRdt[i] == (R[i+1] - R[i]) / model.dt)
        @constraint(model_optimizer, dGdt[i] == (G[i+1] - G[i]) / model.dt)
        @constraint(model_optimizer, dAdt[i] == (A[i+1] - A[i]) / model.dt)
    end

    if obj_option == "net_cost"
        # objective function to minimize
        @NLobjective(model_optimizer, Min, 
            sum(
                (
                    (1 - A[i]) * model.economics.β *
                    ((model.physics.δT_init + 
                        (
                            5.35 * log_JuMP(
                                (model.physics.CO₂_init + cumsum_qMR[i]) /
                                (model.physics.CO₂_init + cumsum_qMR[1])
                            ) * (60. * 60. * 24. * 365.25) + model.physics.γ *
                            (model.physics.τs * model.physics.B)^(-1) *
                            exp( ( model.domain[i] / model.physics.τs )) *
                            cumsum_KFdt[i]
                        ) * (model.physics.B + model.physics.γ)^-1
                    ) * (1. - G[i])
                    )^2 +
                    model.economics.mitigate_cost * f_JuMP(M[i]) +
                    model.economics.remove_cost * f_JuMP(R[i]) +
                    model.economics.geoeng_cost * f_JuMP(G[i]) +
                    model.economics.adapt_cost * f_JuMP(A[i])
                ) *
                discounting_JuMP(model.domain[i]) *
                model.dt
            for i=1:N)
        )
        
    elseif obj_option == "temp"
        @NLobjective(model_optimizer, Min,
            sum(
                (
                    model.economics.mitigate_cost * f_JuMP(M[i]) +
                    model.economics.remove_cost * f_JuMP(R[i]) +
                    model.economics.geoeng_cost * f_JuMP(G[i]) +
                    model.economics.adapt_cost * f_JuMP(A[i])
                ) *
                discounting_JuMP(model.domain[i]) *
                model.dt
            for i=1:N)
        )
        
        for i in 1:N-1
            @NLconstraint(model_optimizer,
                (1 - A[i]) * model.economics.β *
                ((model.physics.δT_init + 
                    (
                        5.35 * log_JuMP(
                            (model.physics.CO₂_init + cumsum_qMR[i]) /
                            (model.physics.CO₂_init + cumsum_qMR[1])
                        ) * (60. * 60. * 24. * 365.25) + model.physics.γ *
                        (model.physics.τs * model.physics.B)^(-1) *
                        exp( ( model.domain[i] / model.physics.τs )) *
                        cumsum_KFdt[i]
                    ) * (model.physics.B + model.physics.γ)^-1
                ) * (1. - G[i])
                )^2 <= (model.economics.β * temp_goal^2)
            )
        end
        i=N
        @NLconstraint(model_optimizer,
            (1 - A[i]) * model.economics.β *
            ((model.physics.δT_init + 
                (
                    5.35 * log_JuMP(
                        (model.physics.CO₂_init + cumsum_qMR[i]) /
                        (model.physics.CO₂_init + cumsum_qMR[1])
                    ) * (60. * 60. * 24. * 365.25) + model.physics.γ *
                    (model.physics.τs * model.physics.B)^(-1) *
                    exp( ( model.domain[i] / model.physics.τs )) *
                    cumsum_KFdt[i]
                ) * (model.physics.B + model.physics.γ)^-1
            ) * (1. - G[i])
            )^2 <= (model.economics.β * temp_final^2)
        )

    elseif obj_option == "budget"
        @NLobjective(model_optimizer, Min,
            sum(
                (1 - A[i]) * model.economics.β *
                (
                    (model.physics.δT_init + 
                        (
                            5.35 * log_JuMP(
                                (model.physics.CO₂_init + cumsum_qMR[i]) /
                                (model.physics.CO₂_init + cumsum_qMR[1])
                            ) * (60. * 60. * 24. * 365.25) + model.physics.γ *
                            (model.physics.τs * model.physics.B)^(-1) *
                            exp( ( model.domain[i] / model.physics.τs )) *
                            cumsum_KFdt[i]
                        ) * (model.physics.B + model.physics.γ)^-1
                    ) * (1. - G[i])
                )^2 *
                discounting_JuMP(model.domain[i]) *
                model.dt
            for i=1:N)
        )
        
        @NLconstraint(model_optimizer,
            sum(
                (
                    model.economics.mitigate_cost * f_JuMP(M[i]) +
                    model.economics.remove_cost * f_JuMP(R[i]) +
                    model.economics.geoeng_cost * f_JuMP(G[i]) +
                    model.economics.adapt_cost * f_JuMP(A[i])
                ) *
                discounting_JuMP(model.domain[i]) *
                model.dt
            for i=1:N) <= budget
        )
        
    elseif obj_option == "expenditure"
        @NLobjective(model_optimizer, Min,
            sum(
                (1 - χ[i]) * model.economics.β *
                (
                    (model.physics.δT_init + 
                        (
                            5.35 * log_JuMP(
                                (model.physics.CO₂_init + cumsum_qMR[i]) /
                                (model.physics.CO₂_init + cumsum_qMR[1])
                            ) * (60. * 60. * 24. * 365.25) + model.physics.γ *
                            (model.physics.τs * model.physics.B)^(-1) *
                            exp( ( model.domain[i] / model.physics.τs )) *
                            cumsum_KFdt[i]
                        ) * (model.physics.B + model.physics.γ)^-1
                    ) * (1. - λ[i])
                )^2 *
                discounting_JuMP(model.domain[i]) *
                model.dt
            for i=1:N)
        )
        
        for i in 1:N
            @NLconstraint(model_optimizer,
                (
                    model.economics.mitigate_cost * f_med_JuMP(M[i]) +
                    model.economics.remove_cost * f_med_JuMP(R[i]) +
                    model.economics.geoeng_cost * f_med_JuMP(G[i]) +
                    model.economics.adapt_cost * f_med_JuMP(A[i])
                ) <= expenditure
            )
        end
    end
    
    optimize!(model_optimizer)
    
    getfield(model.controls, :mitigate)[domain_idx] = value.(M)[domain_idx]
    getfield(model.controls, :remove)[domain_idx] = value.(R)[domain_idx]
    getfield(model.controls, :geoeng)[domain_idx] = value.(G)[domain_idx]
    getfield(model.controls, :adapt)[domain_idx] = value.(A)[domain_idx]
    
end

function step_forward(model::ClimateModel, Δt::Float64, q0::Float64, t0::Float64, Δt0::Float64)

    present_year = deepcopy(model.present_year) + Δt
    present_idx = deepcopy(argmin(abs.(model.domain .- present_year)))
    name = string(Int64(round(present_year)));
    
    controls = Controls(
        deepcopy(model.controls.mitigate),
        deepcopy(model.controls.remove),
        deepcopy(model.controls.geoeng),
        deepcopy(model.controls.adapt)
    )
    
    new_emissions = zeros(size(model.domain))
    new_emissions[model.domain .< model.present_year] = (
        deepcopy(model.economics.baseline_emissions)[model.domain .< model.present_year]
    )
    new_emissions[model.domain .>= model.present_year] = (
        deepcopy(baseline_emissions(model.domain, q0, t0, Δt0))[model.domain .>= model.present_year]
    )
    
    economics = Economics(
        β, utility_discount_rate,
        mitigate_cost, remove_cost, geoeng_cost, adapt_cost,
        0., 0., 0., 0.,
        new_emissions
    )
    model = ClimateModel(
        name, model.domain, model.dt, present_year, economics, model.physics, controls,
    );
    return model
end

step_forward(model::ClimateModel, Δt::Float64, q0::Float64) = step_forward(model::ClimateModel, Δt::Float64, q0::Float64, 2080., 40.)
step_forward(model::ClimateModel, Δt::Float64) = step_forward(model::ClimateModel, Δt::Float64, 5., 2080., 40.)