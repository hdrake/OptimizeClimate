function ∇cost(model::ClimateModel)
    Δcontrol = 1.e-6
    ∇ = zeros((length(model.domain), length(fieldnames(Controls))))
    for (control_idx, controlname) in enumerate(fieldnames(Controls))
        for time_idx in 1:length(model.domain)
            ∇[time_idx,control_idx] = ((
                    discounted_total_cost_constrained(perturbed_model(
                            model, controlname, time_idx, Δcontrol)
                    ) -
                    discounted_total_cost_constrained(model)
                ) / Δcontrol
            )
        end
    end
    
    return ∇
end

function is_converged(∇, tolerance)
    return sum(∇.^2) < tolerance
end

function optimize!(model::ClimateModel, tolerance=1.e-6)
    ∇ = ∇cost(model)
    previous_update_vector = zeros(size(∇))
    
    iterations = 0
    while !is_converged(∇, tolerance)
        ∇ = ∇cost(model)
        learning_rate = 1.e-3
        momentum_fraction = 0.9
        update_vector = (
            ∇ .* learning_rate + 
            previous_update_vector .* momentum_fraction
        )
        
        for (control_idx, controlname) in enumerate(fieldnames(Controls))
            getfield(model.controls, controlname) .-= (
                update_vector[:,control_idx]
            )
        end
        
        previous_update_vector = copy(update_vector)

        if iterations>200
            break
        else
            iterations+=1
        end
        
    end
    print("Converged after $iterations iterations.\n")
end