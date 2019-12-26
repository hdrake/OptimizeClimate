function ∇cost(model::ClimateModel)
    Δcontrol = 1.e-6
    domain_idx = (model.domain .>= model.present_year)
    ∇ = zeros((length(model.domain[domain_idx]), length(fieldnames(Controls))))
    for (control_idx, controlname) in enumerate(fieldnames(Controls))
        for (idx, t_idx) in enumerate((1:length(model.domain))[domain_idx])
            ∇[idx, control_idx] = ((
                    discounted_total_cost_constrained(
                        perturbed_model(model, controlname, t_idx, Δcontrol)
                    ) - discounted_total_cost_constrained(model)
                ) / Δcontrol
            )
        end
    end
    
    return ∇
end

function is_converged(∇, tolerance)
    return sum(∇.^2) < tolerance
end

function optimize!(model::ClimateModel, tolerance=1.e-4)
    domain_idx = (model.domain .>= model.present_year)
    
    ∇ = ∇cost(model)
    previous_update_vector = zeros(size(∇))
    
    iterations = 0
    while !is_converged(∇, tolerance)
        ∇ = ∇cost(model)
        learning_rate = 2.e-3
        momentum_fraction = 0.95
        update_vector = (
            ∇ .* learning_rate + 
            previous_update_vector .* momentum_fraction
        )
        
        for (control_idx, controlname) in enumerate(fieldnames(Controls))
            getfield(model.controls, controlname)[domain_idx] .-= (
                update_vector[:,control_idx]
            )
        end
        
        previous_update_vector = copy(update_vector)
        
        if iterations>5000
            break
        else
            iterations+=1
        end
    end
    print("Converged after $iterations iterations. ")
end