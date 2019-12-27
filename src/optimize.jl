function ∇cost(model::ClimateModel; maxslope=maxslope)
    Δcontrol = 1.e-6
    domain_idx = (model.domain .>= model.present_year)
    ∇ = zeros((length(model.domain[domain_idx]), length(fieldnames(Controls))))
    objectives = [0., 0.]
    for (control_idx, controlname) in enumerate(fieldnames(Controls))
        for (idx, t_idx) in enumerate((1:length(model.domain))[domain_idx])
            perturb_model!(model, controlname, t_idx, Δcontrol)
            objectives[2] = discounted_total_cost_constrained(model, maxslope=maxslope)
            perturb_model!(model, controlname, t_idx, -Δcontrol)
            objectives[1] = discounted_total_cost_constrained(model, maxslope=maxslope)
            
            ∇[idx, control_idx] = diff(objectives)[1] / Δcontrol
        end
    end
    
    return ∇
end

function perturb_model!(model::ClimateModel, controlname::Symbol, time_idx::Int, Δcontrol::Float64)
    getfield(model.controls, controlname)[time_idx] += Δcontrol
end


function is_converged(objectives, tolerance)
    return abs(diff(objectives)[1]/objectives[1]) < tolerance
end

function optimize!(model::ClimateModel; tolerance=1.e-5, maxslope=1. /20.)
    domain_idx = (model.domain .>= model.present_year)
    
    ∇ = ∇cost(model, maxslope=maxslope)
    previous_update_vector = zeros(size(∇))
    
    iterations = 0
    
    objectives = [
        deepcopy(discounted_total_cost_constrained(model, maxslope=maxslope)),
        deepcopy(discounted_total_cost_constrained(model, maxslope=maxslope))+1.
    ]
    while !is_converged(objectives, tolerance)
        
        objectives[2] = objectives[1]
        
        ∇ = ∇cost(model, maxslope=maxslope)
        learning_rate = 1.e-4 / model.dt
        momentum_fraction = 0.90
        update_vector = (
            ∇ .* learning_rate + 
            previous_update_vector .* momentum_fraction
        )
        
        for (control_idx, controlname) in enumerate(fieldnames(Controls))
            getfield(model.controls, controlname)[domain_idx] .-= (
                update_vector[:,control_idx]
            )
        end
        
        previous_update_vector[:,:] = update_vector
        objectives[1] = discounted_total_cost_constrained(model, maxslope=maxslope)
        
#         print(
#             string(
#                 iterations,": ",
#                 round(sqrt(sum(∇.^2))*learning_rate, digits=5)," ",
#                 round(discounted_total_cost_constrained(model), digits=3), " ",
#                 round(discounted_total_cost(model), digits=3), " ",
#                 round(diff(objectives)[1]/objectives[1], digits=6),"\n"
#             )
#         )
        
        if iterations>5000
            break
        else
            iterations+=1
        end
    end
    print("Converged after $iterations iterations. ")
end