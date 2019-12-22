function ensemble_diagnostic_statistics(ensemble::Dict{String, ClimateModel}, diagnostic::Function, domain::Array{Float64,1})
    
    diag_arr = zeros(length(ensemble), length(domain))
    
    count_model = 1
    for (name, model) in ensemble
        diag_arr[count_model,:] = diagnostic(model)
        count_model+=1
    end
    
    myquantile(arr) = quantile(arr, [0.1, 0.5, 0.9])
    
    stats = hcat(map(myquantile, (diag_arr[:,i] for i=1:size(diag_arr,2)))...)
    
    return stats[1,:], stats[2,:], stats[3,:]
end

function ensemble_state_statistics(ensemble::Dict{String, ClimateModel}, symbols::Array{Symbol,1}, domain::Array{Float64,1})
    
    var_arr = zeros(length(ensemble), length(domain))
    
    count_model = 1
    for (name, model) in ensemble
        field = model
        for symbol in symbols
            field = getfield(field, symbol)
        end
        var_arr[count_model,:] = field
        count_model+=1
    end
    
    myquantile(arr) = quantile(arr, [0.1, 0.5, 0.9])
    
    stats = hcat(map(myquantile, (var_arr[:,i] for i=1:size(var_arr,2)))...)
    
    return stats[1,:], stats[2,:], stats[3,:]
end
