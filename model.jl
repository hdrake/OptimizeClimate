struct Controls
    reduce::Array{Float64,1}
    remove::Array{Float64,1}
    geoeng::Array{Float64,1}
    adapt::Array{Float64,1}
end

struct Economics
    β::Float64
    utility_discount_rate::Float64
    
    reduce_cost::Float64
    remove_cost::Float64
    geoeng_cost::Float64
    adapt_cost::Float64
    
    reduce_init::Float64
    remove_init::Float64
    geoeng_init::Float64
    adapt_init::Float64
end

# economic parameters
β = 1.; # damage parameter (10^12 $ / C^2)
utility_discount_rate = 0.014 #2 # utility discount rate (Stern review value)
GWP = 100. # Global World Product (10^12 $/yr)

Economics() = Economics(
    β, utility_discount_rate,
    0.05*GWP, 0.05*GWP, 0.25*GWP, 0.15*GWP,
    0., 0., 0., 0.
)

nondim_linear(t::Array) = (t .- t[1])/(t[end] - t[1]);

function init_linear_controls(t::Array{Float64,1})
    c = Controls(
        nondim_linear(t)/3.,
        nondim_linear(t)/2.,
        nondim_linear(t)/8.,
        nondim_linear(t)/10.
    )
    return c
end

struct ClimateModel
    name::String
    ECS::Float64
    domain::Array{Float64,1}
    controls::Controls
    economics::Economics
    present_year::Float64
    CO₂_init::Float64
    δT_pre::Float64
    
    ϵ::Float64
    
    function ClimateModel(name, ECS, domain, controls, economics, present_year, CO₂_init, δT_pre)
        
        ϵ = ECS/log(2.); # Transient Warming Parameter
        
        return new(
            name, ECS, domain, controls, economics, present_year, CO₂_init, δT_pre,
            ϵ)
    end
end

ClimateModel(name, ECS, domain, controls, economics, present_year) = ClimateModel(
    name, ECS, domain, controls, economics, present_year, 415., 1.1
)