rcParams = PyPlot.PyDict(PyPlot.matplotlib."rcParams")
rcParams["lines.linewidth"] = 3 # Change linewidth

function plot_state(model::ClimateModel)
    figure(figsize=(14,8))
    
    subplot(2,3,1)
    title("emissions scenarios")
    plot(model.domain, model.economics.baseline_emissions, label="no-policy baseline")
    plot(model.domain, effective_emissions(model), label="controlled")
    if model.present_year != model.domain[1]
        plot(
            [model.present_year, model.present_year],
            [-maximum(model.economics.baseline_emissions) * 1.1, maximum(model.economics.baseline_emissions) * 1.1],
            "r--"
        )
    end
    plot(model.domain, zeros(size(model.domain)), "k--", alpha=0.5)
    ylabel(L"CO₂ emissions $q$ (ppm / yr)")
    xlim(model.domain[1],model.domain[end])
    ylim(-maximum(model.economics.baseline_emissions) * 1.1, maximum(model.economics.baseline_emissions) * 1.1)
    xlabel("year")
    grid(true)
    legend()
    annotate(s="a)",xy=(0,1.02),xycoords="axes fraction",fontsize=12)

    subplot(2,3,2)
    title("concentrations scenarios")
    plot(model.domain, CO₂_baseline(model), label=L"$c_{0}(t)$ (no-policy baseline)")
    plot(model.domain, CO₂(model), label=L"$c_{\phi,\varphi}(t)$ (controlled)")
    if model.present_year != model.domain[1]
        plot([model.present_year, model.present_year], [0., maximum(CO₂_baseline(model))*1.05], "r--")
    end
    legend()
    ylabel(L"CO₂ concentration $c$ (ppm)")
    xlabel("year")
    xlim(model.domain[1],model.domain[end])
    ylim([0., maximum(CO₂_baseline(model))*1.05])
    grid(true)
    annotate(s="b)",xy=(0,1.02),xycoords="axes fraction",fontsize=12)
    
    subplot(2,3,3)
    title("optimized control deployments")
    plot(model.domain, model.controls.remove, label=L"$\phi$ (negative emissions)")
    plot(model.domain, model.controls.reduce, label=L"$\varphi$ (emissions reductions)")
    plot(model.domain, model.controls.adapt, label=L"$\chi$ (adaptation)")
    plot(model.domain, model.controls.geoeng, label=L"$\lambda$ (geoengineering)")
    if model.present_year != model.domain[1]
        plot([model.present_year, model.present_year], [0., 1.], "r--")
    end
    ylabel(L"fraction of control technology deployed $\alpha$")
    xlabel("year")
    xlim(model.domain[1],model.domain[end])
    ylim([0,1])
    grid(true)
    legend()
    annotate(s="c)",xy=(0,1.02),xycoords="axes fraction",fontsize=12)
    
    subplot(2,3,4)
    title("costs of deploying climate controls")
    plot(model.domain, f_med(model.controls.remove) * model.economics.remove_cost, label=L"$C_{\phi} f(\phi)$ (negative emissions)")
    plot(model.domain, f_med(model.controls.reduce) * model.economics.reduce_cost, label=L"$C_{\varphi} f(\varphi)$ (emissions reductions)")
    plot(model.domain, f_med(model.controls.adapt) * model.economics.adapt_cost, label=L"$C_{\chi} f(\chi)$ (adaptation)")
    plot(model.domain, f_med(model.controls.geoeng) * model.economics.geoeng_cost, label=L"$C_{\lambda} f(\lambda)$ (geoengineering)")
    ylabel(L"cost of climate controls (10$^{12}$ \$ / year)")
    xlabel("year")
    xlim(model.domain[1],model.domain[end])
    grid(true)
    legend()
    annotate(s="d)",xy=(0,1.02),xycoords="axes fraction",fontsize=12)

    subplot(2,3,5)
    title("temperature change since 1850")
    plot(model.domain,δT_baseline(model), label=L"$\delta T_{0}$ (baseline)")
    plot(model.domain,δT(model), label=L"$\delta T_{\varphi,\phi, \lambda}$ (controlled)")
    plot(model.domain,δT_no_geoeng(model), label=L"$\delta T_{\varphi,\phi}$ (controlled without geoengineering)")
    plot(model.domain,2.0.*ones(size(model.domain)),"k--", alpha=0.5)
    plot(model.domain,1.5.*ones(size(model.domain)),"k--", alpha=0.5)
    if model.present_year != model.domain[1]
        plot([model.present_year, model.present_year], [0., maximum(δT_baseline(model)) * 1.05], "r--")
    end
    ylabel(L"warming $δT$ ($^{\circ}$C)")
    xlabel("year")
    xlim(model.domain[1],model.domain[end])
    ylim([0., maximum(δT_baseline(model)) * 1.05])
    legend()
    grid(true)
    annotate(s="e)",xy=(0,1.02),xycoords="axes fraction",fontsize=12)

    subplot(2,3,6)
    plot(model.domain, damage_cost_baseline(model) .* discounting(model), label="uncontrolled damages")
    plot(model.domain, net_cost(model) .* discounting(model), label="net cost (controlled damages + controls)")
    plot(model.domain, damage_cost(model) .* discounting(model), label="controlled damages")
    plot(model.domain, control_cost(model) .* discounting(model), label="cost of controls")
    if model.present_year != model.domain[1]
        plot(
            [model.present_year, model.present_year],
            [0., maximum(damage_cost_baseline(model) .* discounting(model)) * 1.25],
            "r--"
        )
    end
    ylabel(L"discounted costs (10$^{12}$ \$ / year)")
    xlabel("year")
    xlim(model.domain[1],model.domain[end])
    ylim([0., maximum(damage_cost_baseline(model) .* discounting(model)) * 1.25])
    legend()
    grid(true)
    annotate(s="f)",xy=(0,1.02),xycoords="axes fraction",fontsize=12)

    tight_layout()

end

function plot_ensemble_diagnostic(ensemble::Dict{String, ClimateModel}, symbols::Array{Symbol,1}, domain::Array{Float64,1}, color = "C0", label = nothing)
    first, median, ninth = ensemble_state_statistics(ensemble, symbols, domain)
    fill_between(domain, first, ninth, facecolor=color, alpha=0.4)
    plot(domain, median, "-", color=color, alpha=1.0, label=label)
end

function plot_ensemble_statistics(ensemble::Dict{String, ClimateModel}, diagnostic::Function, domain::Array{Float64,1}, color::String, label::String)
    first, median, ninth = ensemble_diagnostic_statistics(ensemble, diagnostic, domain)
    fill_between(domain, first, ninth, facecolor=color, alpha=0.3)
    plot(domain, median, "-", color=color, alpha=1.0, label=label)
end

function plot_ensemble(ensemble::Dict{String, ClimateModel})
    (model, _) = iterate(values(ensemble))
    domain = model.domain
    
    figure(figsize=(14,8))
    
    subplot(2,3,1)
    title("emissions scenarios")
    plot(model.domain, model.economics.baseline_emissions, label="no-policy baseline")
    plot(model.domain, effective_emissions(model), label="controlled")
    if model.present_year != model.domain[1]
        plot(
            [model.present_year, model.present_year],
            [-maximum(model.economics.baseline_emissions) * 1.1, maximum(model.economics.baseline_emissions) * 1.1],
            "r--"
        )
    end
    plot(model.domain, zeros(size(model.domain)), "k--", alpha=0.5)
    ylabel(L"CO₂ emissions $q$ (ppm / yr)")
    xlim(model.domain[1],model.domain[end])
    ylim(-maximum(model.economics.baseline_emissions) * 1.1, maximum(model.economics.baseline_emissions) * 1.1)
    xlabel("year")
    grid(true)
    legend()
    annotate(s="a)",xy=(0,1.02),xycoords="axes fraction",fontsize=12)
    
    subplot(2,3,2)
    title("concentrations scenarios")
    plot(model.domain, CO₂_baseline(model), label=L"$c_{0}(t)$ (no-policy baseline)")
    plot(model.domain, CO₂(model), label=L"$c_{\phi,\varphi}(t)$ (controlled)")
    if model.present_year != model.domain[1]
        plot([model.present_year, model.present_year], [0., maximum(CO₂_baseline(model))*1.05], "r--")
    end
    legend()
    ylabel(L"CO₂ concentration $c$ (ppm)")
    xlabel("year")
    xlim(model.domain[1],model.domain[end])
    ylim([0., maximum(CO₂_baseline(model))*1.05])
    grid(true)
    annotate(s="b)",xy=(0,1.02),xycoords="axes fraction",fontsize=12)
    
    subplot(2,3,3)
    title("optimized control deployments")
    plot(model.domain, model.controls.remove, label=L"$\phi$ (negative emissions)")
    plot(model.domain, model.controls.reduce, label=L"$\varphi$ (emissions reductions)")
    plot(model.domain, model.controls.adapt, label=L"$\chi$ (adaptation)")
    plot(model.domain, model.controls.geoeng, label=L"$\lambda$ (geoengineering)")
    if model.present_year != model.domain[1]
        plot([model.present_year, model.present_year], [0., 1.], "r--")
    end
    ylabel(L"fraction of control technology deployed $\alpha$")
    xlabel("year")
    xlim(model.domain[1],model.domain[end])
    ylim([0,1])
    grid(true)
    legend()
    annotate(s="c)",xy=(0,1.02),xycoords="axes fraction",fontsize=12)
    
    subplot(2,3,4)
    title("costs of deploying climate controls")
    plot(model.domain, f_med(model.controls.remove) * model.economics.remove_cost, label=L"$C_{\phi} f(\phi)$ (negative emissions)")
    plot(model.domain, f_med(model.controls.reduce) * model.economics.reduce_cost, label=L"$C_{\varphi} f(\varphi)$ (emissions reductions)")
    plot(model.domain, f_med(model.controls.adapt) * model.economics.adapt_cost, label=L"$C_{\chi} f(\chi)$ (adaptation)")
    plot(model.domain, f_med(model.controls.geoeng) * model.economics.geoeng_cost, label=L"$C_{\lambda} f(\lambda)$ (geoengineering)")
    ylabel(L"cost of climate controls (10$^{12}$ \$ / year)")
    xlabel("year")
    xlim(model.domain[1],model.domain[end])
    grid(true)
    legend()
    annotate(s="d)",xy=(0,1.02),xycoords="axes fraction",fontsize=12)

    subplot(2,3,5)
    title("temperature change since 1850")
    plot_ensemble_statistics(
        ensemble, δT, domain,
        "C0", L"$\delta T_{\varphi,\phi, \lambda}$ (controlled)"
    )
    plot_ensemble_statistics(
        ensemble, δT_no_geoeng, domain,
        "C1", L"$\delta T_{\varphi,\phi}$ (controlled without geoengineering)"
    )
    plot_ensemble_statistics(
        ensemble, δT_baseline, domain,
        "C2", L"$\delta T_{0}$ (baseline)"
    )
    plot(domain, 2.0.*ones(size(domain)), "k--", label="Paris Goal", alpha=0.5)
    ylabel(L"warming $δT$ ($^{\circ}$C)")
    xlabel("year")
    xlim([domain[1], domain[end]])
    legend(loc="upper left")
    annotate(s="c)",xy=(0,1.02),xycoords="axes fraction",fontsize=12)

    subplot(2,3,6)
    title("discounted costs and benefits")
    plot_ensemble_statistics(
        ensemble, discounted_damage_cost_baseline, domain,
        "C0", "uncontrolled damages"
    )
    plot_ensemble_statistics(
        ensemble, discounted_net_cost, domain,
        "C1", "net cost (controlled damages + controls)"
    )
    plot_ensemble_statistics(
        ensemble, discounted_damage_cost, domain,
        "C2", "controlled damages"
    )
    plot_ensemble_statistics(
        ensemble, discounted_control_cost, domain,
        "C3", "cost of controls"
    )
    ylabel(L"discounted costs (10$^{12}$ \$ / year)")
    xlabel("year")
    xlim([domain[1], domain[end]])
    legend()
    annotate(s="d)",xy=(0,1.02),xycoords="axes fraction",fontsize=12)
    
    tight_layout()
end

