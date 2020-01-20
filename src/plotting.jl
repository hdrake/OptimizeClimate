function plot_state(model::ClimateModel)
    figure(figsize=(14,8))
    
    subplot(2,3,1)
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
    title("emissions scenarios")
    grid(true)
    legend()
    annotate(s="a)",xy=(0,1.02),xycoords="axes fraction",fontsize=12)

    subplot(2,3,2)
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
    title("concentrations scenarios")
    annotate(s="b)",xy=(0,1.02),xycoords="axes fraction",fontsize=12)
    
    subplot(2,3,3)
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
    title("optimized control deployments")
    grid(true)
    legend()
    annotate(s="c)",xy=(0,1.02),xycoords="axes fraction",fontsize=12)
    
    subplot(2,3,4)
    plot(model.domain, f_med(model.controls.remove) * model.economics.remove_cost, label=L"$C_{\phi} f(\phi)$ (negative emissions)")
    plot(model.domain, f_med(model.controls.reduce) * model.economics.reduce_cost, label=L"$C_{\varphi} f(\varphi)$ (emissions reductions)")
    plot(model.domain, f_med(model.controls.adapt) * model.economics.adapt_cost, label=L"$C_{\chi} f(\chi)$ (adaptation)")
    plot(model.domain, f_med(model.controls.geoeng) * model.economics.geoeng_cost, label=L"$C_{\lambda} f(\lambda)$ (geoengineering)")
    ylabel(L"cost of climate controls (10$^{12}$ \$ / year)")
    xlabel("year")
    xlim(model.domain[1],model.domain[end])
    grid(true)
    title("costs of deploying climate controls")
    legend()
    annotate(s="d)",xy=(0,1.02),xycoords="axes fraction",fontsize=12)

    subplot(2,3,5)
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
    title("warming since 1850")
    annotate(s="e)",xy=(0,1.02),xycoords="axes fraction",fontsize=12)

    subplot(2,3,6)
    plot(model.domain, damage_cost_baseline(model) .* discounting(model), label="uncontrolled damages")
    plot(model.domain, damage_cost(model) .* discounting(model), label="controlled damages")
    plot(model.domain, control_cost(model) .* discounting(model), label="cost of controls")
    plot(model.domain, net_cost(model) .* discounting(model), label="net cost (controlled damages + cost of controls)")
    if model.present_year != model.domain[1]
        plot(
            [model.present_year, model.present_year],
            [0., maximum(damage_cost_baseline(model) .* discounting(model)) * 1.25],
            "r--"
        )
    end
    ylabel(L"discounted costs (10$^{12}$ \$)")
    xlabel("year")
    xlim(model.domain[1],model.domain[end])
    ylim([0., maximum(damage_cost_baseline(model) .* discounting(model)) * 1.25])
    legend()
    grid(true)
    annotate(s="f)",xy=(0,1.02),xycoords="axes fraction",fontsize=12)

    tight_layout()

end

function plot_ensemble_state(ensemble::Dict{String, ClimateModel})
    figure(figsize=(10,8))
    
    first = true
    for (name, model) in ensemble
        subplot(2,2,1)
        plot(model.domain, model.controls.remove, "C0-", alpha=0.05, label=L"$\phi$ (negative emissions)")
        plot(model.domain, model.controls.reduce, "C1-", alpha=0.05, label=L"$\varphi$ (emissions reductions)")
        plot(model.domain, model.controls.adapt, "C2-", alpha=0.05, label=L"$\chi$ (adaptation)")
        plot(model.domain, model.controls.geoeng, "C3-", alpha=0.05, label=L"$\lambda$ (geoengineering)")
        
        if first
            ylabel(L"fraction of control technology deployed $\alpha$")
            xlabel("year")
            xlim([model.domain[1], model.domain[end]])
            title("optimized control deployments")
            legend()
            annotate(s="a)",xy=(0,1.02),xycoords="axes fraction",fontsize=12)
        end


        subplot(2,2,2)
        plot(model.domain,CO₂(model), "C0-", alpha=0.05, label=L"$c_{\phi,\varphi}(t)$")
        plot(model.domain,CO₂_baseline(model), "C1-", alpha=0.05, label=L"$c_{0}(t)$")
        if first
            legend()
            ylabel(L"CO₂ concentration $c$ (ppm)")
            xlabel("year")
            xlim([model.domain[1], model.domain[end]])
            title("concentrations scenarios")
            annotate(s="b)",xy=(0,1.02),xycoords="axes fraction",fontsize=12)
        end
            
        subplot(2,2,3)
        plot(model.domain,δT(model), "C0-", alpha=0.05, label=L"$\delta T_{\varphi,\phi, \lambda}$ (controlled)")
        plot(model.domain,δT_no_geoeng(model), "C1-", alpha=0.05, label=L"$\delta T_{\varphi,\phi}$ (controlled without geoengineering)")
        plot(model.domain,δT_baseline(model), "C2-", alpha=0.05, label=L"$\delta T_{0}$ (baseline)")
        
        if first
            plot(model.domain,2.0.*ones(size(model.domain)), "k--", label="Paris Goal", alpha=0.5)
            ylabel(L"warming $δT$ ($^{\circ}$C)")
            xlabel("year")
            xlim([model.domain[1], model.domain[end]])
            legend(loc="upper left")
            title("warming since 1850")
            annotate(s="c)",xy=(0,1.02),xycoords="axes fraction",fontsize=12)
        end

        subplot(2,2,4)
        discount = discounting(model)
        plot(model.domain, net_cost(model) .* discount, "C0-", alpha=0.05, label="total (controlled) cost")
        plot(model.domain, control_cost(model) .* discount, "C1-", alpha=0.05, label="cost of controls")
        plot(model.domain, damage_cost(model) .* discount, "C2-", alpha=0.05, label="controlled damages")
        plot(model.domain, damage_cost_baseline(model) .* discount, "C3-", alpha=0.05, label="uncontrolled damages")
        if first
            ylabel(L"discounted costs (10$^{12}$ \$)")
            xlabel("year")
            xlim([model.domain[1], model.domain[end]])
            legend()
            annotate(s="d)",xy=(0,1.02),xycoords="axes fraction",fontsize=12)
            first=false
        end
    end
    
    tight_layout()
    savefig("../figures/ensemble_state.png", bbox_inches="tight", dpi=200)
end

function plot_ensemble_stats(ensemble::Dict{String, ClimateModel}, domain::Array{Float64,1})
    figure(figsize=(10,8))

    subplot(2,2,1)
    
    first, median, ninth = ensemble_state_statistics(ensemble, [:controls, :remove], domain)
    fill_between(domain, first, ninth, facecolor="C0", alpha=0.4)
    plot(domain, median, "C0-", alpha=1.0, label=L"$\phi$ (negative emissions)")
    
    first, median, ninth = ensemble_state_statistics(ensemble, [:controls, :reduce], domain)
    fill_between(domain, first, ninth, facecolor="C1", alpha=0.4)
    plot(domain, median, "C1-", alpha=1.0, label=L"$\varphi$ (emissions reductions)")
    
    first, median, ninth = ensemble_state_statistics(ensemble, [:controls, :adapt], domain)
    fill_between(domain, first, ninth, facecolor="C2", alpha=0.4)
    plot(domain, median, "C2-", alpha=1.0, label=L"$\chi$ (adaptation)")
    
    first, median, ninth = ensemble_state_statistics(ensemble, [:controls, :geoeng], domain)
    fill_between(domain, first, ninth, facecolor="C3", alpha=0.4)
    plot(domain, median, "C3-", alpha=1.0, label=L"$\lambda$ (geoengineering)")
    
    ylabel(L"fraction of control technology deployed $\alpha$")
    xlabel("year")
    xlim([domain[1], domain[end]])
    title("optimized control deployments")
    legend()
    annotate(s="a)",xy=(0,1.02),xycoords="axes fraction",fontsize=12)

    subplot(2,2,2)
    
    first, median, ninth = ensemble_diagnostic_statistics(ensemble, CO₂, domain)
    fill_between(domain, first, ninth, facecolor="C0", alpha=0.3)
    plot(domain, median, "C0-", label=L"$c_{\phi,\varphi}(t)$")
    
    first, median, ninth = ensemble_diagnostic_statistics(ensemble, CO₂_baseline, domain)
    fill_between(domain, first, ninth, facecolor="C1", alpha=0.3)
    plot(domain, median, "C1-", label=L"$c_{0}(t)$")
    
    legend()
    ylabel(L"CO₂ concentration $c$ (ppm)")
    xlabel("year")
    xlim([domain[1], domain[end]])
    title("concentrations scenarios")
    annotate(s="b)",xy=(0,1.02),xycoords="axes fraction",fontsize=12)

    subplot(2,2,3)
    
    first, median, ninth = ensemble_diagnostic_statistics(ensemble, δT, domain)
    fill_between(domain, first, ninth, facecolor="C0", alpha=0.3)
    plot(domain, median, "C0-", label=L"$\delta T_{\varphi,\phi, \lambda}$ (controlled)")
    
    first, median, ninth = ensemble_diagnostic_statistics(ensemble, δT_no_geoeng, domain)
    fill_between(domain, first, ninth, facecolor="C1", alpha=0.3)
    plot(domain, median, "C1-", label=L"$\delta T_{\varphi,\phi}$ (controlled without geoengineering)")
    
    first, median, ninth = ensemble_diagnostic_statistics(ensemble, δT_baseline, domain)
    fill_between(domain, first, ninth, facecolor="C2", alpha=0.3)
    plot(domain, median, "C2-", label=L"$\delta T_{0}$ (baseline)")

    plot(domain, 2.0.*ones(size(domain)), "k--", label="Paris Goal", alpha=0.5)
    ylabel(L"warming $δT$ ($^{\circ}$C)")
    xlabel("year")
    xlim([domain[1], domain[end]])
    legend(loc="upper left")
    title("warming since 1850")
    annotate(s="c)",xy=(0,1.02),xycoords="axes fraction",fontsize=12)

    subplot(2,2,4)

    first, median, ninth = ensemble_diagnostic_statistics(ensemble, discounted_net_cost, domain)
    fill_between(domain, first, ninth, facecolor="C0", alpha=0.3)
    plot(domain, median, "C0-", label="total (controlled) cost")
    
    first, median, ninth = ensemble_diagnostic_statistics(ensemble, discounted_control_cost, domain)
    fill_between(domain, first, ninth, facecolor="C1", alpha=0.3)
    plot(domain, median, "C1-", label="cost of controls")
    
    first, median, ninth = ensemble_diagnostic_statistics(ensemble, discounted_damage_cost, domain)
    fill_between(domain, first, ninth, facecolor="C2", alpha=0.3)
    plot(domain, median, "C2-", label="controlled damages")
    
    first, median, ninth = ensemble_diagnostic_statistics(ensemble, discounted_damage_cost_baseline, domain)
    fill_between(domain, first, ninth, facecolor="C3", alpha=0.3)
    plot(domain, median, "C3-", label="uncontrolled damages")
    
    ylabel(L"discounted costs (10$^{12}$ \$)")
    xlabel("year")
    xlim([domain[1], domain[end]])
    legend()
    annotate(s="d)",xy=(0,1.02),xycoords="axes fraction",fontsize=12)
    
    tight_layout()
    savefig("../figures/ensemble_stats.png", bbox_inches="tight", dpi=200)
end

