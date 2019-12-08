function plot_setup(model::ClimateModel)

    figure(figsize=(8,4))
    subplot(1,2,1)
    plot(model.domain, baseline_emissions(model.domain))
    ylabel(L"CO₂ emissions $q$ (ppm / yr)")
    xlabel("year")
    title("baseline emissions")
    annotate(s="a)",xy=(0,1.02),xycoords="axes fraction",fontsize=12)

    subplot(1,2,2)
    plot(model.domain, 1. .-discounting(model.economics, model.domain))
    plot(model.domain, ones(size(model.domain)), "r--")
    xlabel("year")
    ylabel("fraction of cost discounted")
    xlim(model.domain[1],model.domain[end])
    ylim(0,1.05)
    annotate(s="c)",xy=(0,1.02),xycoords="axes fraction",fontsize=12)

    tight_layout()

    savefig("figures/model_setup.png", bbox_inches="tight", dpi=100)
end

function plot_state(model::ClimateModel)
    figure(figsize=(10,8))
    subplot(2,2,1)
    plot(model.domain, model.controls.remove, label=L"$\phi$ (negative emissions)")
    plot(model.domain, model.controls.reduce, label=L"$\varphi$ (emissions reductions)")
    plot(model.domain, model.controls.adapt, label=L"$\chi$ (adaptation)")
    plot(model.domain, model.controls.geoeng, label=L"$\lambda$ (geoengineering)")
    ylabel(L"fraction of control technology deployed $\alpha$")
    xlabel("year")
    title("optimized control deployments")
    legend()
    annotate(s="a)",xy=(0,1.02),xycoords="axes fraction",fontsize=12)


    subplot(2,2,2)
    plot(model.domain,CO₂(model), label=L"$c_{\phi,\varphi}(t)$")
    plot(model.domain,CO₂_baseline(model), label=L"$c_{0}(t)$")
    legend()
    ylabel(L"CO₂ concentration $c$ (ppm)")
    xlabel("year")
    title("concentrations scenarios")
    annotate(s="b)",xy=(0,1.02),xycoords="axes fraction",fontsize=12)

    subplot(2,2,3)
    plot(model.domain,δT(model), label=L"$\delta T_{\varphi,\phi, \lambda}$ (controlled)")
    plot(model.domain,δT_no_geoeng(model), label=L"$\delta T_{\varphi,\phi}$ (controlled without geoengineering)")
    plot(model.domain,δT_baseline(model), label=L"$\delta T_{0}$ (baseline)")
    plot(model.domain,2.0.*ones(size(model.domain)),"k--", label="Paris Goal")
    ylabel(L"warming $δT$ ($^{\circ}$C)")
    xlabel("year")
    ylim([0,4.0])
    legend()
    title("warming since 1850")
    annotate(s="c)",xy=(0,1.02),xycoords="axes fraction",fontsize=12)

    subplot(2,2,4)
    plot(model.domain, net_cost(model) .* discounting(model.economics, model.domain), label="total (controlled) cost")
    plot(model.domain, control_cost(model) .* discounting(model.economics, model.domain), label="cost of controls")
    plot(model.domain, damage_cost(model) .* discounting(model.economics, model.domain), label="controlled damages")
    plot(model.domain, damage_cost_baseline(model) .* discounting(model.economics, model.domain), label="uncontrolled damages")
    ylabel(L"discounted costs (10$^{12}$ \$)")
    xlabel("year")
    legend()
    annotate(s="d)",xy=(0,1.02),xycoords="axes fraction",fontsize=12)

    tight_layout()

    savefig("figures/model_state.png", bbox_inches="tight", dpi=100)
end

function plot_ensemble(ensemble::Dict{String, ClimateModel})
    figure(figsize=(10,8))
    
    first = true
    for (name, model) in ensemble
        subplot(2,2,1)
        plot(model.domain, model.controls.remove, "C0-", label=L"$\phi$ (negative emissions)")
        plot(model.domain, model.controls.reduce, "C1-", label=L"$\varphi$ (emissions reductions)")
        plot(model.domain, model.controls.adapt, "C2-", label=L"$\chi$ (adaptation)")
        plot(model.domain, model.controls.geoeng, "C3-", label=L"$\lambda$ (geoengineering)")
        
        if first
            ylabel(L"fraction of control technology deployed $\alpha$")
            xlabel("year")
            title("optimized control deployments")
            legend()
            annotate(s="a)",xy=(0,1.02),xycoords="axes fraction",fontsize=12)
        end


        subplot(2,2,2)
        plot(model.domain,CO₂(model), "C0-", label=L"$c_{\phi,\varphi}(t)$")
        plot(model.domain,CO₂_baseline(model), "C1-", label=L"$c_{0}(t)$")
        if first
            legend()
            ylabel(L"CO₂ concentration $c$ (ppm)")
            xlabel("year")
            title("concentrations scenarios")
            annotate(s="b)",xy=(0,1.02),xycoords="axes fraction",fontsize=12)
        end
            
        subplot(2,2,3)
        plot(model.domain,δT(model), "C0-", label=L"$\delta T_{\varphi,\phi, \lambda}$ (controlled)")
        plot(model.domain,δT_no_geoeng(model), "C1-", label=L"$\delta T_{\varphi,\phi}$ (controlled without geoengineering)")
        plot(model.domain,δT_baseline(model), "C2-", label=L"$\delta T_{0}$ (baseline)")
        
        if first
            plot(model.domain,2.0.*ones(size(model.domain)), "k--", label="Paris Goal")
            ylabel(L"warming $δT$ ($^{\circ}$C)")
            xlabel("year")
            ylim([0,4.0])
            legend()
            title("warming since 1850")
            annotate(s="c)",xy=(0,1.02),xycoords="axes fraction",fontsize=12)
        end

        subplot(2,2,4)
        discount = discounting(model.economics, model.domain)
        plot(model.domain, net_cost(model) .* discount, "C0-", label="total (controlled) cost")
        plot(model.domain, control_cost(model) .* discount, "C1-", label="cost of controls")
        plot(model.domain, damage_cost(model) .* discount, "C2-", label="controlled damages")
        plot(model.domain, damage_cost_baseline(model) .* discount, "C3-", label="uncontrolled damages")
        if first
            ylabel(L"discounted costs (10$^{12}$ \$)")
            xlabel("year")
            legend()
            annotate(s="d)",xy=(0,1.02),xycoords="axes fraction",fontsize=12)
            first=false
        end
    end
    
    tight_layout()
    savefig("figures/ensemble_state.png", bbox_inches="tight", dpi=100)
end