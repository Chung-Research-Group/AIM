function [] = mixture_plots(ax, num_comp, name_comp, pressure, mixture_loadings,...
                            xlogscale_flag, ylogscale_flag, plot_type, unit_pressure, unit_loading)    
    %% Function to plot isotherms
    
    %% Function variables
    % colors = ["green", "blue", "red"];
    line_width = 2.5;
    scatter_marker_size = 90;
    line_marker_size = 3.0;
    fontsize = 18;
    fontweight = 'bold';
    % fontname = 'AvantGrande';
    fontname = 'Calibri';
    colororder(ax, "gem12");
    %% Plotting isotherm data
    cla (ax, "reset");
    hold(ax, 'on');        % For multiple plots
    if ~plot_type
        for k=1:num_comp
            q = plot(ax, pressure, mixture_loadings(:, k));
            q.LineWidth = line_width;

            y_label_string = "Gas Uptake";
        end
    else
        adsorbed_MF = mixture_loadings ./ sum(mixture_loadings, 2);
        for k=1:num_comp
            q = plot(ax, pressure, adsorbed_MF(:, k));
            q.LineWidth = line_width;

            y_label_string = "Adsorbed Mole Fractions";
            ylim(ax, [0, 1.00])
        end
    end

    hold(ax, 'off');

    % Checking for log scale requirement
    if xlogscale_flag
        xscale(ax, "log");
    end

    if ylogscale_flag
        yscale(ax, "log");
    end

    ax.Box = "on";
    ax.Color = [1 1 1];

    if ~isempty(unit_pressure)
        x_label = sprintf("Pressure (%s)", unit_pressure);
    else
        x_label = sprintf("Pressure");
    end

    if ~isempty(unit_loading) && ~plot_type
        y_label = sprintf("%s (%s)", y_label_string, unit_loading);
    elseif plot_type
        y_label = sprintf("%s (-)", y_label_string);
    end
    
    xlabel(ax, x_label, FontSize = fontsize, FontName=fontname, FontWeight=fontweight);
    ylabel(ax, y_label, FontSize = fontsize, FontName=fontname, FontWeight=fontweight);
    
    % ax.Title.FontWeight =fontweight;
    % ax.Title.FontSize = fontsize+2;
    % ax.Title.FontName=fontname;
    legend(ax, name_comp, Location="best");
    
    grid (ax, "on")
    ax.GridAlpha = 0.05;

    ax.FontName = fontname;
    ax.FontWeight = fontweight;
    ax.FontSize = fontsize;
end