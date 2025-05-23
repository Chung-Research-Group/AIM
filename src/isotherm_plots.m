function [] = isotherm_plots(ax, pressure_data, loading_data, T_flag, T_array, isotherm_struc,...
                             xlogscale_flag, ylogscale_flag, unit_pressure, unit_loading, unit_temperature, P_saturation, vir_flag, vir_num_a)    
    % P_saturation is required for Klotz, Dubinin-Astakhov, and Do-Do model
    if nargin<11
        P_saturation = 1.0;
    end
    %% Function variables    
    line_width = 2.5;
    scatter_marker_size = 90;
    fontsize = 18;
    fontweight = 'bold';
    % fontname = 'AvantGrande';
    fontname = 'Calibri';
    
    legend_array = {};      % Array to hold legend entries
    % color_array = ["blue", "green", "red"];       % Array to hold colors
    %% Plotting isotherm data
    cla (ax, "reset");
    colororder(ax, "gem12");

    hold(ax, 'on');        % For multiple plots
    if ~T_flag
        % Plot isotherm for single temperature
        s_plots = scatter(ax, pressure_data, loading_data, scatter_marker_size, 'filled', 'o', 'MarkerEdgeColor','k');
        s_plots.SizeData = scatter_marker_size;
        legend_array = [legend_array, 'Data'];

    else
        % Check the dimensions of pressure data, loading data are consistent with temperature
        if ~((size(pressure_data, 2) == length(T_array)) ...
             && (size(loading_data, 2) == length(T_array)))
            error("The of pressure/loading data must be specified for every temperature value");
        end
        
        % Plotting isotherms for different temperatures
        for i=1:length(T_array)
            s_plots = scatter(ax, pressure_data(:, i), loading_data(:, i), 'filled', 'o', 'MarkerEdgeColor','k');
            s_plots.SizeData = scatter_marker_size;
            if ~isempty(unit_temperature)
                legend_entry = {sprintf('Data %.1f(%s)', T_array(i), unit_temperature)};
            else
                legend_entry = {sprintf('Data %.1f', T_array(i))};
            end
            legend_array = [legend_array, legend_entry];
        end
    end
    
    %% Plotting isotherm curve using fitted isotherm function 
    if ~(isempty(isotherm_struc))               % If isotherm struc is specified with all the parameters
        isotherm_fun = isotherm_struc.fun;
        isotherm_params = isotherm_struc.fitted_params;
        p_sat_flag = isotherm_struc.p_sat;

        if ~T_flag
            if p_sat_flag
                loading_pred = isotherm_fun(isotherm_params, pressure_data./P_saturation);
            else
                loading_pred = isotherm_fun(isotherm_params, pressure_data);
            end
            q = plot(ax, pressure_data, loading_pred);
            q.LineWidth = line_width;
            % q.Color = char(color_array);
            legend_array = [legend_array, 'Isotherm Fit'];
        else
            if ~vir_flag
                dH = isotherm_struc.dH;
                T_ref = isotherm_struc.T_ref;
                norm_pressure_data = pressure_data .* exp(-1.0 * dH/8.3144 .* (1./T_array - 1/T_ref));
                predicted_loadings = zeros(size(pressure_data, 1), size(T_array, 2));
 
                for j=1:length(T_array)      
                    predicted_loadings(:, j) = isotherm_fun(isotherm_params, norm_pressure_data(:, j)); 
                    q = plot(ax, pressure_data(:, j), predicted_loadings(:, j));
                    q.LineWidth = line_width;  
                end
            else
                pressure_pred = exp(isotherm_fun(isotherm_params, loading_data, T_array, vir_num_a));
                for k=1:length(T_array)
                    q = plot(ax, pressure_pred(:, k), loading_data(:, k));
                    q.LineWidth = line_width;           
                end
            end
            
            for m=1:length(T_array)  
                if ~isempty(unit_temperature)
                    legend_entry = {sprintf('Isotherm Fit %.1f(%s)', T_array(m), unit_temperature)};
                else
                    legend_entry = {sprintf('Isotherm Fit %.1f', T_array(m))};
                end
                legend_array = [legend_array, legend_entry];
            end
        end
    end
    
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

    if ~isempty(unit_loading)
        y_label = sprintf("Gas Uptake (%s)", unit_loading);
    else
        y_label = sprintf("Gas Uptake");
    end

    xlabel(ax, x_label, FontSize = fontsize, FontName=fontname, FontWeight=fontweight);
    ylabel(ax, y_label, FontSize = fontsize, FontName=fontname, FontWeight=fontweight);
    
    legend(ax, legend_array, Location="best");
    
    grid (ax, "on")
    ax.GridAlpha = 0.05;

    ax.FontName = fontname;
    ax.FontWeight = fontweight;
    ax.FontSize = fontsize;
    hold(ax, 'off');

    % Ensuring consistent color scheme for plots
    % Grab handles of scatter plots
    scatter_handles = findobj(ax, 'Type', 'Scatter');
    % if isempty(scatter_handles)
    %     scatter_handles = 0;
    % end

    % Grab handles of line plots
    line_handles = findobj(ax, 'Type', 'Line');
    % if isempty(line_handles)
    %     line_handles = 0;
    % end

    % Get color order from the axes
    colors = ax.ColorOrder;
    color_idx = 1;

    for i=1:length(scatter_handles)
        % Setting the corresponding scatter and line pair to have the same
        % color index        
        if ~isempty(scatter_handles)
            set(scatter_handles(i), 'CData', colors(color_idx, :));
        end

        if ~isempty(line_handles)
            set(line_handles(i), 'Color', colors(color_idx, :));
        end

        % Max colors in the paletter are 12, hence, need to cycle through
        % them.
        if rem(i, length(colors)) == 0
            color_idx = 1;
        else
            color_idx = color_idx+1;
        end
    end
end