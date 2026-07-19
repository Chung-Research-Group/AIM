function [] = heat_plot(ax, pressure_data, loading_data, T_flag, T_array, isotherm_struc, ...
                             xlogscale_flag, ylogscale_flag, unit_loading)
%HEAT_PLOT Plot heat of adsorption with a 95% pointwise confidence band.

    line_width = 2.5;
    scatter_marker_size = 60;
    fontsize = 18;
    fontweight = 'bold';
    fontname = 'Calibri';
    gas_constant = 8.314;
    alpha = 0.05;

    if isempty(isotherm_struc)
        return;
    end

    loading_vector = loading_data(:);
    valid_loading = isfinite(loading_vector) & loading_vector ~= 0;
    loading_vector = loading_vector(valid_loading);

    dH = isotherm_struc.dH(:);
    if isscalar(dH)
        dH = repmat(dH, size(loading_vector));
    elseif numel(dH) ~= numel(loading_vector)
        error('Inconsistent size of loading and heat of adsorption data.');
    end

    cla(ax, 'reset');
    colororder(ax, 'gem12');
    color = ax.ColorOrder(1, :);
    hold(ax, 'on');

    ci_lower = NaN(size(dH));
    ci_upper = NaN(size(dH));

    if isfield(isotherm_struc, 'name') && strcmpi(isotherm_struc.name, 'Virial')
        num_a_params = infer_num_a_parameters( ...
            isotherm_struc.fitted_params, loading_vector, dH, gas_constant);
        uncertainty = reconstruct_virial_covariance( ...
            isotherm_struc, pressure_data, loading_data, T_array, num_a_params);
        heat_model = @(params, loading) gas_constant .* ...
            polyval(flip(params(1:num_a_params)), loading);
        heat_ci = delta_method_ci(heat_model, isotherm_struc.fitted_params, ...
            uncertainty.covariance, loading_vector, alpha);
        dH = heat_ci.estimate;
        ci_lower = heat_ci.lower;
        ci_upper = heat_ci.upper;

    elseif isfield(isotherm_struc, 'dH_std_error') && ...
            isfinite(isotherm_struc.dH_std_error)
        standard_error = abs(isotherm_struc.dH_std_error);
        critical_value = -sqrt(2) .* erfcinv(alpha);
        ci_lower = dH - critical_value .* standard_error;
        ci_upper = dH + critical_value .* standard_error;
    end

    [loading_sorted, sort_idx] = sort(loading_vector);
    plotted_heat = -dH(sort_idx) ./ 1e03;
    plotted_lower = -ci_upper(sort_idx) ./ 1e03;
    plotted_upper = -ci_lower(sort_idx) ./ 1e03;

    valid_ci = isfinite(plotted_lower) & isfinite(plotted_upper);
    if nnz(valid_ci) >= 2
        fill(ax, ...
            [loading_sorted(valid_ci); flipud(loading_sorted(valid_ci))], ...
            [plotted_lower(valid_ci); flipud(plotted_upper(valid_ci))], ...
            color, 'FaceAlpha', 0.16, 'EdgeColor', 'none');
    end

    plot(ax, loading_sorted, plotted_heat, 'Color', color, 'LineWidth', line_width);
    scatter(ax, loading_sorted, plotted_heat, scatter_marker_size, color, ...
        'filled', 'o', 'MarkerEdgeColor', 'k');

    if nnz(valid_ci) >= 2
        legend(ax, {'95% pointwise CI', 'Heat of adsorption', 'Calculated values'}, ...
            'Location', 'best');
    else
        legend(ax, {'Heat of adsorption', 'Calculated values'}, 'Location', 'best');
    end

    if xlogscale_flag
        xscale(ax, 'log');
    end
    if ylogscale_flag
        yscale(ax, 'log');
    end

    ax.Box = 'on';
    ax.Color = [1 1 1];
    y_label = strcat(char(8722), '\DeltaH_{ads} (kJ/mol)');
    if ~isempty(unit_loading)
        x_label = sprintf('Gas Uptake (%s)', unit_loading);
    else
        x_label = 'Gas Uptake';
    end

    xlabel(ax, x_label, FontSize=fontsize, FontName=fontname, FontWeight=fontweight);
    ylabel(ax, y_label, 'Interpreter', 'tex', ...
        FontSize=fontsize, FontName=fontname, FontWeight=fontweight);
    grid(ax, 'on');
    ax.GridAlpha = 0.05;
    ytickformat(ax, '%.1f');
    ax.FontName = fontname;
    ax.FontWeight = fontweight;
    ax.FontSize = fontsize;
    hold(ax, 'off');
end

function num_a_params = infer_num_a_parameters(params, loading, dH, gas_constant)
% Select the smallest leading Virial polynomial that reproduces stored dH.
    max_num_a = numel(params);
    scale = max(norm(dH), 1);
    tolerance = 1e-9 .* scale;
    num_a_params = max_num_a;

    for candidate = 1:max_num_a
        prediction = gas_constant .* polyval(flip(params(1:candidate)), loading);
        if norm(prediction(:) - dH(:)) <= tolerance
            num_a_params = candidate;
            return;
        end
    end
end
