function [] = isotherm_plots(ax, pressure_data, loading_data, T_flag, T_array, isotherm_struc, ...
                             xlogscale_flag, ylogscale_flag, unit_pressure, unit_loading, unit_temperature, ...
                             P_saturation, vir_flag, vir_num_a)
%ISOTHERM_PLOTS Plot adsorption data, fitted curves, and 95% confidence bands.

    if nargin < 12 || isempty(P_saturation)
        P_saturation = 1.0;
    end
    if nargin < 13 || isempty(vir_flag)
        vir_flag = false;
    end
    if nargin < 14
        vir_num_a = [];
    end

    line_width = 2.5;
    scatter_marker_size = 90;
    fontsize = 18;
    fontweight = 'bold';
    fontname = 'Calibri';
    alpha = 0.05;
    legend_array = {};

    cla(ax, 'reset');
    colororder(ax, 'gem12');
    colors = ax.ColorOrder;
    hold(ax, 'on');

    if ~T_flag
        data_color = colors(1, :);
        scatter(ax, pressure_data, loading_data, scatter_marker_size, data_color, ...
            'filled', 'o', 'MarkerEdgeColor', 'k');
        legend_array{end+1} = 'Data';
    else
        if ~((size(pressure_data, 2) == length(T_array)) && ...
             (size(loading_data, 2) == length(T_array)))
            error('The pressure/loading data must be specified for every temperature value.');
        end

        for temperature_idx = 1:length(T_array)
            current_color = colors(mod(temperature_idx - 1, size(colors, 1)) + 1, :);
            scatter(ax, pressure_data(:, temperature_idx), loading_data(:, temperature_idx), ...
                scatter_marker_size, current_color, 'filled', 'o', 'MarkerEdgeColor', 'k');
            legend_array{end+1} = temperature_label( ...
                'Data', T_array(temperature_idx), unit_temperature);
        end
    end

    if ~isempty(isotherm_struc)
        isotherm_fun = isotherm_struc.fun;
        isotherm_params = isotherm_struc.fitted_params;
        p_sat_flag = isotherm_struc.p_sat;

        if ~T_flag
            valid = isfinite(pressure_data) & isfinite(loading_data);
            pressure_fit = pressure_data(valid);
            loading_fit = loading_data(valid);
            [pressure_grid, sort_idx] = sort(pressure_fit(:));
            loading_fit = loading_fit(sort_idx); %#ok<NASGU>

            if p_sat_flag
                model_fun = @(params, pressure) isotherm_fun(params, pressure ./ P_saturation);
            else
                model_fun = @(params, pressure) isotherm_fun(params, pressure);
            end

            uncertainty = reconstruct_isotherm_covariance( ...
                isotherm_struc, pressure_fit, loading_data(valid), P_saturation);
            ci = delta_method_ci(model_fun, isotherm_params, ...
                uncertainty.covariance, pressure_grid, alpha);

            current_color = colors(1, :);
            add_vertical_band(ax, pressure_grid, ci.lower, ci.upper, current_color);
            q = plot(ax, pressure_grid, ci.estimate, 'Color', current_color, ...
                'LineWidth', line_width);
            q.DisplayName = 'Isotherm Fit';
            legend_array{end+1} = 'Isotherm Fit';
            if all(isfinite(ci.lower), 'all')
                legend_array{end+1} = '95% pointwise CI';
                add_ci_legend_proxy(ax, current_color);
            end

        elseif ~vir_flag
            dH = isotherm_struc.dH;
            T_ref = isotherm_struc.T_ref;
            gas_constant = 8.3144;
            norm_pressure_data = pressure_data .* ...
                exp(-dH ./ gas_constant .* (1 ./ T_array - 1 ./ T_ref));

            covariance_fit = reconstruct_isotherm_covariance( ...
                isotherm_struc, norm_pressure_data, loading_data, 1.0);
            parameter_covariance = covariance_fit.covariance;
            dH_standard_error = NaN;
            if isfield(isotherm_struc, 'dH_std_error')
                dH_standard_error = isotherm_struc.dH_std_error;
            end

            for temperature_idx = 1:length(T_array)
                valid = isfinite(pressure_data(:, temperature_idx)) & ...
                        isfinite(loading_data(:, temperature_idx));
                pressure_grid = sort(pressure_data(valid, temperature_idx));
                current_temperature = T_array(temperature_idx);
                current_color = colors(mod(temperature_idx - 1, size(colors, 1)) + 1, :);

                if isfinite(dH_standard_error)
                    combined_params = [isotherm_params(:).', dH];
                    combined_covariance = blkdiag(parameter_covariance, dH_standard_error.^2);
                    model_fun = @(params, pressure) isotherm_fun( ...
                        params(1:end-1), pressure .* exp(-params(end) ./ gas_constant .* ...
                        (1 ./ current_temperature - 1 ./ T_ref)));
                else
                    combined_params = isotherm_params(:).';
                    combined_covariance = parameter_covariance;
                    model_fun = @(params, pressure) isotherm_fun( ...
                        params, pressure .* exp(-dH ./ gas_constant .* ...
                        (1 ./ current_temperature - 1 ./ T_ref)));
                end

                ci = delta_method_ci(model_fun, combined_params, ...
                    combined_covariance, pressure_grid, alpha);
                add_vertical_band(ax, pressure_grid, ci.lower, ci.upper, current_color);
                plot(ax, pressure_grid, ci.estimate, 'Color', current_color, ...
                    'LineWidth', line_width);
                legend_array{end+1} = temperature_label( ...
                    'Isotherm Fit', current_temperature, unit_temperature);
            end

            if any(isfinite(parameter_covariance), 'all')
                add_ci_legend_proxy(ax, colors(1, :));
                legend_array{end+1} = '95% pointwise CI';
            end

        else
            virial_uncertainty = reconstruct_virial_covariance( ...
                isotherm_struc, pressure_data, loading_data, T_array, vir_num_a);

            for temperature_idx = 1:length(T_array)
                valid = isfinite(loading_data(:, temperature_idx)) & ...
                        isfinite(pressure_data(:, temperature_idx)) & ...
                        pressure_data(:, temperature_idx) > 0;
                loading_grid = sort(loading_data(valid, temperature_idx));
                current_temperature = T_array(temperature_idx);
                current_color = colors(mod(temperature_idx - 1, size(colors, 1)) + 1, :);
                model_fun = @(params, loading) isotherm_fun( ...
                    params, loading, repmat(current_temperature, size(loading)), vir_num_a);
                log_pressure_ci = delta_method_ci(model_fun, isotherm_params, ...
                    virial_uncertainty.covariance, loading_grid, alpha);

                pressure_estimate = exp(log_pressure_ci.estimate);
                pressure_lower = exp(log_pressure_ci.lower);
                pressure_upper = exp(log_pressure_ci.upper);
                add_horizontal_band(ax, loading_grid, pressure_lower, pressure_upper, current_color);
                plot(ax, pressure_estimate, loading_grid, 'Color', current_color, ...
                    'LineWidth', line_width);
                legend_array{end+1} = temperature_label( ...
                    'Isotherm Fit', current_temperature, unit_temperature);
            end

            if any(isfinite(virial_uncertainty.covariance), 'all')
                add_ci_legend_proxy(ax, colors(1, :));
                legend_array{end+1} = '95% pointwise CI';
            end
        end
    end

    if xlogscale_flag
        xscale(ax, 'log');
    end
    if ylogscale_flag
        yscale(ax, 'log');
    end

    ax.Box = 'on';
    ax.Color = [1 1 1];
    xlabel(ax, axis_label('Pressure', unit_pressure), ...
        FontSize=fontsize, FontName=fontname, FontWeight=fontweight);
    ylabel(ax, axis_label('Gas Uptake', unit_loading), ...
        FontSize=fontsize, FontName=fontname, FontWeight=fontweight);

    if ~isempty(legend_array)
        legend(ax, legend_array, 'Location', 'best');
    end
    grid(ax, 'on');
    ax.GridAlpha = 0.05;
    ax.FontName = fontname;
    ax.FontWeight = fontweight;
    ax.FontSize = fontsize;
    hold(ax, 'off');
end

function add_vertical_band(ax, x, lower, upper, color)
    valid = isfinite(x) & isfinite(lower) & isfinite(upper);
    x = x(valid);
    lower = lower(valid);
    upper = upper(valid);
    if numel(x) < 2
        return;
    end
    fill(ax, [x; flipud(x)], [lower; flipud(upper)], color, ...
        'FaceAlpha', 0.16, 'EdgeColor', 'none', 'HandleVisibility', 'off');
end

function add_horizontal_band(ax, y, lower_x, upper_x, color)
    valid = isfinite(y) & isfinite(lower_x) & isfinite(upper_x) & lower_x > 0;
    y = y(valid);
    lower_x = lower_x(valid);
    upper_x = upper_x(valid);
    if numel(y) < 2
        return;
    end
    fill(ax, [lower_x; flipud(upper_x)], [y; flipud(y)], color, ...
        'FaceAlpha', 0.16, 'EdgeColor', 'none', 'HandleVisibility', 'off');
end

function add_ci_legend_proxy(ax, color)
    patch(ax, NaN, NaN, color, 'FaceAlpha', 0.16, 'EdgeColor', 'none');
end

function label = temperature_label(prefix, temperature, unit_temperature)
    if ~isempty(unit_temperature)
        label = sprintf('%s %.1f(%s)', prefix, temperature, unit_temperature);
    else
        label = sprintf('%s %.1f', prefix, temperature);
    end
end

function label = axis_label(quantity, unit)
    if ~isempty(unit)
        label = sprintf('%s (%s)', quantity, unit);
    else
        label = quantity;
    end
end
