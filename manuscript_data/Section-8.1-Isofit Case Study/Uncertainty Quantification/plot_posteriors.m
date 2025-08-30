clear;
marg_pos = readmatrix("CALF_20_marg_pos.csv");
marg_pos_bin_edges = readmatrix("CALF_20_marg_pos_bin_edge.csv");
bounds = readmatrix("CALF_20_params_bounds.csv");
lb = bounds(1, :);
ub = bounds(2, :);
np = size(marg_pos, 1);
param_names = ["q_{sat,1}", "b_1", "q_{sat,2}", "b_2"];
param_units = ["mol/kg", "Pa^{-1}", "mol/kg", "Pa^{-1}", "-"];
leg_size = 12;
label_size = 24;
idx_scale = [2, 4];
% idx_scale = [];
scale_values = [10^4, 10^6];

% theta_nominal = [2.756811, 8.026405e-04, 3.191820, 4.656369e-06];
theta_nominal = [2.81704, 7.65077e-04, 3.20816, 4.08352e-06];

hist_values = marg_pos;
dist_params_norm = zeros(size(marg_pos, 1), 2);
dist_params_skew_norm = zeros(size(marg_pos, 1), 4);

% for n=1:np-1
%     if ismember(n, idx_scale)
%         idx = find(n==idx_scale);
%         scale_value = scale_values(idx);
%         dist_params_norm(n, :) = fit_dist(marg_pos_bin_edges(n, :)*scale_value,...
%                                             marg_pos(n, :)./scale_value, 2);
%         dist_params_skew_norm(n, :) = fit_dist(marg_pos_bin_edges(n, :)*scale_value,...
%                                             marg_pos(n, :)./scale_value, 3);
%     else
%         if n==2
%         s = 2;
%         end
%         % dist_params_norm(n, :) = fit_dist(marg_pos_bin_edges(n, :),...
%         %                                    marg_pos(n, :), 2, theta_nominal(n));
%         % 
%         % dist_params_skew_norm(n, :) = fit_dist(marg_pos_bin_edges(n, :),...
%                                             % marg_pos(n, :), 3, theta_nominal(n));
%     end
% end


for p=1:np-1
    fig1 = figure(p);
    fig1.Position = [100, 50, 500, 375];
    % subplot(3, 2, p);
    hold on
    % Plot marginal posterior histogram
    if ismember(p, idx_scale)
        idx = find(p==idx_scale);
        scale_value = scale_values(idx);
        h = histogram('BinEdges', marg_pos_bin_edges(p, :).*scale_value,...
              'BinCounts', marg_pos(p, :), 'Normalization', 'pdf');
        
        AIM_fit_line = xline(theta_nominal(p).*scale_value);
        
        x = get_centroid(marg_pos_bin_edges(p, :).*scale_value);
        dist_params_norm(p, :) = fit_dist(marg_pos_bin_edges(p, :).*scale_value, ...
                                            h.Values, 2, theta_nominal(p));
        % dist_params_skew_norm(p, :) = fit_dist(marg_pos_bin_edges(p, :).*scale_value,...
                                            % h.Values, 3, theta_nominal(p));
        
        fitted_dist_norm = fun_norm(dist_params_norm(p, :), x);
        % fitted_dist_skewed = skew_normal_pdf(dist_params_skew_norm(p, :), x);
        % 
        % plot(x, fitted_dist_norm, 'Color', 'red', 'LineWidth', 2.0);
        % plot(x, fitted_dist_skewed, 'Color', 'blue', 'LineWidth', 2.0);

        xlabel(gca, sprintf("%s %c 10^{%d} (%s)", param_names(p), char(215),...
                int16(log10(scale_value)), param_units(p)), "FontSize",label_size);
        xlim(gca, [lb(p), ub(p)].*scale_value);

        % Get axis limits
        ax = gca;
        xLim = ax.XLim;
        yLim = ax.YLim;
        mu = dist_params_norm(p, 1);
        sigma = dist_params_norm(p, 2);
        % Annotate the figure with μ (red) and σ (green), including values in respective colors
        annotationText = sprintf('\\bf\\color{red}\\mu = %.2e\n\\bf\\color{blue}\\sigma = %.2e', mu/scale_value, sigma/scale_value);
        text(xLim(2) - 0.1 * range(xLim), yLim(2) - 0.1 * range(yLim), ...
            annotationText, 'FontSize', label_size-4,  'BackgroundColor', 'w', ...
            'HorizontalAlignment', 'right', 'VerticalAlignment', 'top', ...
            'Interpreter', 'tex');  % Enable TeX interpreter for color formatting
    else
        h = histogram('BinEdges', marg_pos_bin_edges(p, :),...
              'BinCounts', marg_pos(p, :), 'Normalization', 'pdf');
        AIM_fit_line = xline(theta_nominal(p));
          
        x = get_centroid(marg_pos_bin_edges(p, :));
        dist_params_norm(p, :) = fit_dist(marg_pos_bin_edges(p, :), h.Values, 2, theta_nominal(p));
        % dist_params_skew_norm(p, :) = fit_dist(marg_pos_bin_edges(p, :), h.Values, 3, theta_nominal(p));

        % fitted_dist_norm = fun_norm(dist_params_norm(p, :), x);
        % fitted_dist_skewed = skew_normal_pdf(dist_params_skew_norm(p, :), x);
        % plot(x, fitted_dist_norm, 'Color', 'red', 'LineWidth', 2.0);
        % plot(x, fitted_dist_skewed, 'Color', 'blue', 'LineWidth', 2.0);
        
        xlabel(gca, sprintf("%s (%s)", param_names(p), param_units(p)),...
            "FontSize",label_size);
        xlim(gca, [lb(p), ub(p)]);

        % Get axis limits
        ax = gca;
        xLim = ax.XLim;
        yLim = ax.YLim;
        mu = dist_params_norm(p, 1);
        sigma = dist_params_norm(p, 2);

        % Annotate the figure with μ (red) and σ (green), including values in respective colors
        annotationText = sprintf('\\bf\\color{red}\\mu = %.2f\n\\bf\\color{blue}\\sigma = %.2f', mu, sigma);
        text(xLim(2) - 0.1 * range(xLim), yLim(2) - 0.1 * range(yLim), ...
            annotationText, 'FontSize', label_size-4, 'BackgroundColor', 'w', ...
            'HorizontalAlignment', 'right', 'VerticalAlignment', 'top', ...
            'Interpreter', 'tex');  % Enable TeX interpreter for color formatting
    end
    % h.FaceColor = "#0072BD";
    h.FaceColor = "#D95319";
    h.FaceAlpha = 0.7;
    AIM_fit_line.Layer = "top";
    % AIM_fit_line.Color = "#0072BD";
    AIM_fit_line.Color = "black";
    AIM_fit_line.LineStyle = "--";
    AIM_fit_line.LineWidth = 3.0;
    
    ax = gca;
    % xlabel(ax, sprintf("%s, (%s)", param_names(p), param_units(p)), "FontSize",label_size);
    ylabel(ax, sprintf("P(%s|D)", param_names(p)), "FontSize",label_size);
    set(ax, 'YTick', []);
    xtickformat('%.1f')
    hold off
    % grid on
    ax.FontName = 'Deja Vu Sans';
    ax.FontSize = label_size;
    ax.GridAlpha = 0.3;
    box on;
    saveas(fig1, sprintf('CALF-20-CO2-%s.svg', param_names(p)))
end

% close all;

%% Miscellenous Function
function fitted_params = fit_dist(bin_edges, bin_count, num_params, nominal_val)

    x = get_centroid(bin_edges);
    options = optimoptions('lsqcurvefit', 'FunctionTolerance', 1e-12, 'MaxIterations',1e05,...
                        'Display', 'off', 'StepTolerance', 1e-10);
    if num_params == 2
        % x0 = [x(round(end/2)), 0.1*x(round(end/2))];
        x0 = [nominal_val, 1e-03];
        lb = [min(x), 0];       % Lower bound (mu, sigma)
        ub = [max(x), 10];      % Upper bound (mu, sigma)
        % fit_function = @(params) residual_norm(params, x, bin_count);

        problem = createOptimProblem('lsqcurvefit', 'objective', @fun_norm,...
                    'xdata', x, 'ydata', bin_count, 'lb', lb, 'ub', ub, 'x0', x0,...
                    'options', options);
    else
        x0 = [nominal_val, 0.1*x(round(end/2)), 0];
        lb = [min(x), 0, -50];          % Lower bound (location, scale, shape)
        ub = [max(x), max(x), 50];      % Upper bound (location, scale, shape)
        % fit_function = @(params) residual_sk_norm(params, x, bin_count);
        problem = createOptimProblem('lsqcurvefit', 'objective', @skew_normal_pdf,...
                    'xdata', x, 'ydata', bin_count, 'lb', lb, 'ub', ub, 'x0', x0,...
                    'options', options);
    end
    
    % Solver: lsqcurvefit
    ms = MultiStart('Display', 'off', 'UseParallel', true);
    rs = RandomStartPointSet('NumStartPoints', 2000);
    sol = run(ms, problem, rs);
    fitted_params = sol;
    % options = optimoptions("ga", "Display", "off", "FunctionTolerance", 1e-14,...
    %                             "MaxTime", 5000, "MaxGenerations", 1000000);
    % [fitted_params, residual, exitflag] = ga(fit_function, num_params, [], [], [], [], lb, ub,...
    %              [], options);
    
    if num_params > 2
        delta = fitted_params(3)/sqrt(1+fitted_params(3)^2);
        mu_calc = fitted_params(1) + fitted_params(2)*delta*sqrt(2/pi());
        fitted_params = [fitted_params, mu_calc];
    end

    % function residual = residual_norm(params, x, y)
    %     f_pred = fun_norm(params, x);
    %     residual = sum(abs(f_pred - y).^2);
    % end
    % 
    % function residual = residual_sk_norm(params, x, y)
    %     f_pred = skew_normal_pdf(params, x);
    %     residual = sum(abs(f_pred - y).^2);
    % end
end

function x = get_centroid(bin_edges)
    % Calculate bin centroid locations
    x = zeros(1, length(bin_edges) - 1); % Bin centroid locations
    for i = 2:length(bin_edges)
	    x(i-1) = mean([bin_edges(i), bin_edges(i-1)]);
    end
end

function f = fun_norm(params, x)
    % N_p = length(x);
    N_p = 1;
    mu = params(1);
    sigma = params(2);
    SSE = (x-mu).^2;
    multiplier = 1/(sigma*sqrt(2*pi))^N_p;
    exp_part = exp(-SSE/2/sigma^2);
    f = multiplier.*exp_part;
end

function [f] = skew_normal_pdf(params, x)
    % Unpack parameter values	
    e = params(1); % Location parameter
    w = params(2); % Scale parameter
    a = params(3); % Shape parameter
    
    % Caluclate vector of skew-normal PDF values (f)
    phi_pdf = (1/sqrt(2*pi))*exp(-((x-e).^2)/2/(w^2));
    phi_cdf = 0.5*(1 + erf(a*(x-e)/w/sqrt(2)));
    f = 2*phi_pdf.*phi_cdf/w;
end

