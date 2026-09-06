function uncertainty = reconstruct_virial_covariance(isotherm_struc, pressure_data, loading_data, temperature_data, num_a_params)
%RECONSTRUCT_VIRIAL_COVARIANCE Recover the Virial-fit parameter covariance.

    uncertainty = empty_result();
    if isempty(isotherm_struc) || ~isfield(isotherm_struc, 'std_error') || ...
            ~isfield(isotherm_struc, 'fitted_params') || ~isfield(isotherm_struc, 'fun')
        return;
    end

    pressure_columns = size(pressure_data, 2);
    if numel(temperature_data) ~= pressure_columns
        return;
    end

    pressure = [];
    loading = [];
    temperature = [];
    for column_idx = 1:pressure_columns
        current_pressure = pressure_data(:, column_idx);
        current_loading = loading_data(:, column_idx);
        valid = isfinite(current_pressure) & isfinite(current_loading) & ...
                current_pressure > 0 & current_loading ~= 0;
        pressure = [pressure; current_pressure(valid)]; %#ok<AGROW>
        loading = [loading; current_loading(valid)]; %#ok<AGROW>
        temperature = [temperature; repmat(temperature_data(column_idx), sum(valid), 1)]; %#ok<AGROW>
    end

    params = isotherm_struc.fitted_params(:).';
    n_params = numel(params);
    dof = numel(loading) - n_params;
    if dof <= 0 || num_a_params < 1 || num_a_params > n_params
        return;
    end

    observed_log_pressure = log(pressure);
    stored_se = isotherm_struc.std_error(:);
    model_fun = @(candidate_params, ~) isotherm_struc.fun( ...
        candidate_params, loading, temperature, num_a_params);

    weight_names = {'Uniform', 'Biased'};
    weight_vectors = {
        ones(size(loading)), ...
        1 ./ (1 + loading) ...
    };

    best_score = Inf;
    best_covariance = NaN(n_params);
    best_rank = 0;
    best_condition = Inf;
    best_weighting = 'Unavailable';

    for candidate_idx = 1:numel(weight_vectors)
        weights = weight_vectors{candidate_idx};
        prediction = model_fun(params, []);
        residual = weights .* (observed_log_pressure - prediction(:));
        residual_jacobian = numerical_residual_jacobian(model_fun, params, weights);

        information = residual_jacobian.' * residual_jacobian;
        covariance = (sum(residual.^2) ./ dof) .* pinv(information);
        covariance = (covariance + covariance.') ./ 2;
        candidate_se = sqrt(max(real(diag(covariance)), 0));

        comparable = isfinite(stored_se) & stored_se > 0 & ...
                     isfinite(candidate_se) & candidate_se > 0;
        if any(comparable)
            score = norm(log(candidate_se(comparable) ./ stored_se(comparable)));
        else
            score = Inf;
        end

        if score < best_score
            best_score = score;
            best_covariance = covariance;
            best_rank = rank(information);
            best_condition = cond(information);
            best_weighting = weight_names{candidate_idx};
        end
    end

    if any(~isfinite(best_covariance), 'all')
        return;
    end

    uncertainty.covariance = best_covariance;
    uncertainty.degrees_of_freedom = dof;
    uncertainty.jacobian_rank = best_rank;
    uncertainty.condition_number = best_condition;
    uncertainty.weighting = best_weighting;
    uncertainty.diagonal_match_score = best_score;
    uncertainty.identifiable = best_rank == n_params;
end

function jacobian = numerical_residual_jacobian(model_fun, params, weights)
    n_params = numel(params);
    base_prediction = model_fun(params, []);
    jacobian = zeros(numel(base_prediction), n_params);
    step_scale = eps(class(params)).^(1/3);

    for parameter_idx = 1:n_params
        step = step_scale .* max(abs(params(parameter_idx)), 1);
        params_plus = params;
        params_minus = params;
        params_plus(parameter_idx) = params_plus(parameter_idx) + step;
        params_minus(parameter_idx) = params_minus(parameter_idx) - step;
        prediction_plus = model_fun(params_plus, []);
        prediction_minus = model_fun(params_minus, []);
        derivative = (prediction_plus(:) - prediction_minus(:)) ./ (2 .* step);
        jacobian(:, parameter_idx) = -weights .* derivative;
    end
end

function result = empty_result()
    result = struct( ...
        'covariance', NaN(0), ...
        'degrees_of_freedom', NaN, ...
        'jacobian_rank', 0, ...
        'condition_number', Inf, ...
        'weighting', 'Unavailable', ...
        'diagonal_match_score', Inf, ...
        'identifiable', false);
end
