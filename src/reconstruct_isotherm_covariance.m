function uncertainty = reconstruct_isotherm_covariance(isotherm_struc, pressure_data, loading_data, P_saturation)
%RECONSTRUCT_ISOTHERM_COVARIANCE Recover the full covariance used by IsoFit.
%   The current AIM fit structure stores parameter standard errors but not
%   off-diagonal covariance. This function reconstructs the local residual
%   Jacobian for both supported weighting schemes and selects the covariance
%   whose diagonal most closely reproduces the stored standard errors.

    if nargin < 4 || isempty(P_saturation)
        P_saturation = 1.0;
    end

    uncertainty = empty_result();
    if isempty(isotherm_struc) || ~isfield(isotherm_struc, 'std_error') || ...
            ~isfield(isotherm_struc, 'fitted_params') || ~isfield(isotherm_struc, 'fun')
        return;
    end

    pressure = pressure_data(:);
    loading = loading_data(:);
    valid = isfinite(pressure) & isfinite(loading);
    pressure = pressure(valid);
    loading = loading(valid);

    params = isotherm_struc.fitted_params(:).';
    n_params = numel(params);
    n_observations = numel(loading);
    dof = n_observations - n_params;
    if dof <= 0 || n_observations == 0
        return;
    end

    if isfield(isotherm_struc, 'p_sat') && isotherm_struc.p_sat
        predictor = pressure ./ P_saturation;
    else
        predictor = pressure;
    end

    model_fun = isotherm_struc.fun;
    stored_se = isotherm_struc.std_error(:);
    weight_names = {'Uniform', 'Biased'};
    weight_vectors = {
        ones(size(loading)), ...
        mean(loading) ./ max(abs(loading), sqrt(eps(class(loading)))) ...
    };

    best_score = Inf;
    best_covariance = NaN(n_params);
    best_rank = 0;
    best_condition = Inf;
    best_weighting = 'Unavailable';

    for candidate_idx = 1:numel(weight_vectors)
        weights = weight_vectors{candidate_idx};
        prediction = model_fun(params, predictor);
        residual = weights .* (loading - prediction(:));
        residual_jacobian = numerical_residual_jacobian( ...
            model_fun, params, predictor, weights);

        information = residual_jacobian.' * residual_jacobian;
        information_rank = rank(information);
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
            best_rank = information_rank;
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

function jacobian = numerical_residual_jacobian(model_fun, params, predictor, weights)
    n_params = numel(params);
    base_prediction = model_fun(params, predictor);
    jacobian = zeros(numel(base_prediction), n_params);
    step_scale = eps(class(params)).^(1/3);

    for parameter_idx = 1:n_params
        step = step_scale .* max(abs(params(parameter_idx)), 1);
        params_plus = params;
        params_minus = params;
        params_plus(parameter_idx) = params_plus(parameter_idx) + step;
        params_minus(parameter_idx) = params_minus(parameter_idx) - step;
        prediction_plus = model_fun(params_plus, predictor);
        prediction_minus = model_fun(params_minus, predictor);
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
