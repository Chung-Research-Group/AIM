function ci = delta_method_ci(model_fun, fitted_params, covariance, predictor, alpha)
%DELTA_METHOD_CI Pointwise confidence intervals by first-order propagation.
%   CI = DELTA_METHOD_CI(MODEL_FUN, PARAMS, COVARIANCE, PREDICTOR, ALPHA)
%   evaluates MODEL_FUN(PARAMS, PREDICTOR) and propagates the full parameter
%   covariance through a central-difference prediction Jacobian. The result
%   is an asymptotic, pointwise confidence interval for the fitted mean.
%
%   The returned structure contains estimate, standard_error, lower, upper,
%   confidence_level, and method fields.

    if nargin < 5 || isempty(alpha)
        alpha = 0.05;
    end

    validateattributes(alpha, {'numeric'}, {'scalar', 'real', '>', 0, '<', 1});
    validateattributes(fitted_params, {'numeric'}, {'vector', 'real', 'finite'});

    params = fitted_params(:).';
    n_params = numel(params);
    estimate = model_fun(params, predictor);
    estimate_shape = size(estimate);
    estimate_vector = estimate(:);

    ci = struct( ...
        'estimate', estimate, ...
        'standard_error', NaN(estimate_shape), ...
        'lower', NaN(estimate_shape), ...
        'upper', NaN(estimate_shape), ...
        'confidence_level', 1 - alpha, ...
        'method', 'Delta method (asymptotic pointwise)');

    if ~isequal(size(covariance), [n_params, n_params]) || ...
            any(~isfinite(covariance), 'all')
        return;
    end

    covariance = (covariance + covariance.') ./ 2;
    [eigenvectors, eigenvalues] = eig(covariance, 'vector');
    eigenvalues = max(real(eigenvalues), 0);
    covariance = real(eigenvectors * diag(eigenvalues) * eigenvectors.');

    prediction_jacobian = zeros(numel(estimate_vector), n_params);
    step_scale = eps(class(params)).^(1/3);

    for parameter_idx = 1:n_params
        step = step_scale .* max(abs(params(parameter_idx)), 1);
        params_plus = params;
        params_minus = params;
        params_plus(parameter_idx) = params_plus(parameter_idx) + step;
        params_minus(parameter_idx) = params_minus(parameter_idx) - step;

        prediction_plus = model_fun(params_plus, predictor);
        prediction_minus = model_fun(params_minus, predictor);
        prediction_jacobian(:, parameter_idx) = ...
            (prediction_plus(:) - prediction_minus(:)) ./ (2 .* step);
    end

    propagated_variance = sum((prediction_jacobian * covariance) .* ...
                              prediction_jacobian, 2);
    propagated_variance = max(real(propagated_variance), 0);
    standard_error = sqrt(propagated_variance);

    % Standard-normal critical value avoids adding Statistics Toolbox as a
    % runtime dependency. This is the usual large-sample delta-method band.
    critical_value = -sqrt(2) .* erfcinv(2 .* (1 - alpha ./ 2));
    margin = critical_value .* standard_error;

    ci.standard_error = reshape(standard_error, estimate_shape);
    ci.lower = reshape(estimate_vector - margin, estimate_shape);
    ci.upper = reshape(estimate_vector + margin, estimate_shape);
end
