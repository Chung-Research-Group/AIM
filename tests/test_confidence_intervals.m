function tests = test_confidence_intervals
%TEST_CONFIDENCE_INTERVALS Tests for delta-method uncertainty propagation.
    tests = functiontests(localfunctions);
end

function setupOnce(testCase)
    repo_root = fileparts(fileparts(mfilename('fullpath')));
    src_dir = fullfile(repo_root, 'src');
    addpath(src_dir);
    testCase.TestData.src_dir = src_dir;
end

function teardownOnce(testCase)
    rmpath(testCase.TestData.src_dir);
end

function testDeltaMethodMatchesLinearAnalyticVariance(testCase)
    model = @(params, x) params(1) + params(2) .* x;
    params = [2.0, 0.5];
    covariance = [0.04, -0.012; -0.012, 0.01];
    x = [0; 1; 3];

    ci = delta_method_ci(model, params, covariance, x, 0.05);

    design = [ones(size(x)), x];
    expected_variance = sum((design * covariance) .* design, 2);
    verifyEqual(testCase, ci.standard_error, sqrt(expected_variance), ...
        'RelTol', 2e-6, 'AbsTol', 1e-10);
    verifyEqual(testCase, ci.estimate, model(params, x), 'AbsTol', eps);
    verifyGreaterThan(testCase, ci.upper, ci.estimate);
    verifyLessThan(testCase, ci.lower, ci.estimate);
end

function testOffDiagonalCovarianceChangesBandWidth(testCase)
    model = @(params, x) params(1) + params(2) .* x;
    params = [1.0, 1.0];
    x = 2.0;
    covariance_full = [1.0, -0.45; -0.45, 0.25];
    covariance_diagonal = diag(diag(covariance_full));

    ci_full = delta_method_ci(model, params, covariance_full, x, 0.05);
    ci_diagonal = delta_method_ci(model, params, covariance_diagonal, x, 0.05);

    verifyNotEqual(testCase, ci_full.standard_error, ci_diagonal.standard_error);
    verifyLessThan(testCase, ci_full.standard_error, ci_diagonal.standard_error);
end

function testInvalidCovarianceReturnsUnavailableBand(testCase)
    model = @(params, x) params(1) .* x;
    ci = delta_method_ci(model, 2.0, NaN, [1; 2], 0.05);

    verifyTrue(testCase, all(isnan(ci.lower), 'all'));
    verifyTrue(testCase, all(isnan(ci.upper), 'all'));
    verifyEqual(testCase, ci.estimate, [2; 4]);
end

function testIsothermCovarianceReconstructionUsesStoredStandardErrors(testCase)
    pressure = linspace(0.1, 2.0, 20)';
    model = @(params, x) params(1) + params(2) .* x;
    params = [1.2, 0.8];
    loading = model(params, pressure) + 0.01 .* sin((1:numel(pressure))');

    design = [ones(size(pressure)), pressure];
    residual = loading - model(params, pressure);
    dof = numel(loading) - numel(params);
    expected_covariance = (sum(residual.^2) ./ dof) .* pinv(design.' * design);

    fit_struct = struct( ...
        'fun', model, ...
        'fitted_params', params, ...
        'std_error', sqrt(diag(expected_covariance)), ...
        'p_sat', false);

    uncertainty = reconstruct_isotherm_covariance( ...
        fit_struct, pressure, loading, 1.0);

    verifyEqual(testCase, uncertainty.weighting, 'Uniform');
    verifyTrue(testCase, uncertainty.identifiable);
    verifyEqual(testCase, uncertainty.covariance, expected_covariance, ...
        'RelTol', 5e-6, 'AbsTol', 1e-12);
end

function testVirialCovarianceReconstruction(testCase)
    loading = linspace(0.1, 2.0, 12)';
    temperatures = [280, 320];
    params = [0.7, -0.15, 90.0];
    num_a = 2;
    virial_fun = @(candidate, q, T, number_a) ...
        polyval(flip(candidate(1:number_a)), q) + candidate(3) ./ T;

    pressure = NaN(numel(loading), numel(temperatures));
    loading_matrix = repmat(loading, 1, numel(temperatures));
    observed_log_pressure = [];
    design = [];
    for idx = 1:numel(temperatures)
        log_pressure = virial_fun(params, loading, ...
            repmat(temperatures(idx), size(loading)), num_a) + ...
            0.002 .* cos((1:numel(loading))' + idx);
        pressure(:, idx) = exp(log_pressure);
        observed_log_pressure = [observed_log_pressure; log_pressure]; %#ok<AGROW>
        design = [design; ones(size(loading)), loading, ...
            repmat(1 ./ temperatures(idx), size(loading))]; %#ok<AGROW>
    end

    prediction = design * params(:);
    residual = observed_log_pressure - prediction;
    dof = numel(residual) - numel(params);
    expected_covariance = (sum(residual.^2) ./ dof) .* pinv(design.' * design);

    fit_struct = struct( ...
        'fun', virial_fun, ...
        'fitted_params', params, ...
        'std_error', sqrt(diag(expected_covariance)));

    uncertainty = reconstruct_virial_covariance( ...
        fit_struct, pressure, loading_matrix, temperatures, num_a);

    verifyEqual(testCase, uncertainty.weighting, 'Uniform');
    verifyTrue(testCase, uncertainty.identifiable);
    verifyEqual(testCase, uncertainty.covariance, expected_covariance, ...
        'RelTol', 2e-5, 'AbsTol', 1e-11);
end
