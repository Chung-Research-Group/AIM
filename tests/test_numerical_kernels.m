function tests = test_numerical_kernels
%TEST_NUMERICAL_KERNELS Regression tests for core AIM numerical routines.
    tests = functiontests(localfunctions);
end

function setupOnce(testCase)
    repo_root = fileparts(fileparts(mfilename('fullpath')));
    src_dir = fullfile(repo_root, 'src');
    addpath(src_dir);
    testCase.TestData.src_dir = src_dir;
end

function teardown(~)
    clear global cached_p0 Henry_Coeff cached_p0_local;
end

function teardownOnce(testCase)
    rmpath(testCase.TestData.src_dir);
    clear global cached_p0 Henry_Coeff cached_p0_local;
end

function testWenoPreservesConstantField(testCase)
    flux = 3.25 .* ones(12, 3);

    upwind = WENO(flux, 'upwind');
    downwind = WENO(flux, 'downwind');

    verifyEqual(testCase, upwind, 3.25 .* ones(11, 3), 'AbsTol', 10 * eps);
    verifyEqual(testCase, downwind, 3.25 .* ones(11, 3), 'AbsTol', 10 * eps);
end

function testWenoReconstructsLinearInteriorExactly(testCase)
    flux = (0:11)';

    upwind = WENO(flux, 'upwind');
    downwind = WENO(flux, 'downwind');

    verifyEqual(testCase, upwind(3:10), (2.5:9.5)', 'AbsTol', 100 * eps);
    verifyEqual(testCase, downwind(2:9), (1.5:8.5)', 'AbsTol', 100 * eps);
end

function testWenoIsMirrorSymmetric(testCase)
    flux = [0.0; 0.2; 0.9; 0.1; -0.4; 0.7; 1.3; 0.4];

    downwind = WENO(flux, 'downwind');
    mirrored_upwind = flipud(WENO(flipud(flux), 'upwind'));

    verifyEqual(testCase, downwind, mirrored_upwind, ...
                'RelTol', 1e-13, 'AbsTol', 1e-13);
end

function testWenoAvoidsSignedEpsilonSingularity(testCase)
    % The former formulation used (delta_flux + 1e-10)^4 and became
    % singular when delta_flux was exactly -1e-10.
    flux = [0; 0; -1e-10; -2e-10; -3e-10; -4e-10];

    reconstructed = WENO(flux, 'upwind');

    verifyTrue(testCase, all(isfinite(reconstructed), 'all'));
end

function testWenoRejectsInvalidDirection(testCase)
    verifyError(testCase, @() WENO(ones(6, 1), 'sideways'), ...
                'AIM:WENO:InvalidFlowDirection');
end

function testWenoRejectsNonfiniteFlux(testCase)
    verifyError(testCase, @() WENO([0; 1; NaN; 2], 'upwind'), ...
                'AIM:WENO:NonFiniteFlux');
end

function testLangmuirGrandPotentialDerivative(testCase)
    % q(P) = d psi / d ln(P) = P * d psi / dP.
    iso = singleSiteLangmuir(2.4, 0.8);
    pressure = logspace(-5, 2, 30)';

    loading = Isotherm_functions(1, iso, pressure, 298.15, 0, 1);

    h = 1e-6;
    psi_plus = Isotherm_functions(1, iso, pressure .* exp(h), 298.15, 0, 0);
    psi_minus = Isotherm_functions(1, iso, pressure .* exp(-h), 298.15, 0, 0);
    numerical_derivative = (psi_plus - psi_minus) ./ (2 * h);

    verifyEqual(testCase, numerical_derivative, loading, ...
                'RelTol', 2e-7, 'AbsTol', 1e-10);
end

function testAllSupportedGrandPotentialsAreThermodynamicallyConsistent(testCase)
    % For every supported model, q(P) = d(psi)/d(log(P)).  This catches a
    % scientifically invalid IAST input even when the Newton solver itself
    % converges cleanly.
    models = {
        [1; 2.4; 0.8; 0; 0; 0; 0], ...
        [2; 1.2; 0.8; 0.9; 0.05; 0; 0], ...
        [3; 2.4; 0.8; 0.7; 0; 0; 0], ...
        [4; 1.2; 0.8; 0.7; 0.9; 0.05; 1.3], ...
        [5; 2.4; 0.8; 0.1; 0; 0; 0], ...
        [6; 2.4; 0.8; 0.4; 0; 0; 0], ...
        [7; 2.4; 0.8; 0.02; 0; 0; 0], ...
        [8; 2.4; 0.8; 1.4; 0; 0; 0], ...
        [9; 2.4; 0.8; 0.55; 0; 0; 0] ...
    };
    pressure = logspace(-6, 1, 35)';
    h = 1e-5;

    for model_idx = 1:numel(models)
        iso = models{model_idx};
        loading = Isotherm_functions(1, iso, pressure, 298.15, 0, 1);
        psi_plus = Isotherm_functions( ...
            1, iso, pressure .* exp(h), 298.15, 0, 0);
        psi_minus = Isotherm_functions( ...
            1, iso, pressure .* exp(-h), 298.15, 0, 0);
        numerical_derivative = (psi_plus - psi_minus) ./ (2 .* h);

        verifyEqual(testCase, numerical_derivative, loading, ...
            'RelTol', 2e-6, 'AbsTol', 2e-9, ...
            sprintf('Model flag %d violates q=dpsi/dlnP.', model_idx));
    end
end

function testTothGrandPotentialAtHighReducedPressure(testCase)
    iso = zeros(7, 1);
    iso(1:4) = [9; 2.4; 0.8; 0.5];
    pressure = logspace(2, 10, 25)';
    h = 1e-4;

    loading = Isotherm_functions(1, iso, pressure, 298.15, 0, 1);
    psi_plus = Isotherm_functions( ...
        1, iso, pressure .* exp(h), 298.15, 0, 0);
    psi_minus = Isotherm_functions( ...
        1, iso, pressure .* exp(-h), 298.15, 0, 0);

    verifyEqual(testCase, (psi_plus - psi_minus) ./ (2 .* h), loading, ...
        'RelTol', 2e-6, 'AbsTol', 2e-8);
end

function testLangmuirLowAndHighPressureLimits(testCase)
    capacity = 2.4;
    affinity = 0.8;
    iso = singleSiteLangmuir(capacity, affinity);
    pressure = [0; 1e-12; 1e12];

    loading = Isotherm_functions(1, iso, pressure, 298.15, 0, 1);

    verifyEqual(testCase, loading(1), 0, 'AbsTol', eps);
    verifyEqual(testCase, loading(2) ./ pressure(2), ...
                capacity .* affinity, 'RelTol', 1e-10);
    verifyEqual(testCase, loading(3), capacity, 'RelTol', 1e-10);
end

function testExothermicTemperatureScaling(testCase)
    iso = singleSiteLangmuir(2.4, 0.8);
    iso(end-1) = -20000;  % J/mol
    iso(end) = 298.15;    % K

    pressure = ones(3, 1);
    temperature = [250; 298.15; 350];
    loading = Isotherm_functions(1, iso, pressure, temperature, 1, 1);

    verifyGreaterThan(testCase, loading(1), loading(2));
    verifyGreaterThan(testCase, loading(2), loading(3));
end

function testUnsupportedIsothermFailsExplicitly(testCase)
    iso = zeros(7, 1);
    iso(1) = 10;

    verifyError(testCase, ...
        @() Isotherm_functions(1, iso, 1.0, 298.15, 0, 1), ...
        'AIM:Isotherm:UnsupportedModel');
end

function testNegativeIsothermParameterIsRejected(testCase)
    iso = singleSiteLangmuir(2.4, -0.8);

    verifyError(testCase, ...
        @() Isotherm_functions(1, iso, 1.0, 298.15, 0, 1), ...
        'AIM:Isotherm:InvalidParameters');
end

function testBetDomainIsChecked(testCase)
    iso = zeros(7, 1);
    iso(1) = 7;
    iso(2:4) = [1.0; 2.0; 0.5];

    verifyError(testCase, ...
        @() Isotherm_functions(1, iso, 2.0, 298.15, 0, 1), ...
        'AIM:Isotherm:BETDomain');
end

function testSingleComponentIastMatchesPureLoading(testCase)
    global cached_p0 Henry_Coeff cached_p0_local;
    cached_p0 = [];
    cached_p0_local = [];
    Henry_Coeff = 2.4 * 0.8;

    iso = singleSiteLangmuir(2.4, 0.8);
    pressure = logspace(-5, 1, 20)';
    composition = ones(size(pressure));
    temperature = 298.15 .* ones(size(pressure));

    expected = Isotherm_functions(1, iso, pressure, temperature, 0, 1);
    actual = IAST_func_NR(1, iso, pressure, composition, temperature, 0, []);

    verifyEqual(testCase, actual, expected, 'RelTol', 1e-12, 'AbsTol', 1e-12);
end

function testIdenticalBinaryIastPreservesComposition(testCase)
    global cached_p0 Henry_Coeff cached_p0_local;
    cached_p0 = [];
    cached_p0_local = [];

    capacity = 2.4;
    affinity = 0.8;
    Henry_Coeff = [capacity * affinity, capacity * affinity];

    pure_iso = singleSiteLangmuir(capacity, affinity);
    mixture_iso = repmat(pure_iso, 1, 2);
    pressure = logspace(-5, 2, 20)';
    composition = repmat([0.25, 0.75], numel(pressure), 1);
    temperature = 298.15 .* ones(size(pressure));

    pure_loading = Isotherm_functions( ...
        1, pure_iso, pressure, temperature, 0, 1);
    actual = IAST_func_NR( ...
        2, mixture_iso, pressure, composition, temperature, 0, []);
    expected = pure_loading .* composition;

    verifyEqual(testCase, actual, expected, ...
                'RelTol', 1e-10, 'AbsTol', 1e-12);
    verifyEqual(testCase, sum(actual, 2), pure_loading, ...
                'RelTol', 1e-10, 'AbsTol', 1e-12);
end

function testZeroPressureIastReturnsZero(testCase)
    global cached_p0 Henry_Coeff cached_p0_local;
    cached_p0 = [];
    cached_p0_local = [];
    Henry_Coeff = [2.4 * 0.8, 1.5 * 0.3];

    iso = [singleSiteLangmuir(2.4, 0.8), ...
           singleSiteLangmuir(1.5, 0.3)];
    pressure = zeros(4, 1);
    composition = repmat([0.4, 0.6], 4, 1);
    temperature = 298.15 .* ones(size(pressure));

    actual = IAST_func_NR( ...
        2, iso, pressure, composition, temperature, 0, []);

    verifyEqual(testCase, actual, zeros(4, 2), 'AbsTol', eps);
end

function testIastAcceptsScalarTemperature(testCase)
    global cached_p0 Henry_Coeff cached_p0_local;
    cached_p0 = [];
    cached_p0_local = [];
    Henry_Coeff = [2.4 * 0.8, 1.5 * 0.3];

    iso = [singleSiteLangmuir(2.4, 0.8), ...
           singleSiteLangmuir(1.5, 0.3)];
    pressure = logspace(-4, 1, 8)';
    composition = repmat([0.4, 0.6], numel(pressure), 1);

    scalar_temperature = IAST_func_NR( ...
        2, iso, pressure, composition, 298.15, 0, []);
    vector_temperature = IAST_func_NR( ...
        2, iso, pressure, composition, ...
        298.15 .* ones(size(pressure)), 0, []);

    verifyEqual(testCase, scalar_temperature, vector_temperature, ...
        'RelTol', 1e-12, 'AbsTol', 1e-12);
end

function testIastRejectsInvalidComposition(testCase)
    global cached_p0 Henry_Coeff cached_p0_local;
    cached_p0 = [];
    cached_p0_local = [];
    Henry_Coeff = [2.4 * 0.8, 1.5 * 0.3];
    iso = [singleSiteLangmuir(2.4, 0.8), ...
           singleSiteLangmuir(1.5, 0.3)];

    verifyError(testCase, ...
        @() IAST_func_NR(2, iso, 1.0, [0.7, 0.6], 298.15, 0, []), ...
        'AIM:IAST:InvalidComposition');
end

function testMixPredUsesSharedIastKernelAndRestoresGlobalState(testCase)
    global cached_p0 Henry_Coeff cached_p0_local;
    cached_p0 = 11;
    Henry_Coeff = 22;
    cached_p0_local = 33;

    iso = [singleSiteLangmuir(2.4, 0.8), ...
           singleSiteLangmuir(1.5, 0.3)];
    pressure = logspace(-4, 1, 8)';
    composition = repmat([0.4, 0.6], numel(pressure), 1);

    actual = IAST_func_mix_pred_app_NR( ...
        2, iso, pressure, composition, 298.15, []);

    local_henry = [2.4 * 0.8, 1.5 * 0.3];
    cached_p0 = [];
    cached_p0_local = [];
    Henry_Coeff = local_henry;
    expected = IAST_func_NR( ...
        2, iso, pressure, composition, 298.15, 0, []);

    verifyEqual(testCase, actual, expected, ...
        'RelTol', 1e-8, 'AbsTol', 1e-11);

    % Re-establish sentinels and verify that the adapter's cleanup path is
    % independent of the comparison call above.
    cached_p0 = 11;
    Henry_Coeff = 22;
    cached_p0_local = 33;
    IAST_func_mix_pred_app_NR( ...
        2, iso, pressure, composition, 298.15, []);
    verifyEqual(testCase, cached_p0, 11);
    verifyEqual(testCase, Henry_Coeff, 22);
    verifyEqual(testCase, cached_p0_local, 33);
end

function testNonisothermalJacobianHasCompleteStateSize(testCase)
    n_cells = 5;
    jacobian_pattern = Jacobian(n_cells, 1);
    state_size = 12 .* n_cells + 24;

    verifySize(testCase, jacobian_pattern, [state_size, state_size]);

    component_one_outlet = 2 .* n_cells + 4;
    component_one_interior = component_one_outlet - 1;
    expected_row = jacobian_pattern(component_one_interior, :);
    expected_row(component_one_outlet) = 0;
    verifyEqual(testCase, jacobian_pattern(component_one_outlet, :), expected_row);
end

function iso = singleSiteLangmuir(capacity, affinity)
    iso = zeros(7, 1);
    iso(1) = 1;
    iso(2) = capacity;
    iso(3) = affinity;
end
