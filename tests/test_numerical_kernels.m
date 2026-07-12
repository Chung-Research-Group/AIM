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

function testLangmuirGrandPotentialDerivative(testCase)
    % q(P) = d psi / d ln(P) = P * d psi / dP.
    iso = zeros(7, 1);
    iso(1) = 1;
    iso(2) = 2.4;
    iso(3) = 0.8;
    pressure = logspace(-5, 2, 30)';

    loading = Isotherm_functions(1, iso, pressure, 298.15, 0, 1);

    h = 1e-6;
    psi_plus = Isotherm_functions(1, iso, pressure .* exp(h), 298.15, 0, 0);
    psi_minus = Isotherm_functions(1, iso, pressure .* exp(-h), 298.15, 0, 0);
    numerical_derivative = (psi_plus - psi_minus) ./ (2 * h);

    verifyEqual(testCase, numerical_derivative, loading, ...
                'RelTol', 2e-7, 'AbsTol', 1e-10);
end

function testUnsupportedIsothermFailsExplicitly(testCase)
    iso = zeros(7, 1);
    iso(1) = 10;

    verifyError(testCase, ...
        @() Isotherm_functions(1, iso, 1.0, 298.15, 0, 1), ...
        'AIM:Isotherm:UnsupportedModel');
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

    iso = zeros(7, 1);
    iso(1) = 1;
    iso(2) = 2.4;
    iso(3) = 0.8;
    pressure = logspace(-5, 1, 20)';
    composition = ones(size(pressure));
    temperature = 298.15 .* ones(size(pressure));

    expected = Isotherm_functions(1, iso, pressure, temperature, 0, 1);
    actual = IAST_func_NR(1, iso, pressure, composition, temperature, 0, []);

    verifyEqual(testCase, actual, expected, 'RelTol', 1e-12, 'AbsTol', 1e-12);
end
