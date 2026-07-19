function tests = test_buildfile
%TEST_BUILDFILE Regression tests for the AIM build plan configuration.
    tests = functiontests(localfunctions);
end

function setupOnce(testCase)
    repo_root = fileparts(fileparts(mfilename('fullpath')));
    testCase.TestData.repo_root = repo_root;
    addpath(repo_root);
end

function teardownOnce(testCase)
    rmpath(testCase.TestData.repo_root);
end

function testTestTaskSourceFilesMatchesExpectedCoverageList(testCase)
    % Pins the exact list of files whose coverage is tracked by
    % `buildtool test`, guarding against the newly added coupling module
    % being dropped from, or misspelled within, that list.
    plan = buildfile();
    test_task = plan("test");

    expected = ["src/WENO.m"
                "src/Isotherm_functions.m"
                "src/IAST_func_NR.m"
                "src/IAST_func_mix_pred_app_NR.m"
                "src/Jacobian.m"
                "src/couple_pressure_temperature_rates.m"];

    verifyEqual(testCase, test_task.SourceFiles(:), expected(:));
end

function testTestTaskCoversNewPressureTemperatureCouplingModule(testCase)
    plan = buildfile();
    test_task = plan("test");

    verifyTrue(testCase, ...
        any(test_task.SourceFiles == "src/couple_pressure_temperature_rates.m"));
end

function testTestTaskSourceFilesAllExistOnDisk(testCase)
    plan = buildfile();
    test_task = plan("test");
    repo_root = testCase.TestData.repo_root;

    for idx = 1:numel(test_task.SourceFiles)
        source_file = fullfile(repo_root, test_task.SourceFiles(idx));
        verifyTrue(testCase, isfile(source_file), ...
            sprintf('Missing covered source file: %s', test_task.SourceFiles(idx)));
    end
end

function testTestTaskRunsWithStrictFailureMode(testCase)
    plan = buildfile();
    test_task = plan("test");

    verifyTrue(testCase, test_task.Strict);
end