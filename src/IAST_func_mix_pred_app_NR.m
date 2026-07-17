function partial_loadings = IAST_func_mix_pred_app_NR( ...
        num_components, isotherm_params_array, Pressure, ...
        gas_phase_mol_fraction, T, ~)
%IAST_FUNC_MIX_PRED_APP_NR Evaluate MixPred IAST with the shared solver.
%   Keeping a second Newton implementation caused the GUI and BreakLab to
%   apply different convergence criteria.  This adapter computes the local
%   initialization data, calls the validated shared kernel, and restores the
%   global cache state on exit.

    global cached_p0 Henry_Coeff cached_p0_local;

    previous_cached_p0 = cached_p0;
    previous_henry = Henry_Coeff;
    previous_cached_p0_local = cached_p0_local;

    try
        pressure_range = [1e-7; 1e-6];
        component_pressures = repmat(pressure_range, 1, num_components);
        low_pressure_loading = Isotherm_functions( ...
            num_components, isotherm_params_array, component_pressures, ...
            T, 0, 1);
        Henry_Coeff = diff(low_pressure_loading, 1, 1) ./ diff(pressure_range);

        cached_p0 = [];
        cached_p0_local = [];
        partial_loadings = IAST_func_NR( ...
            num_components, isotherm_params_array, Pressure, ...
            gas_phase_mol_fraction, T, 0, []);
    catch solver_error
        cached_p0 = previous_cached_p0;
        Henry_Coeff = previous_henry;
        cached_p0_local = previous_cached_p0_local;
        rethrow(solver_error);
    end

    cached_p0 = previous_cached_p0;
    Henry_Coeff = previous_henry;
    cached_p0_local = previous_cached_p0_local;
end
