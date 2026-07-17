function partial_loadings = IAST_func_NR(num_components, isotherm_params_array, Pressure, gas_phase_mol_fraction, T_array, T_flag, ~)
    global cached_p0;
    global Henry_Coeff;
    global cached_p0_local;

    validateattributes(num_components, {'numeric'}, ...
        {'scalar', 'integer', 'positive'});
    validateattributes(Pressure, {'numeric'}, ...
        {'column', 'real', 'finite', 'nonnegative'});
    validateattributes(gas_phase_mol_fraction, {'numeric'}, ...
        {'2d', 'real', 'finite', 'nonnegative'});
    if size(gas_phase_mol_fraction, 1) ~= size(Pressure, 1) || ...
            size(gas_phase_mol_fraction, 2) ~= num_components
        error('AIM:IAST:CompositionShape', ...
              'gas_phase_mol_fraction must be nPressure-by-num_components.');
    end
    if any(sum(gas_phase_mol_fraction, 2) > 1 + 1e-10)
        error('AIM:IAST:InvalidComposition', ...
              'Adsorbing-component mole fractions cannot sum to more than one.');
    end
    validateattributes(T_array, {'numeric'}, {'real', 'finite', 'positive'});
    if isscalar(T_array)
        T_array = repmat(T_array, size(Pressure));
    elseif ~isequal(size(T_array), size(Pressure))
        error('AIM:IAST:TemperatureShape', ...
              'T_array must be scalar or have the same size as Pressure.');
    end

    size_of_pressure_vector = size(Pressure, 1);
    partial_pressures = Pressure .* gas_phase_mol_fraction;
    fictitious_pressure = zeros(size_of_pressure_vector, num_components);

    pressure_tol = 1e-13;
    loading_floor = 1e-15;
    max_iter = 100;

    if isempty(Henry_Coeff) || numel(Henry_Coeff) < num_components
        error('AIM:IAST:InvalidHenryCoefficient', ...
              'Henry coefficients must be initialized for every adsorbing component.');
    end

    active_henry = Henry_Coeff(1:num_components);
    if any(~isfinite(active_henry)) || any(active_henry <= 0)
        error('AIM:IAST:InvalidHenryCoefficient', ...
              'Henry coefficients must be finite and strictly positive.');
    end

    sum_Hcoeff = sum(partial_pressures .* active_henry, 2);
    H_guess = zeros(size(partial_pressures));
    for component_idx = 1:num_components
        H_guess(:, component_idx) = min(sum_Hcoeff ./ active_henry(component_idx), Pressure);
    end
    H_guess = max(H_guess, pressure_tol);

    use_cache = true;
    if use_cache && ~isempty(cached_p0) && isequal(size(cached_p0), size(H_guess))
        guess = cached_p0;
        invalid_cache = ~isfinite(guess) | guess <= 0;
        guess(invalid_cache) = H_guess(invalid_cache);
    else
        guess = H_guess;
    end

    for node_idx = 1:size_of_pressure_vector
        if sum(partial_pressures(node_idx, :), 2) <= pressure_tol
            continue;
        end

        active_components = find(partial_pressures(node_idx, :) > pressure_tol);
        if isempty(active_components)
            continue;
        elseif isscalar(active_components)
            fictitious_pressure(node_idx, active_components) = ...
                partial_pressures(node_idx, active_components);
        else
            try
                solution = NR_function( ...
                    guess(node_idx, active_components), ...
                    isotherm_params_array(:, active_components), ...
                    partial_pressures(node_idx, active_components), ...
                    T_array(node_idx, 1), ...
                    length(active_components), ...
                    T_flag, ...
                    max_iter);
                fictitious_pressure(node_idx, active_components) = solution;
            catch solver_error
                error('AIM:IAST:ConvergenceFailure', ...
                      'IAST failed at pressure node %d: %s', ...
                      node_idx, solver_error.message);
            end
        end
    end

    cached_p0_local = fictitious_pressure;

    loading_array = Isotherm_functions( ...
        num_components, isotherm_params_array, fictitious_pressure, T_array, T_flag, 1);
    loading_array = max(loading_array, loading_floor);

    safe_fictitious_pressure = max(fictitious_pressure, loading_floor);
    adsorbed_mole_fractions = partial_pressures ./ safe_fictitious_pressure;

    inverse_loading = sum(adsorbed_mole_fractions ./ loading_array, 2);
    loading_total = zeros(size(inverse_loading));
    valid_rows = isfinite(inverse_loading) & inverse_loading > 0;
    loading_total(valid_rows) = 1 ./ inverse_loading(valid_rows);

    partial_loadings = adsorbed_mole_fractions .* loading_total;
    if any(~isfinite(partial_loadings), 'all')
        error('AIM:IAST:NonFiniteLoading', ...
              'IAST produced a non-finite partial loading.');
    end

    function converged_sol = NR_function(p0_guess, local_iso_params, partial_pressure, T, ncomp, local_T_flag, local_max_iter)
        grand_potential_tol = 1e-10;
        mole_fraction_tol = 1e-10;
        pressure_floor = 1e-15;

        partial_pressure = partial_pressure(:);
        p0_guess = max(p0_guess(:), pressure_floor);

        G_vec = zeros(ncomp, 1);
        Jac = zeros(ncomp, ncomp);
        delta = zeros(ncomp, 1);
        update_guess = p0_guess;

        for iteration = 1:local_max_iter
            grand_potentials = Isotherm_functions( ...
                ncomp, local_iso_params, p0_guess', T, local_T_flag, 0);
            loadings = Isotherm_functions( ...
                ncomp, local_iso_params, p0_guess', T, local_T_flag, 1);
            grand_potentials = grand_potentials(:);
            loadings = loadings(:);

            if any(~isfinite(grand_potentials)) || any(~isfinite(loadings)) || ...
                    any(loadings <= 0)
                error('Non-finite or non-positive isotherm values during Newton iteration %d.', iteration);
            end

            G_vec(1:end-1) = grand_potentials(1:end-1) - grand_potentials(end);
            G_vec(end) = 1 - sum(partial_pressure ./ p0_guess);

            Jac(:) = 0;
            diag_idx = sub2ind([ncomp, ncomp], 1:ncomp-1, 1:ncomp-1);
            Jac(diag_idx) = loadings(1:end-1) ./ p0_guess(1:end-1);
            Jac(1:ncomp-1, ncomp) = -loadings(end) / p0_guess(end);
            Jac(ncomp, :) = (partial_pressure ./ p0_guess.^2)';

            diag_values = Jac(diag_idx);
            if any(~isfinite(diag_values)) || any(abs(diag_values) < eps)
                error('Singular IAST Jacobian during Newton iteration %d.', iteration);
            end

            elimination = Jac(ncomp, 1:ncomp-1)' ./ diag_values(:);
            Jac(ncomp, ncomp) = Jac(ncomp, ncomp) - ...
                sum(elimination .* Jac(1:ncomp-1, ncomp));
            G_vec(ncomp) = G_vec(ncomp) - ...
                sum(elimination .* G_vec(1:ncomp-1));

            if ~isfinite(Jac(ncomp, ncomp)) || abs(Jac(ncomp, ncomp)) < eps
                error('Singular reduced IAST Jacobian during Newton iteration %d.', iteration);
            end

            delta(ncomp) = G_vec(ncomp) / Jac(ncomp, ncomp);
            delta(1:ncomp-1) = ...
                (G_vec(1:ncomp-1) - delta(ncomp) .* Jac(1:ncomp-1, ncomp)) ./ diag_values(:);

            step_scale = 1.0;
            while any(p0_guess - step_scale .* delta <= pressure_floor) && step_scale > 2^-20
                step_scale = step_scale / 2;
            end
            if step_scale <= 2^-20
                error('Newton step could not maintain positive fictitious pressures.');
            end

            update_guess = p0_guess - step_scale .* delta;
            updated_potential = Isotherm_functions( ...
                ncomp, local_iso_params, update_guess', T, local_T_flag, 0);
            updated_potential = updated_potential(:);

            gp_residual = max(abs(updated_potential - mean(updated_potential)));
            mf_residual = abs(sum(partial_pressure ./ update_guess) - 1);

            if gp_residual <= grand_potential_tol && mf_residual <= mole_fraction_tol
                converged_sol = reshape(update_guess, 1, []);
                return;
            end

            p0_guess = update_guess;
        end

        error(['Maximum Newton iterations reached. Final spreading-potential residual = %.3e, ' ...
               'mole-fraction residual = %.3e.'], gp_residual, mf_residual);
    end
end
