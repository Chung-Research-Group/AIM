function solution = Isotherm_functions(num_comp, iso_params, P, T, T_flag, calc_mode)
%ISOTHERM_FUNCTIONS Evaluate pure-component loading or reduced grand potential.
%   calc_mode = 1 returns loading; calc_mode = 0 returns reduced grand potential.

    validateattributes(num_comp, {'numeric'}, {'scalar', 'integer', 'positive'});
    validateattributes(P, {'numeric'}, {'2d', 'real', 'finite', 'nonnegative'});
    validateattributes(T_flag, {'numeric', 'logical'}, {'scalar'});
    validateattributes(calc_mode, {'numeric', 'logical'}, {'scalar'});

    if size(P, 2) < num_comp
        error('AIM:Isotherm:PressureShape', ...
              'P must contain at least num_comp columns.');
    end
    if size(iso_params, 2) < num_comp || size(iso_params, 1) < 2
        error('AIM:Isotherm:ParameterShape', ...
              'iso_params has insufficient rows or columns.');
    end
    if ~(calc_mode == 0 || calc_mode == 1)
        error('AIM:Isotherm:InvalidCalculationMode', ...
              'calc_mode must be 0 (grand potential) or 1 (loading).');
    end

    iso_type_flag = iso_params(1, 1:num_comp);
    supported_flags = 1:9;
    unsupported = ~ismember(iso_type_flag, supported_flags);
    if any(unsupported)
        bad_flags = unique(iso_type_flag(unsupported));
        error('AIM:Isotherm:UnsupportedModel', ...
              'Unsupported isotherm model flag(s): %s. IAST/BreakLab support flags 1-9 only.', ...
              mat2str(bad_flags));
    end

    R = 8.31446261815324;
    if T_flag
        validateattributes(T, {'numeric'}, {'real', 'finite', 'positive'});
        if size(iso_params, 1) < 3
            error('AIM:Isotherm:MissingTemperatureParameters', ...
                  'Temperature-dependent evaluation requires dH and T_ref rows.');
        end
        dH = iso_params(end-1, 1:num_comp);
        T_ref = iso_params(end, 1:num_comp);
        if any(~isfinite(dH)) || any(~isfinite(T_ref)) || any(T_ref <= 0)
            error('AIM:Isotherm:InvalidTemperatureParameters', ...
                  'dH must be finite and T_ref must be finite and positive.');
        end
        P = P .* exp(-dH ./ R .* (1 ./ T - 1 ./ T_ref));
    end

    solution = zeros(size(P), 'like', P);

    for i = 1:num_comp
        flag = iso_type_flag(i);
        p = P(:, i);

        switch flag
            case {1, 2} % SS/DS Langmuir
                require_rows(5, flag);
                q1 = iso_params(2, i); b1 = iso_params(3, i);
                q2 = iso_params(4, i); b2 = iso_params(5, i);
                require_nonnegative([q1, b1, q2, b2], flag);
                if calc_mode == 1
                    solution(:, i) = q1 .* b1 .* p ./ (1 + b1 .* p) + ...
                                     q2 .* b2 .* p ./ (1 + b2 .* p);
                else
                    solution(:, i) = q1 .* log1p(b1 .* p) + q2 .* log1p(b2 .* p);
                end

            case {3, 4} % SS/DS Langmuir-Freundlich
                require_rows(7, flag);
                q1 = iso_params(2, i); b1 = iso_params(3, i); n1 = iso_params(4, i);
                q2 = iso_params(5, i); b2 = iso_params(6, i); n2 = iso_params(7, i);
                require_nonnegative([q1, b1, q2, b2], flag);
                if n1 <= 0 || (q2 > 0 && n2 <= 0)
                    error('AIM:Isotherm:InvalidExponent', ...
                          'Langmuir-Freundlich exponents must be positive.');
                end
                term1 = b1 .* p.^n1;
                if q2 == 0
                    term2 = zeros(size(p), 'like', p);
                else
                    term2 = b2 .* p.^n2;
                end
                if calc_mode == 1
                    solution(:, i) = q1 .* term1 ./ (1 + term1) + ...
                                     q2 .* term2 ./ (1 + term2);
                else
                    solution(:, i) = q1 ./ n1 .* log1p(term1);
                    if q2 > 0
                        solution(:, i) = solution(:, i) + q2 ./ n2 .* log1p(term2);
                    end
                end

            case 5 % Quadratic
                require_rows(4, flag);
                q = iso_params(2, i); b = iso_params(3, i); c = iso_params(4, i);
                require_nonnegative([q, b, c], flag);
                denominator = 1 + b .* p + c .* p.^2;
                if calc_mode == 1
                    solution(:, i) = q .* (b .* p + 2 .* c .* p.^2) ./ denominator;
                else
                    solution(:, i) = q .* log(denominator);
                end

            case 6 % Temkin approximation
                require_rows(4, flag);
                m = iso_params(2, i); b = iso_params(3, i); theta = iso_params(4, i);
                require_nonnegative([m, b], flag);
                lang_term = b .* p ./ (1 + b .* p);
                if calc_mode == 1
                    solution(:, i) = m .* lang_term + ...
                                     m .* theta .* lang_term.^2 .* (lang_term - 1);
                else
                    solution(:, i) = m .* (log1p(b .* p) - 0.5 .* theta .* lang_term.^2);
                end

            case 7 % BET
                require_rows(4, flag);
                m = iso_params(2, i); b_surface = iso_params(3, i); b_layers = iso_params(4, i);
                require_nonnegative([m, b_surface, b_layers], flag);
                if any(b_layers .* p >= 1)
                    error('AIM:Isotherm:BETDomain', ...
                          'BET evaluation requires b_layers * P < 1.');
                end
                if calc_mode == 1
                    solution(:, i) = m .* b_surface .* p ./ ...
                        ((1 - b_layers .* p) .* ...
                         (1 - b_layers .* p + b_surface .* p));
                else
                    solution(:, i) = m .* log((1 + b_surface .* p - b_layers .* p) ./ ...
                                              (1 - b_layers .* p));
                end

            case 8 % Sips
                require_rows(4, flag);
                m = iso_params(2, i); b = iso_params(3, i); n = iso_params(4, i);
                require_nonnegative([m, b], flag);
                if n <= 0
                    error('AIM:Isotherm:InvalidExponent', 'Sips exponent must be positive.');
                end
                term = (b .* p).^(1 ./ n);
                if calc_mode == 1
                    solution(:, i) = m .* term ./ (1 + term);
                else
                    solution(:, i) = m .* n .* log1p(term);
                end

            case 9 % Toth
                require_rows(4, flag);
                m = iso_params(2, i); b = iso_params(3, i); n = iso_params(4, i);
                require_nonnegative([m, b], flag);
                if n <= 0
                    error('AIM:Isotherm:InvalidExponent', 'Toth exponent must be positive.');
                end
                temp = b .* p;
                theta1 = temp ./ (1 + temp.^n).^(1 ./ n);
                if calc_mode == 1
                    solution(:, i) = m .* theta1;
                else
                    solution(:, i) = toth_grand_potential(temp, m, n);
                end
        end
    end

    if any(~isfinite(solution), 'all')
        error('AIM:Isotherm:NonFiniteResult', ...
              'Isotherm evaluation produced NaN or Inf. Check parameters and pressure domain.');
    end

    function require_rows(last_row, model_flag)
        if size(iso_params, 1) < last_row
            error('AIM:Isotherm:MissingParameters', ...
                  'Model flag %d requires parameters through row %d.', ...
                  model_flag, last_row);
        end
    end

    function require_nonnegative(values, model_flag)
        if any(~isfinite(values)) || any(values < 0)
            error('AIM:Isotherm:InvalidParameters', ...
                  'Model flag %d requires finite nonnegative capacity/affinity parameters.', ...
                  model_flag);
        end
    end

    function grand_potential = toth_grand_potential(reduced_pressure, capacity, exponent)
        % q = d(psi)/d(log(P)).  With
        % theta=t/(1+t^n)^(1/n), direct integration gives
        % psi/m = sum_{k=0}^inf theta^(n*k+1)/(n*k+1).
        % The former fixed 100-term rearrangement lost accuracy as theta
        % approached one.  Sum the well-conditioned leading terms and use
        % an Euler-Maclaurin tail only where it is needed.
        grand_potential = zeros(size(reduced_pressure), 'like', reduced_pressure);
        positive = reduced_pressure > 0;
        if ~any(positive)
            return;
        end

        t = reduced_pressure(positive);
        if exponent == 1
            grand_potential(positive) = capacity .* log1p(t);
            return;
        end

        log_z = zeros(size(t), 'like', t);
        small_t = t <= 1;
        log_z(small_t) = exponent .* log(t(small_t)) - ...
            log1p(t(small_t).^exponent);
        log_z(~small_t) = -log1p(t(~small_t).^(-exponent));
        theta = exp(log_z ./ exponent);

        n_terms = 160;
        series_sum = zeros(size(t), 'like', t);
        for k = 0:n_terms-1
            series_sum = series_sum + ...
                theta .* exp(k .* log_z) ./ (exponent .* k + 1);
        end

        % For z <= 0.8, the omitted tail is below double-precision scale
        % after 160 terms for practical Toth exponents.  Near z=1, use the
        % first Euler-Maclaurin corrections to retain high-pressure accuracy
        % without an unbounded iteration count.
        tail_rows = log_z > log(0.8);
        if any(tail_rows)
            local_log_z = log_z(tail_rows);
            lambda = -local_log_z;
            a = 1 ./ exponent;
            f_k = theta(tail_rows) .* exp(n_terms .* local_log_z) ./ ...
                (exponent .* n_terms + 1);
            integral_tail = theta(tail_rows) .* exp(lambda .* a) ./ exponent .* ...
                expint(lambda .* (n_terms + a));
            derivative_f = f_k .* ...
                (local_log_z - exponent ./ (exponent .* n_terms + 1));
            series_sum(tail_rows) = series_sum(tail_rows) + ...
                integral_tail + 0.5 .* f_k - derivative_f ./ 12;
        end

        grand_potential(positive) = capacity .* series_sum;
    end
end
