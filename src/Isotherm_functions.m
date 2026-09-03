function solution = Isotherm_functions(num_comp, iso_params, P, T, T_flag, calc_mode)
%ISOTHERM_FUNCTIONS Evaluate pure-component loading or reduced grand potential.
%   calc_mode = 1 returns loading; calc_mode = 0 returns reduced grand potential.
%
%% Extracting isotherm model flags for components 
    iso_type_flag = iso_params(1, 1:num_comp);
%
%% Non-isothermal partial pressures based on CC-equation
    R = 8.31446261815324;
    if T_flag
        dH = iso_params(end-1, 1:num_comp);
        T_ref = iso_params(end, 1:num_comp);
        P = P .* exp(-dH ./ R .* (1 ./ T - 1 ./ T_ref));
    end
%
%% Loading OR reduced grand potential calculation
    solution = zeros(size(P), 'like', P);    
    % Loop over the components
    for i = 1:num_comp
        flag = iso_type_flag(i);
        p = P(:, i);

        switch flag
            case {1, 2} % SS/DS Langmuir
                q1 = iso_params(2, i); b1 = iso_params(3, i);
                q2 = iso_params(4, i); b2 = iso_params(5, i);
                if calc_mode == 1
                    solution(:, i) = q1 .* b1 .* p ./ (1 + b1 .* p) + ...
                                     q2 .* b2 .* p ./ (1 + b2 .* p);
                else
                    solution(:, i) = q1 .* log1p(b1 .* p) + q2 .* log1p(b2 .* p);
                end

            case {3, 4} % SS/DS Langmuir-Freundlich
                q1 = iso_params(2, i); b1 = iso_params(3, i); n1 = iso_params(4, i);
                q2 = iso_params(5, i); b2 = iso_params(6, i); n2 = iso_params(7, i);
                
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
                q = iso_params(2, i); b = iso_params(3, i); c = iso_params(4, i);
                denominator = 1 + b .* p + c .* p.^2;
                if calc_mode == 1
                    solution(:, i) = q .* (b .* p + 2 .* c .* p.^2) ./ denominator;
                else
                    solution(:, i) = q .* log(denominator);
                end

            case 6 % Temkin approximation
                m = iso_params(2, i); b = iso_params(3, i); theta = iso_params(4, i);
                lang_term = b .* p ./ (1 + b .* p);
                if calc_mode == 1
                    solution(:, i) = m .* lang_term + ...
                                     m .* theta .* lang_term.^2 .* (lang_term - 1);
                else
                    solution(:, i) = m .* (log1p(b .* p) - 0.5 .* theta .* lang_term.^2);
                end

            case 7 % BET
                m = iso_params(2, i); b_surface = iso_params(3, i); b_layers = iso_params(4, i);
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
                m = iso_params(2, i); b = iso_params(3, i); n = iso_params(4, i);
                term = (b .* p).^(1 ./ n);
                if calc_mode == 1
                    solution(:, i) = m .* term ./ (1 + term);
                else
                    solution(:, i) = m .* n .* log1p(term);
                end

            case 9 % Toth
                m = iso_params(2, i); b = iso_params(3, i); n = iso_params(4, i);
                temp = b .* p;
                theta1 = temp ./ (1 + temp.^n).^(1 ./ n);
                if calc_mode == 1
                    solution(:, i) = m .* theta1;
                else
                    solution(:, i) = toth_grand_potential(temp, m, n);
                end
        end
    end
    %
    %% Validation
    if any(~isfinite(solution), 'all')
        error('AIM:Isotherm:NonFiniteResult', ...
              'Isotherm evaluation produced NaN or Inf. Check parameters and pressure domain.');
    end
    %
    %% Miscellenous functions
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
%
end
