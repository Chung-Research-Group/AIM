function partial_loadings = IAST_func_NR(num_components, isotherm_params_array, Pressure, gas_phase_mol_fraction, T_array, T_flag, ~)
    global cached_p0;
    global Henry_Coeff;
    global cached_p0_local;
    %% Parameter Unpacking and Decalarations
    size_of_pressure_vector = length(Pressure);
    partial_pressures = Pressure .* gas_phase_mol_fraction;
    fictitious_pressure = zeros(size_of_pressure_vector, num_components);
    %   
    %% Numerical Parameters
    tolerance = 1e-13;
    eps = 1e-15;
    max_iter = 100;
    %
    %% Initialization of Guess for IAST Equations Solution
    % If Initial Guess has not been provided in input arguments, 
    % it will be guessed proportional to pure compoenent loading.    
    % loading_array = Isotherm_functions(num_components, isotherm_params_array, partial_pressures, T_array, T_flag, 1);
    % loading_array(loading_array == 0) = eps;
    % initial_loading_guess = loading_array ./ (sum(loading_array, 2));
    % init_guess_p0 = partial_pressures ./ initial_loading_guess;
    % init_guess_p0(isnan(init_guess_p0)) = 0;
    % init_guess_p0(isinf(init_guess_p0)) = 0;
    sum_Hcoeff = sum(partial_pressures.*Henry_Coeff, 2);
    H_guess = zeros(size(partial_pressures));
    for l = 1:num_components
        H_guess(:, l) = min(sum_Hcoeff./Henry_Coeff(l), Pressure);
    end
    % end
    
    %% Whether to use cached ficititious pressures or not
    useCache = 1;
    if ~(isempty(cached_p0)) && useCache
        % idx_zeros = cached_p0<=0;  
        cached_p0(cached_p0<=0) = H_guess(cached_p0<=0);
        guess = cached_p0;
    else
        % guess = init_guess_p0(:, 1:num_components);         
        guess = H_guess;
    end
    %
    %% Calculate IAST loadings for all the nodes
    for k = 1:size_of_pressure_vector
        try
            % In case of no adsorbing component is present
            if sum(partial_pressures(k, :), 2) <= tolerance
                continue;
            else               
                % Find the idx of all adsorbing comp which have partial
                % pressure greater than tolerance value
                idx = find(partial_pressures(k, :)>tolerance);
             
                if isempty(idx)
                    % fictitious_pressure(k, idx) =  init_guess_p0(k, idx);
                    continue;

                elseif isscalar(idx)
                    % In case of single adsorbing component, the initial guess
                    % is the correct solution
                    % fictitious_pressure(k, idx) =  init_guess_p0(k, idx);
                    fictitious_pressure(k, idx) =  partial_pressures(k, idx);
                
                else
                    % Calculate IAST loadings 
                    solution = NR_function(guess(k, idx), isotherm_params_array(:, idx), partial_pressures(k, idx), T_array(k, 1), length(idx), T_flag, max_iter);
                    fictitious_pressure(k, idx) = solution;
                end
            end
        catch ME1
                error(['Failed to solve the IAST equations...!', ...
                       ' Error:%s'], ME1.message);    
        end
    end
    %
    %% Cache the solution
    % cached_p0 = fictitious_pressure;
    cached_p0_local = fictitious_pressure;
    %
    %% Calculate pure component loading at fictitious pressure
    loading_array = Isotherm_functions(num_components, isotherm_params_array, fictitious_pressure, T_array, T_flag, 1);
    loading_array(loading_array==0) = eps;
    %
    %% Calculate adsorbed mole fractions
    fictitious_pressure(fictitious_pressure==0) = eps;
    adsorbed_mole_fractions = partial_pressures ./ fictitious_pressure;
    % adsorbed_mole_fractions(isnan(adsorbed_mole_fractions)) = 0;
    % adsorbed_mole_fractions(isinf(adsorbed_mole_fractions)) = 0;
    %
    %% Calculate partial loadings 
    inverse_loading = sum(adsorbed_mole_fractions./loading_array, 2);
    loading_total = 1 ./ inverse_loading;
    loading_total(isnan(loading_total)) = 0;
    loading_total(isinf(loading_total)) = 0;
    
    partial_loadings = (adsorbed_mole_fractions .* loading_total);
    %
    %% Checks
    if any(isnan(partial_loadings))
        error(" NaN partial loadings...")
    end

    if any(isinf(partial_loadings))
        error("Inf partial loadings...")
    end
    %% NR Function
    function converged_sol = NR_function(p0_guess, isotherm_params_array, partial_pressure, T, ncomp, T_flag, max_iter)   
        % s = 1e-13;
        convergence_tol = 1e-10;
        G_vec = zeros(ncomp, 1);
        Jac = zeros(ncomp, ncomp);
        delta = zeros(ncomp, 1);
        grand_potentials = zeros(ncomp, 1);
        loadings = zeros(ncomp, 1);
        convergence_flag = 0;
        iter_count = 0;

        % Casting p0_guess and partial pressures as column vector
        partial_pressure = partial_pressure(:);
        p0_guess = p0_guess(:);
        update_guess = zeros(size(p0_guess));

        % NR-Iterations
        for m=1:max_iter
            % To avoid undefined values
            % p0_guess(p0_guess<s) = s;
            % p0_guess(p0_guess<0) = s;
    
            grand_potentials(:) = Isotherm_functions(ncomp, isotherm_params_array, p0_guess', T, T_flag, 0);         
            loadings(:) = Isotherm_functions(ncomp, isotherm_params_array, p0_guess', T, T_flag, 1); 
            
            % Residual vector
            G_vec(1:end-1, 1) = grand_potentials(1:ncomp-1, 1) - grand_potentials(ncomp, 1);       
            G_vec(end, 1) = 1 - sum(partial_pressure ./ p0_guess);
            
            % Jacobian Diagonal elements
            Jac(sub2ind([ncomp, ncomp], 1:ncomp-1, 1:ncomp-1)) = loadings(1:end-1)./p0_guess(1:end-1);            
            
            % Jacobian Last column
            Jac(1:ncomp-1, ncomp) = -1.0 * loadings(ncomp) / p0_guess(ncomp);
            
            % Jacobian Last row
            Jac(ncomp, :) = 1.0 .* partial_pressure ./ p0_guess.^2;
    
            % Update the Jac(ncomp, ncomp) entry based on modified fastIAST
            temp = (Jac(ncomp, 1:ncomp-1)...
                    ./...
                    Jac(sub2ind([ncomp, ncomp], 1:ncomp-1, 1:ncomp-1)));
            temp = temp(:);

            Jac(ncomp, ncomp) = Jac(ncomp, ncomp)...
                                - sum(temp .* Jac(1:ncomp-1, ncomp));

            Jac(ncomp, 1:ncomp-1) = 0;

            G_vec(ncomp) = G_vec(ncomp) - sum(temp .* G_vec(1:ncomp-1));
            
            % if any(isinf(temp))
            %     ss = 1;
            % end
            % 
            % if any(p0_guess<=0)
            %     ss = 1;
            % end

            % Calculate delta
            delta(ncomp) = G_vec(ncomp) / Jac(ncomp, ncomp);

            Jac_diag = Jac(sub2ind([ncomp, ncomp], 1:ncomp-1, 1:ncomp-1));
            Jac_diag = Jac_diag(:);

            delta(ncomp-1:-1:1) = (G_vec(ncomp-1:-1:1) - delta(ncomp) .* Jac(ncomp-1:-1:1, ncomp))...
                                    ./ Jac_diag(ncomp-1:-1:1);
            
            % delta(:) = Jac \ G_vec;
            
            % Update guess
            idx_check = (p0_guess - delta) < 0;       
            scaling = 1.0;
            update_guess(~idx_check) = p0_guess(~idx_check) - scaling.*delta(~idx_check);
            update_guess(idx_check) = 1/2 .* p0_guess(idx_check);

            % Convergence criteria
            mol_frac_conv_criteria = min(...
                                        max(1/abs(min(p0_guess)), 1e-10),...
                                        1e-8);

            mol_frac_check = ismembertol(sum(abs(partial_pressure)./abs(update_guess)), 1, mol_frac_conv_criteria);
            grand_potentials(:) = Isotherm_functions(ncomp, isotherm_params_array, update_guess', T, T_flag, 0);
            % delta_p0_ratio = sum(abs(delta)./abs(update_guess));
            
            if (std(grand_potentials) < convergence_tol) && (mol_frac_check || iter_count > 60)
                convergence_flag = 1;
                break;
            % if delta_p0_ratio < convergence_tol && mol_frac_check
            %     convergence_flag = 1;
            %     break;
            % elseif (sum(abs(delta)) < convergence_tol) && mol_frac_check
            %     convergence_flag = 1;
            %     break;
            else
                p0_guess(:) = update_guess;
                iter_count = iter_count + 1;
            end  
        end

        % Return solution
        if convergence_flag
            % Reshape and return the solution
            converged_sol = reshape(update_guess, 1, []);
        else
            error("MAX iteration reached NR-Method failed to converge");
        end
    end
end