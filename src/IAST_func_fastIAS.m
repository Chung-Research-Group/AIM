function partial_loadings = IAST_func_fastIAS(num_components, isotherm_params_array, Pressure, gas_phase_mol_fraction, T_array, T_flag, initial_guess)
global cached_p0; 
%% Parameter Unpacking and Decalarations
    size_of_pressure_vector = length(Pressure);     % Length of mole fraction vector
    partial_pressures = Pressure .* gas_phase_mol_fraction;
    adsorbed_mole_fractions = zeros(size_of_pressure_vector, num_components);
    fictitious_pressure = zeros(size_of_pressure_vector, num_components);
    %   
    %% Initialization of Guess for IAST Equations Solution
    % Solve the non-linear equations returened by the 'func' functions using fsolve method.
    % If Initial Guess has not been provided in input arguments, 
    % it will be guessed proportional to pure compoenent loading.
    
    if nargin<7 || ~(any(sum(initial_guess, 2)))
        loading_array = Isotherm_functions(num_components, isotherm_params_array, partial_pressures, T_array, T_flag, 1);
        loading_array(loading_array == 0) = 1e-07;
        initial_guess = loading_array ./ (sum(loading_array, 2));
        init_guess_p0 = partial_pressures ./ initial_guess;
    end
    
    % Use to avoid undefined values when the column is filled 
    % with only non-adsorbing gas or when the gas loading is very low
    tolerance = 1e-10;
    
    % guess = initial_guess(:, 1:num_components-1);
    useCache = 1;
    if ~(isempty(cached_p0)) && useCache
        guess = cached_p0;
    else
        guess = init_guess_p0(:, 1:num_components);
        cached_p0 = zeros(size(guess));
    end
    
%     options = optimset('Display','off', 'TolFun', 1e-10, 'MaxIter', 400, 'JacobPattern', jp);   % Options struc for fsolve
    % options = optimset('Display','on', 'TolFun', 1e-10, 'MaxIter', 10e5, 'OptimalityTolerance', 1e-10);   % Options struc for fsolve
    % options = optimset('TolFun', 1e-8, 'Display', 'on', 'Diagnostics', 'off', 'TolX', 1e-09, 'FinDiffRelStep', 1e-9, 'Algorithm','levenberg-marquardt',...
    %     'MaxIter', 10e5);
    options = optimoptions("fsolve", "Display","off", "FunctionTolerance",1e-10, "MaxIterations",10e5, "OptimalityTolerance",1e-12, "StepTolerance",1e-12, ...
                            "FiniteDifferenceType","central", "FiniteDifferenceStepSize",1e-9, "MaxFunctionEvaluations",10e5);
    
    for k = 1:size_of_pressure_vector
        try
            % In case of noadsorbing component
            if sum(partial_pressures(k, :), 2) <= tolerance
                % adsorbed_mole_fractions(k, :) = zeros(1, num_components)+tolerance;
                fictitious_pressure(k, :) = zeros(1, num_components);
                cached_p0(k, :) = fictitious_pressure(k, :);
                % cached_p0(k, :) = partial_pressures(k, :)./adsorbed_mole_fractions(k, :);
                
            % In case of single adsorbing component, the initial guess will give the correct loading
            elseif num_components == 1
                % adsorbed_mole_fractions(k, :) = initial_guess(k, :);
                fictitious_pressure(k, :) = init_guess_p0(k, :);

            % Calculate the IAST loadings 
            else               
                function_handle = @(x)residual_func(x, isotherm_params_array, partial_pressures(k, :), T_array(k, 1));
                
                % Use default set of options and cache initial guess
                [solution, ~, exitflag] = fsolve(function_handle, guess(k, :), options);
                
                % Use default set of options and runtime loading based calculated initial guess
                if exitflag <= 0
                    [solution, ~, exitflag] = fsolve(function_handle, init_guess_p0(k, :), options);
                end
                
                % Change algorithm and use cached initial guess
                if exitflag <= 0
                    options.Algorithm = 'levenberg-marquardt';
                    % options.Display = 'final';
                    [solution, ~, exitflag] = fsolve(function_handle, guess(k, :), options);
                end
                
                % Change algorithm and use runtime loading based calculated initial guess
                if exitflag <= 0
                    options.Algorithm = 'levenberg-marquardt';
                    % options.Display = 'final';
                    [solution, ~, exitflag] = fsolve(function_handle, init_guess_p0(k, :), options);                 
                end
                
                % Check exit flags
                if exitflag == 0
                    error("Maximum Iterations for fsolve reached! Terminating the solution.");
                elseif exitflag < 0
                    error("fsolve failed to solve the equations!")
                end
                
                % Cache the solution
                cached_p0(k, :) = solution;
                
                % % calculate adsorbed mole fractions
                % adsorbed_mole_fractions(k, :) = partial_pressures(k, :)./solution;
                fictitious_pressure(k, :) = solution;
                % 
                % % Tolerance value check that sum of adsorbed mole fraction is unity
                % ads_mol_frac_unity_tol = 1e-5;    
                % if ~(ismembertol(sum(adsorbed_mole_fractions(k, :)), 1, ads_mol_frac_unity_tol))
                %     error("The sum of adsorbed mole fractions is not unity, fsolve failed to solve the equations!")
                % end
            end
        catch ME1
                error(['Failed to solve the IAST equations...!', ...
                       ' Error:%s'], ME1.message);    
        end
    end
    
    % adsorbed_mole_fractions(adsorbed_mole_fractions == 0) = tolerance;
    % fictitious_pressure(fictitious_pressure == 0) = tolerance;

    adsorbed_mole_fractions = partial_pressures ./ fictitious_pressure;
    % adsorbed_mole_fractions(adsorbed_mole_fractions == 0) = tolerance;
    adsorbed_mole_fractions(isnan(adsorbed_mole_fractions)) = 0;
    adsorbed_mole_fractions(isinf(adsorbed_mole_fractions)) = 0;

    % Tolerance value check that sum of adsorbed mole fraction is unity
    ads_mol_frac_unity_tol = 1e-5;    
    if ~all(ismembertol(sum(adsorbed_mole_fractions, 2), 1, ads_mol_frac_unity_tol))
        disp('hassan');
        % error("The sum of adsorbed mole fractions is not unity, fsolve failed to solve the equations!")
    end

    % pressure_0 = partial_pressures ./ adsorbed_mole_fractions;
    loading_array = Isotherm_functions(num_components, isotherm_params_array, fictitious_pressure, T_array, T_flag, 1);
    % loading_array(loading_array==0) = tolerance;
    
    inverse_loading = sum(adsorbed_mole_fractions./loading_array, 2);
    % inverse_loading(inverse_loading==0) = tolerance;
    loading_total = 1 ./ inverse_loading;
    loading_total(isnan(loading_total)) = tolerance;
    loading_total(isinf(loading_total)) = tolerance;

    partial_loadings = (adsorbed_mole_fractions .* loading_total);
    
    if any(find(isnan(partial_loadings))) 
        s=5;
    elseif any(find(isinf(partial_loadings)))
        s=8;
    end
    function my_res = residual_func(p0, isotherm_params_array, partial_pressure, T)
        % Calculate the difference in spreading pressure for (N-1) components.
        % The spreading pressures are all calculated at the fictious pressure calculated based on the
        % adsorbed mole fraction and partial pressure of the components.
        % Fictitous Pressure = p_* = p_i / x_i
        
        s = 1e-08;       
        % x = [(adsorbed_MF), (1-sum(adsorbed_MF))];
        % x(x==0) = s;
        
        if any(p0 < 0)
            my_res = 1e05.*ones(1, size(p0, 2));
            return
        end
        
        % To avoid undefined values
        p0(p0==0) = s;

        spreading_pressures = Isotherm_functions(num_components, isotherm_params_array, p0, T, T_flag, 0);
        
        sp_pressure_end = spreading_pressures(1, end);

        spreading_pressures_diff = spreading_pressures(:, 1:num_components-1) - sp_pressure_end;
        close_cond = 1 - sum(partial_pressure ./ p0);

        if ~isreal(spreading_pressures_diff)
            my_res = 1e05.*ones(1, size(p0, 2));
        else
            my_res = [spreading_pressures_diff, close_cond];
        end
        if any(isnan(my_res))
            error("NaN values in IAST: residual function...")
        end
    end
end