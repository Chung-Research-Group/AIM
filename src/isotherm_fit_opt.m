function opt = isotherm_fit_opt(isotherm_model, loading_data, Pressure)
    %% Function variables
    num_param_vector = [2, 4, 3, 6, 3, 3, 3, 3, 3, 6, 3, 4, 6];
    UL_COMMON = 1e03;           % it was 1e03 initially, changed after ding STA model
    UL_LANG_CONSTANT = 1e03;    % Changed from 1 after trying high affinity xylene cases
    LL_theta = -100;        % Specific to Temkin model
    UL_theta = 100;         % Specific to Temkin model

    SAT_LOADING_GUESS = max(loading_data, [], "all");
    LANG_CONST_GUESS = 1e-10;
    EXPONENT_GUESS = 1.0;
    THETA_GUESS = 0.0;
    BET_TOL = 1e-8;         % Need tolerance to avoid getting undefined values in BET  

    % Get guesses for saturation loading and langmuir constant
    if (nargin==3) && length(loading_data)>1 
        [Q_SAT, K_LANG] = get_guess(loading_data, Pressure);
        if ~isnan(Q_SAT)
            SAT_LOADING_GUESS = Q_SAT;
        end
        if ~isnan(K_LANG)
            LANG_CONST_GUESS = K_LANG;
        end
    end

    % Get inflection pressure applicable only in STA model
    if strcmpi(isotherm_model, "Structural-Transition-Adsorption")    
        if nargin<3
            P_tr_guess = 1.0;    
        else
            try
                % Find inflection point pressure
                P_tr_guess = get_P_inflection(Pressure, loading_data);
            catch
                P_tr_guess = mean(Pressure);
            end
    
            if isnan(any(P_tr_guess)) || isempty(P_tr_guess)
                P_tr_guess = mean(Pressure);
            end
        end
    end

    %% Isotherm options structure based on the chosen isotherm model
    switch(isotherm_model)
        case 'SS-Langmuir'
            opt.fun = @SS_Langmuir;
            opt.num_params = num_param_vector(1); 
            % Lower Bounds
            opt.lb = zeros(1, opt.num_params);
            % Upper Bounds
            opt.ub = UL_COMMON .* ones(1, opt.num_params);
            opt.ub(2) = UL_LANG_CONSTANT;       % Upper bound for langmuir constant b 
            opt.guess = [SAT_LOADING_GUESS, LANG_CONST_GUESS];
            opt.T_flag = 1;
            opt.p_sat = 0;

        case 'DS-Langmuir'
            opt.fun = @DS_Langmuir;
            opt.num_params = num_param_vector(2);
            % Lower Bounds
            opt.lb = zeros(1, opt.num_params); 
            % Upper Bounds
            opt.ub = UL_COMMON .* ones(1, opt.num_params);
            opt.ub([2, 4]) = UL_LANG_CONSTANT;  % Upper bound for langmuir constant b and d
            opt.guess = [0.90*SAT_LOADING_GUESS, LANG_CONST_GUESS,... 
                         0.10*SAT_LOADING_GUESS, 0.5*LANG_CONST_GUESS];
            opt.T_flag = 1;
            opt.p_sat = 0;

        case 'SS-Langmuir-Freundlich'
            opt.fun = @SS_Langmuir_Freundlich;
            opt.num_params = num_param_vector(3);
            % Lower Bounds
            opt.lb = zeros(1, opt.num_params);
            % Upper Bounds
            opt.ub = UL_COMMON .* ones(1, opt.num_params);
            opt.ub(2) = UL_LANG_CONSTANT;       % Upper bound for langmuir constant b
            opt.guess = [SAT_LOADING_GUESS, LANG_CONST_GUESS, EXPONENT_GUESS];
            opt.T_flag = 1;
            opt.p_sat = 0;

        case 'DS-Langmuir-Freundlich'
            opt.fun = @DS_Langmuir_Freundlich;
            opt.num_params = num_param_vector(4);
            % Lower Bounds
            opt.lb = zeros(1, opt.num_params);
            % Upper Bounds
            opt.ub = UL_COMMON .* ones(1, opt.num_params);
            opt.ub([2, 5]) = UL_LANG_CONSTANT;  % Upper bound for langmuir constant b
            opt.guess = [0.75*SAT_LOADING_GUESS, LANG_CONST_GUESS, EXPONENT_GUESS ... 
                         0.25*SAT_LOADING_GUESS, 0.5*LANG_CONST_GUESS, EXPONENT_GUESS];
            opt.T_flag = 1;
            opt.p_sat = 0;

        case 'Quadratic'
            opt.fun = @Quadratic;
            opt.num_params = num_param_vector(5);
            % Lower Bounds
            opt.lb = zeros(1, opt.num_params);
            % Upper Bounds
            opt.ub = UL_COMMON .* ones(1, opt.num_params);
            opt.ub([2, 3]) = UL_LANG_CONSTANT;
            opt.guess = [0.5*SAT_LOADING_GUESS, LANG_CONST_GUESS, LANG_CONST_GUESS^2];
            opt.T_flag = 0;
            opt.p_sat = 0;

        case 'Temkin'
            opt.fun = @Temkin;
            opt.num_params = num_param_vector(6);
            % Lower Bounds
            opt.lb = zeros(1, opt.num_params); 
            opt.lb(3) = LL_theta; % Because theta parameter can be negative
            % Upper Bounds
            opt.ub = UL_COMMON .* ones(1, opt.num_params);
            opt.ub(2) = UL_LANG_CONSTANT;
            opt.ub(3) = UL_theta;
            opt.guess = [SAT_LOADING_GUESS, LANG_CONST_GUESS, THETA_GUESS];
            opt.T_flag = 0;
            opt.p_sat = 0;

        case 'BET'
            opt.fun = @BET;
            opt.num_params = num_param_vector(7);
            opt.lb = zeros(1, opt.num_params);
            opt.ub = UL_COMMON .* ones(1, opt.num_params);
            opt.ub([2, 3]) = UL_LANG_CONSTANT .* BET_TOL;
            opt.guess = [SAT_LOADING_GUESS, LANG_CONST_GUESS, LANG_CONST_GUESS^2];
            % opt.guess = [sat_loading_guess, sat_loading_guess, langmuir_k_guess^2];
            opt.T_flag = 0;
            opt.p_sat = 0;

        case 'Sips'
            opt.fun = @Sips;
            opt.num_params = num_param_vector(8);
            opt.lb = zeros(1, opt.num_params);
            opt.ub = UL_COMMON .* ones(1, opt.num_params);
            opt.ub(2) = UL_LANG_CONSTANT;
            opt.guess = [SAT_LOADING_GUESS, LANG_CONST_GUESS, EXPONENT_GUESS];
            % opt.guess = [sat_loading_guess, sat_loading_guess, langmuir_k_guess^2];
            opt.T_flag = 0;
            opt.p_sat = 0;

         case 'Toth'
            opt.fun = @Toth;
            opt.num_params = num_param_vector(9);
            opt.lb = zeros(1, opt.num_params);
            opt.ub = UL_COMMON .* ones(1, opt.num_params);
            opt.ub(2) = UL_LANG_CONSTANT;
            opt.guess = [SAT_LOADING_GUESS, LANG_CONST_GUESS, EXPONENT_GUESS];
            % opt.guess = [sat_loading_guess, sat_loading_guess, langmuir_k_guess^2];
            opt.T_flag = 0;
            opt.p_sat = 0;

         case 'Structural-Transition-Adsorption'
            opt.fun = @STA;
            opt.num_params = num_param_vector(10);
            opt.lb = zeros(1, opt.num_params);
            opt.ub = UL_COMMON .* ones(1, opt.num_params);
            % opt.ub([2, 3]) = UL_langmuir_constant .* BET_tol;
            % opt.guess = [sat_loading_guess, K_guess, C_guess, n_guess];
            % The guess for q_LP should be higher than q_NP
            opt.guess = [SAT_LOADING_GUESS.*0.25, LANG_CONST_GUESS,...
                         SAT_LOADING_GUESS.*0.75, LANG_CONST_GUESS,...
                         EXPONENT_GUESS, P_tr_guess];
            opt.T_flag = 0;
            opt.p_sat = 0;
        
        case 'Dubinin-Astakhov'
            opt.fun = @DA_Isotherm;
            opt.num_params = num_param_vector(11);
            opt.lb = zeros(1, opt.num_params);
            opt.lb(2) = 1e-05;          % TO avoid inifnitny 
            opt.ub = UL_COMMON .* ones(1, opt.num_params);
            % opt.ub([2, 3]) = UL_langmuir_constant .* BET_tol;
            opt.guess = [SAT_LOADING_GUESS, 1, EXPONENT_GUESS];
            opt.T_flag = 0;
            opt.p_sat = 1;

        case 'Klotz'
            opt.fun = @Klotz;
            opt.num_params = num_param_vector(12);
            opt.lb = zeros(1, opt.num_params);
            opt.ub = UL_COMMON .* ones(1, opt.num_params);
            % opt.ub([2, 3]) = UL_langmuir_constant .* BET_tol;
            opt.guess = [SAT_LOADING_GUESS, 1, 1, EXPONENT_GUESS];
            opt.T_flag = 0;
            opt.p_sat = 1;

        case 'Do-Do'
            opt.fun = @DD_Isotherm;
            opt.num_params = num_param_vector(13);
            opt.lb = zeros(1, opt.num_params);
            opt.ub = UL_COMMON .* ones(1, opt.num_params);
            opt.ub(2) = 1.0;        % f parameter is a fraction ranging from o to 1 in DD model
            % opt.ub([2, 3]) = UL_langmuir_constant .* BET_tol;
            % opt.guess = [sat_loading_guess, K_guess, C_guess, n_guess];
            opt.guess = [SAT_LOADING_GUESS, 0.0, 0.0, 0.0,...
                         EXPONENT_GUESS, 2*EXPONENT_GUESS];
            opt.T_flag = 0;
            opt.p_sat = 1;
            
        case 'Virial'
            opt.fun = @Virial_func;
            opt.num_a = 7;      % Refers to number of "a" parameters in Virial equation
            opt.num_b = 3;      % Refers to number of "b" parameters in Virial equation
            opt.num_params = opt.num_a + opt.num_b;
            % opt.lb = zeros(1, opt.num_params);
            opt.lb = -Inf.* ones(1, opt.num_params);
            opt.ub = Inf .* ones(1, opt.num_params);
            % opt.ub([2, 3]) = UL_langmuir_constant .* BET_tol;
            opt.guess = zeros(1, opt.num_params);
            opt.T_flag = 0;
            opt.p_sat = 0;

        otherwise
            error('Invalid isotherm model: %s', isotherm_model);
    end

    opt.name = isotherm_model;

    %% Checking the initial guess values
    % if nargin > 2  % Initial guess specified by user
    %     %% Validation of initial guess
    %     % Checking if number of initial guesses/lower bounds/upper bounds is equal to number of isotherm parameter
    %     if size(initial_guess, 2) ~= opt.num_params
    %         warning("Number of initial guesses/lower bounds/upper bounds/ does not match the number of parameters for the chosen isotherm model. " + ...
    %                 "Using default initialization...")
    %         % opt.guess = default_initial_guess .* ones(1, opt.num_params);
    % 
    %     % Specifying the new bounds and checking if any initial guess value is not within bounds
    %     else
    %         % opt.guess = initial_guess;
    %         is_in_range = initial_guess >= opt.lb & initial_guess <= opt.ub;
    %         if ~all(is_in_range)
    %             warning("Some of the initial guess values provided are " + ...
    %                 "not within the bounds. Those values will be initialized using default initial guesses.");
    % 
    %             % index_out_range = find(~is_in_range);
    %             % opt.guess(index_out_range) = default_initial_guess .* ones(1, length(index_out_range));
    %         end
    % 
    %         index_in_range = find(is_in_range);
    %         opt.guess(index_in_range) = initial_guess(is_in_range);
    %     end      
    % end
    %
    %% Miscelleneous Functions
    % SS-Langmuir
    function loading = SS_Langmuir(params, P)
        loading = (params(1)*params(2).*P) ./ (1+params(2).*P);
    end
    
    % DS-Langmuir
    function loading = DS_Langmuir(params, P)
        loading = (params(1)*params(2).*P) ./ (1+params(2).*P)... 
                + (params(3)*params(4).*P) ./ (1+params(4).*P);
    end

    % SS-Langmuir-Freundlich
    function loading = SS_Langmuir_Freundlich(params, P)
        loading = (params(1)*params(2).*P.^params(3)) ./ (1+params(2).*P.^params(3));
    end

    % DS-Langmuir-Freundlich
    function loading = DS_Langmuir_Freundlich(params, P)
        loading = (params(1)*params(2).*P.^params(3)) ./ (1+params(2).*P.^params(3))...
                + (params(4)*params(5).*P.^params(6)) ./ (1+params(5).*P.^params(6));
    end

    % Quadratic Isotherm
    function loading = Quadratic(params, P)
        loading = params(1) .* (params(2).*P + 2*params(3).*P.^2) ...
                ./ (1 + params(2).*P + params(3).*P.^2);
    end

    % Temkin Isotherm
    function loading = Temkin(params, P)
        m = params(1);
        b = params(2);
        theta = params(3);

        lang_term = b.*P ./ (1+b.*P); 
        first_term = m .* lang_term;
        first_multiplier = m.*theta .* (lang_term).^2;
        second_multiplier = lang_term - 1;
        loading = first_term + first_multiplier .* second_multiplier;
    end
    
    % BET Isotherm
    function loading = BET(params, P)
        m = params(1);
        b_surface = params(2);
        b_layers = params(3);

        loading = m.*b_surface.*P ./ (1 - b_layers.*P) ./ (1 - b_layers.*P + b_surface.*P);
    end

    % Sips Isotherm
    function loading = Sips(params, P)
        m_sat = params(1);
        b = params(2);
        n = params(3);
        num = m_sat .* (b.*P).^(1/n);
        den = 1 + (b.*P).^(1/n);
        loading = num ./ den;
    end
    
    % Toth Isotherm
    function loading = Toth(params, P)
        m_sat = params(1);
        b = params(2);
        n = params(3);
        num = m_sat .* b.*P;
        den = (1 + (b.*P).^(n)).^(1/n);
        loading = num ./ den;
    end

    % Klotz Isotherm
    function loading = Klotz(params, P)
        % Here P = relative pressure
        m_sat = params(1);
        K = params(2);
        C = params(3);

        % n can only have integer value
        % n = round(params(4));
        n = params(4);

        q = K.*P;
    
        num = m_sat.*C.*q .* (1 - (1+n).*q.^n + n.*q.^(n+1));
        den = (1-q) .* (1 + (C-1).*q - C.*q.^(n+1));        
        loading = num ./ den;
    end
    
    % DD-Isotherm
    function loading = DD_Isotherm(params, P)
        % Here P = relative pressure
        m_sat = params(1);
        f = params(2);
        K_f = params(3);
        K_mu = params(4);

        % Alpha and beta can only have integer value
        % alpha = round(params(5));
        % beta = round(params(6));
        alpha = params(5);
        beta = params(6);
        

        % Based on DD model, beta is always greater than alpha
        if beta<alpha
            loading = zeros(size(P));
        else        
            num = f .* (K_f.*P) .* (1 - (1+beta).*P.^beta + beta.*P.^(beta+1));
            den = (1-P) .* (1 + (K_f-1).*P - K_f.*P.^(beta+1));
    
            first_term = num ./ den;
            second_term = (1-f) .* (K_mu.*P.^alpha) ./ (1 + K_mu.*P.^alpha);
    
            loading = m_sat.*(first_term + second_term);
        end
    end

    % DA isotherm
    function loading = DA_Isotherm(params, P)
        % Here P = relative pressure
        m = params(1);
        epsilon = params(2);
        n = params(3);

        loading = m .* exp(-(1/epsilon .* log(1./P)).^n);
    end

    % STA Isotherm
    function loading = STA(params, P)
        q_NP = params(1);
        b_NP = params(2);
        q_LP = params(3);
        b_LP = params(4);
        s = params(5);
        P_tr = params(6);
        
        den_NP = (1+b_NP.*P);
        den_LP = (1+b_LP.*P);

        % y = ((1+b_NP*P_tr)./(den_NP)).^q_NP .* ((1+b_LP*P_tr)./(den_LP)).^q_LP;
        y = ((1+b_NP*P_tr)./(den_NP)).^q_NP .* ((1+b_LP*P)./(1+b_LP.*P_tr)).^q_LP;

        sigma = y.^s ./ (1+y.^s);
        
        lang_NP = q_NP.*b_NP.*P./den_NP;
        lang_LP = q_LP.*b_LP.*P./den_LP;

        loading = (1-sigma).*lang_NP + sigma.*lang_LP;
    end
    
    % Virial polynomial function
    function pressure = Virial_func(params, loading, temperature, num_a_params)
        % returns log of pressure values
        % Flip function is required because the polyval accepts the
        % coefficients in the order of decreasing degree.
        a_terms = 1./temperature .* polyval(flip(params(1:num_a_params)), loading);
        b_terms = polyval(flip(params(num_a_params+1:end)), loading);
        pressure = a_terms + b_terms + log(loading);
    end
    %
    %% Function to calculate initial guesses for saturation loading and Langmuir constant
    function [sat_loading_guess, langmuir_constant_guess] = get_guess(loading_data, pressure_data)
        try 
            % Check whether loading data is of single isotherm or multiple
            % isotherm
            [~, col] = size(loading_data);
            if col > 1
                [~, idx] = max(loading_data, [], "all");
                [~, max_val_col] = ind2sub(size(loading_data), idx);

                % Get loading column corresponding to max value
                loading_data = loading_data(:, max_val_col);
                pressure_data = pressure_data(:, max_val_col);    
            end
            
            % Remvoing zero-values entries
            [idx_ret, ~] = find(loading_data~=0.0);
            loading_data = loading_data(idx_ret);
            pressure_data = pressure_data(idx_ret);
            
            data_matrix = sortrows(horzcat(loading_data, pressure_data), 1, "ascend");
            loading_data = data_matrix(:, 1);
            pressure_data = data_matrix(:, 2);
            [~, idx_max] = max(loading_data);
             
            sat_loading_guess = mean(loading_data(idx_max-4:idx_max));
            % K_H = (loading_data(4)-loading_data(1))./ (pressure_data(4)-pressure_data(1));
            K_H = loading_data(1:4)'/pressure_data(1:4)';
            langmuir_constant_guess = K_H / sat_loading_guess;
        catch
            sat_loading_guess = NaN;
            langmuir_constant_guess = NaN;
        end
    end
    %% Function to find the inflection point for STA isotherm
    function P_inflection = get_P_inflection(Pressure, loading)
        % smoothing the Pressure values and using spline interpolation
        size_p = length(Pressure);
        resolve_factor = 2;
        pressure_smoothed = linspace(min(Pressure), max(Pressure), resolve_factor*size_p);
        
        % Smoothing the loading using 'loess' model
        loading = smooth(loading, 'loess');
        loading_smoothed = spline(Pressure, loading, pressure_smoothed);
        
        % Derivative Calculation using Finite Differences
        dqdP = diff(loading_smoothed) ./ diff(pressure_smoothed);
        d2qdp2 = diff(dqdP) ./ diff(pressure_smoothed(2:end));
        
        % Idx where decond derivative changes sign, in case of multiple
        % inflection points, only the first one is selected
        inflection_idx = find(diff(sign(d2qdp2)), 1, "first");
        
        P_inflection = pressure_smoothed(inflection_idx+2); % +2 because of finite differences
    end
end