function solution = Adsorption_Breakthrough_case(parameter_set, progress_bar)
    %% Global Vars
    global cached_p0;
    global Henry_Coeff;
    global cached_p0_local;
    cached_p0 = [];

    %% Prameter loading and unpacking
    num_component = parameter_set.CompNum;
    isotherm_params = parameter_set.IsothermParams;
    isotherm_params(isnan(isotherm_params)) = 0;        % All of the unspecified/NaN values are converted to zero 
    mixture_predict_method = parameter_set.MixtureMode;
    calc_type = parameter_set.CalcMode;
    numerical_params = parameter_set.NumPar;

    % Calling input_params_parsing function
    params = ProcessInputParameters(parameter_set);
    N   = params(1);
    L =   params(2);    
    T_0 = params(4)         ;
    mu   = params(5)        ;
    R = params(6)           ;
    v_0 = params(7)         ;
    q_s0 = params(8)        ;
    P_0 = params(14)        ;
    dia_in = params(16)     ;
    MW_1 = params(19)       ;
    MW_2 = params(20)       ;
    MW_3 = params(21)       ;
    MW_4 = params(22)       ;
    MW_5 = params(23)       ;
    y_0_1 = params(24)      ;
    y_0_2 = params(25)      ;
    y_0_3 = params(26)      ;
    y_0_4 = params(27)      ;
    epsilon = params(28)    ;
    r_p = params(31)        ;
    ro_ads = params(32)     ;
    time =   params(36) * L/v_0;
    n_dot0 = params(38);
    y_init_1 = params(40);
    y_init_2 = params(41);
    y_init_3 = params(42);
    y_init_4 = params(43);
    T_amb    = params(59);
%   timespan = (0:0.01:time) .* v_0/L;
    timespan = (0:1:time+0.1) .* v_0/L;
    
    AbsTol = numerical_params(3);
    RelTol = numerical_params(4);
    
    % Calculating Henry coefficient
    p_range = linspace(1e-7, 1e-6, 2).';
    p_comp = repmat(p_range, 1, num_component-1);
    loading = Isotherm_functions(num_component-1, isotherm_params, p_comp, 298.15, 0, 1);
    Henry_Coeff = (loading(end, :) - loading(end-1, :)) ./ (p_range(end)-p_range(end-1));

    %% Initial Conditions
    x0 = zeros(11*N+22, 1);
    
    % Pressure
    x0(1) = 1;
    x0(2:N+2) = 1;
    
    % 1st Component
    x0(N+3)  = y_0_1;
    x0(N+4:2*N+4) = y_init_1;
    
    % 2nd Component
    x0(2*N+5) = y_0_2;
    x0(2*N+6:3*N+6) =y_init_2;
 
    % 3rd Component
    x0(3*N+7) = y_0_3;
    x0(3*N+8:4*N+8) = y_init_3;
    
    % 4th Component
    x0(4*N+9) = y_0_4;
    x0(4*N+10:5*N+10) = y_init_4;
    
    % Molar loadings
    x0(5*N+11:6*N+12) = 0;
    x0(6*N+13:7*N+14) = 0;
    x0(7*N+15:8*N+16) = 0;
    x0(8*N+17:9*N+18) = 0;
    x0(9*N+19:10*N+20) = 0;
    
    % Gas Temperature
    x0(10*N+21:11*N+22) = T_0/T_0;

    % Wall Temperature
    x0(11*N+23:12*N+24) = T_amb/T_0;
 %
    %% Function Alias for isothermal/nonisothermal calculations
    switch calc_type
        case 0
            isotherm_params(end-1, :) = 0;   % dH=0
            isotherm_params(end, :) = T_0;   % T = T_feed
            Adsorption_fxn_IAST     = @(t, x) Adsorption_step_IAST_Iso(t, x, params, isotherm_params, num_component)          ;
            Adsorption_fxn_Ext_Lang = @(t, x) Adsorption_step_EDSL_Iso(t, x, params, isotherm_params, num_component)          ;
            x0 = x0(1:10*N+20, 1);    % Not solving for Gas and Wall Temperature values
        case 1
            T_ref = isotherm_params(end, 1:num_component-1);
            isotherm_params(end-1, :) = isotherm_params(end-1, :).*1e03; % dH conversion from kJ/mole to J/mole
            if ~all(T_ref)
                error(['Reference Temperatures can not be zero or unspecified for non-isothermal calculations.' ...
                        'Please specify valid reference temperatures for all components.']);
            end
            Adsorption_fxn_IAST     = @(t, x) Adsorption_step_IAST_Non_Iso(t, x, params, isotherm_params, num_component)          ;
            Adsorption_fxn_Ext_Lang = @(t, x) Adsorption_step_EDSL_Non_Iso(t, x, params, isotherm_params, num_component)      ;

        otherwise
            error('Invalid Calculation type...!');
    end

    event_fxn = @(t, x) event_function(t, x) ;
    %
    %% ODE numerical options and calling function
    Jac_patt = Jacobian(N, calc_type);
    % options = odeset('JPattern',Jac_patt, 'RelTol', RelTol, 'AbsTol', AbsTol, 'Events', event_fxn, 'MaxStep', 5);
    options = odeset('JPattern',Jac_patt, 'RelTol', RelTol, 'AbsTol', AbsTol, 'Events', event_fxn);
     % options = odeset('RelTol', RelTol, 'AbsTol', AbsTol, 'Events', event_fxn);

    if mixture_predict_method == 0
        if any(isotherm_params(1, 1:num_component) > 2)
            error(['In case of Extended Langmuir method, isotherm type of' ...
                ' all components must be either SS-Langmuir or DS-Langmuir!'])
        else
            [t_1, ode_sol] = ode15s(Adsorption_fxn_Ext_Lang, timespan, x0, options);
            % sol_struc = ode15s(Adsorption_fxn_Ext_Lang, [0 time].*v_0/L, x, options);
        end

    elseif mixture_predict_method == 1
        % tic;
        [t_1, ode_sol] = ode15s(Adsorption_fxn_IAST, timespan, x0, options);
        % sol_struc = ode15s(Adsorption_fxn_IAST, [0 time].*v_0/L, x, options);
        % toc;
    else
        error("Please specify the mixture isotherm prediction type!")
    end
%
    %% Solution Processing
    % Interpolating solution variables
    % t_1 = timespan;
    % ads_reac_sol = deval(sol_struc, t_1);
    % 
    % % Transporsing solution variables such that column ==> nodes, rows ==>
    % % time values
    % ads_reac_sol = ads_reac_sol';

    % Solution Correction
    ode_sol = Correction_func(ode_sol);

    % Pressure Correction
    ode_sol = pressure_correction(ode_sol, n_dot0, "Bottom");

    % Unpacking solution in respective variables
    % solution.t = t_1'.*L/v_0;
    solution.t = t_1.*L/v_0;
    % solution.P    = ads_reac_sol(:, 1: N+2) .* P_0/1000;     % Pressure in kPa
    solution.P    = ode_sol(:, 1: N+2) .* P_0;            % Pressure in Pa
    solution.C1 = max(ode_sol(:, N+3: 2*N+4),  0);          
    solution.C2  = max(ode_sol(:, 2*N+5: 3*N+6), 0);
    solution.C3   = max(ode_sol(:, 3*N+7: 4*N+8), 0);
    solution.C4  = max(ode_sol(:, 4*N+9: 5*N+10), 0);
    solution.C5 = max(1 - solution.C1 - solution.C2 - solution.C3 - solution.C4, 0);
    solution.x1C1 = ode_sol(:, 5*N+11:6*N+12) .* q_s0  ;
    solution.x2C2 = ode_sol(:, 6*N+13:7*N+14) .* q_s0  ;
    solution.x3C3 = ode_sol(:, 7*N+15:8*N+16) .* q_s0  ;
    solution.x4C4 = ode_sol(:, 8*N+17:9*N+18) .* q_s0  ;
    solution.x5C5 = ode_sol(:, 9*N+19:10*N+20).* q_s0  ;
    
    % Temperature in K
    if calc_type
        solution.T    = (ode_sol(:, 10*N+21:11*N+22) .* T_0);
        solution.Tw   = (ode_sol(:, 11*N+23:12*N+24) .* T_0);
    else
        solution.T = ones(size(ode_sol(:, 1:N+2))).*T_0;
        solution.Tw = ones(size(ode_sol(:, 1:N+2))).*T_0;
    end

    % g = adsorbed_amount(solution);
%%  Supplementary Functions
    %% Solution Correction Function
    function x_new = Correction_func(x)
    % This function checks and correct the boundary values of solution array
    % and also checks the bounds of mole fraction.
    % INPUT:
    % x  = Solution matrix from ode solver.
    
    % OUTPUT:
    % x_new = Corrected Solution matrix
    
        idx             = find(x(:, N+1) < 1)               ;  %  if P_N+1 < 1 then P_N+2 = P_N+1
        x(idx, N+2)     = x(idx, N+1)                       ;  % P_N+2 = P_N+1
        
        x(:, N+3:2*N+4) = max(min(x(:, N+3:2*N+4), 1), 0)       ;  % 0 <= y => 1
        x(:, 2*N+5:3*N+6) = max(min(x(:, 2*N+5:3*N+6), 1), 0)   ;  % 0 <= y => 1
        x(:, 3*N+7:4*N+8) = max(min(x(:, 3*N+7:4*N+8), 1), 0)   ;  % 0 <= y => 1
        x(:, 4*N+9:5*N+10) = max(min(x(:, 4*N+9:5*N+10), 1), 0) ;  % 0 <= y => 1
        
        % Mole fraction at outlet boundary
        x(:, 2*N+4) = x(:, 2*N+3);
        x(:, 3*N+6) = x(:, 3*N+5);
        x(:, 4*N+8) = x(:, 4*N+7);
        x(:, 5*N+10) = x(:, 5*N+9);
        
        % Molar Loading at inlet boundary
        x(:, 5*N+11)     = x(:, 5*N+12)                     ;  % x1_1 = x1_2
        x(:, 6*N+13)     = x(:, 6*N+14)                     ;  % x2_1 = x2_2
        x(:, 7*N+15)     = x(:, 7*N+16)                     ;  % x3_1 = x3_2
        x(:, 8*N+17)     = x(:, 8*N+18)                     ;  % x4_1 = x4_2
        x(:, 9*N+19)     = x(:, 9*N+20)                     ;  % x5_1 = x5_2
    
        % Molar Loading at outlet boundary
        x(:, 6*N+12)     = x(:, 6*N+11)                     ;  % x1_N+2 = x1_N+1
        x(:, 7*N+14)     = x(:, 7*N+13)                     ;  % x2_N+2 = x2_N+1
        x(:, 8*N+16)     = x(:, 8*N+15)                     ;  % x3_N+2 = x3_N+1
        x(:, 9*N+18)     = x(:, 9*N+17)                     ;  % x4_N+2 = x4_N+1
        x(:, 10*N+20)    = x(:, 10*N+19)                    ;  % x5_N+2 = x5_N+1
        
        if calc_type
            % Temperature at outlet boundary
            x(:, 11*N+22) = x(:, 11*N+21);
        end
        
        x_new = x   ;    
    end
%
    %% Pressure Correction Function
    function corrected_solution = pressure_correction(solution, molarflux, column_end)
    % This function is to be used when the coulmn inlet is specified as
    % velocity inlet. The function uses Ergun's equation to reverse
    % calculate the Pressure at inlet boundary.
    % INPUT:
    % solution= Solution array from ODE step calculation
    % molarflux = Molar flux into the column
    % column_end = The end of the column where the material enters the
    % column.
    %
    % OUTPUT:
    % corrected_solution = The solution array with corrrected pressure at
    % the boundary.
    % In our cycle the material is introduced always from the bottom,
    % hence, the default column_end is bottom.

        if nargin<3
        column_end = 'Bottom';
        end
        
        % Non-Dimensional step length
        dz  = 1/N ;

        % Unpacking data
        if strcmpi(column_end, "Top") == 1
            P = solution(:, N+1) .* P_0 ;
            if calc_type
                T = solution(:, 11*N+22) .* T_0  ;
            else
                T = ones(size(P)).*T_0;
            end
            y_1 = solution(:, 2*N+4);
            y_2 = solution(:, 3*N+6);
            y_3 = solution(:, 4*N+8);
            y_4 = solution(:, 5*N+10);
            y_5  = 1 - y_1 - y_2 - y_3- y_4 ;
        end
        
        if strcmpi(column_end, "Bottom") == 1
            P = solution(:, 2) .* P_0;
            if calc_type
                T = solution(:, 10*N+21) .* T_0;
            else
                T = ones(size(P)).*T_0;
            end
            y_1 = solution(:, N+3);
            y_2 = solution(:, 2*N+5);
            y_3 = solution(:, 3*N+7);
            y_4 = solution(:, 4*N+9);
            y_5  = 1 - y_1 - y_2 - y_3- y_4 ;
        end
        
        % Calculate the Molar Mass of the gas at 1st node [kg/m^3]
        MW_t = y_1 .* MW_1 + y_2 .* MW_2 + ...
                y_3 .* MW_3 + y_4 .* MW_4 + y_5 .* MW_5;
       
        a_1   = 150*mu*(1-epsilon)^2*dz*L/2/4/r_p^2/epsilon^3./T./R ;
        a_2_1 = 1.75*(1-epsilon)/2/r_p/epsilon/epsilon/epsilon*dz*L/2;
        a_2   = a_2_1/R./T.*molarflux.*MW_t;
    
        a =  a_1+a_2;
        b =  P./T./R;
        c = -molarflux;
    
        vh = (-b+sqrt(b.^2-4.*a.*c))/2./a./v_0 ;

        a_p = a_1.*T*R      ;
        b_p = a_2_1.*MW_t./R./T ;
        
        if strcmpi(column_end, "Top") == 1
            solution(:, N+2)  = ((a_p.*vh.*v_0+P)./(1-b_p.*vh.*v_0*vh.*v_0))./P_0;
        end

        if strcmpi(column_end, "Bottom") == 1
            solution(:, 1)  = ((a_p.*vh.*v_0+P)./(1-b_p.*vh.*v_0.*vh.*v_0))./P_0;
        end    
        corrected_solution = solution;  
    end
%
    %% Breakthrough Area calculation
    function [n_ads_tot, n_loading] = adsorbed_amount(solution)
        %% Function variables
        R = 8.314;
        area_cross_sect = pi()*(dia_in/2)^2;
        dz = L/N;
        t = solution.t;
        V_tot_col = area_cross_sect .* L;
        mass_ads = ro_ads * V_tot_col;
        %
        %% Inlet
        P_in = solution.P(:, 1).*1e03;
        T_in = solution.T(:, 1)+273.15;
        y_in = [solution.C1(:, 1), solution.C2(:, 1), solution.C3(:, 1),...
                solution.C4(:, 1), solution.C5(:, 1)];
        C_tot_0 = P_in./T_in./R;
        ndot_tot_0 = v_0.*C_tot_0.*area_cross_sect;
        %
        %% Outlet
        P_out = solution.P(:, end-1:end).*1e03;
        T_out = solution.T(:, end)+273.15;
        y_out = [solution.C1(:, end), solution.C2(:, end), solution.C3(:, end),...
                 solution.C4(:, end), solution.C5(:, end)];        
        C_tot_out = P_out(:, 2)./T_out./R;
        
        % Velocity calculation
        ro_g = sum([MW_1, MW_2, MW_3, MW_4, MW_5].*y_out, 2) .* C_tot_out;
        dPdz = 2*(P_out(:, 2)-P_out(:, 1))./dz ;

        viscous_term =  150*mu*(1-epsilon)^2/4/r_p^2/epsilon^3             ;
        kinetic_term = (1.75*(1-epsilon)/2/r_p/epsilon^3) * ro_g           ;
        v_out            = -sign(dPdz).*(-viscous_term+(abs(viscous_term^2+...  
                        4*kinetic_term.*abs(dPdz))).^(.5))/2./kinetic_term ;
        
        ndot_tot_out = abs(v_out).*C_tot_out.*area_cross_sect;
        %
        %% Amount present at the start of simulation/initialization gas phase
        V_gas = V_tot_col*epsilon;
        y_init = [y_init_1, y_init_2, y_init_3, y_init_4];
        y_init = [y_init, 1-sum(y_init)];
        n_tot_init = P_0/R/T_0 .* V_gas .* y_init;
        %
        %% Amount present at the end of simulation gas phase
        P_sim_end = solution.P(end, 2:end-1).*1e03;
        T_sim_end = solution.T(end, 2:end-1)+273.15;
        y_sim_end = [solution.C1(end, 2:end-1); solution.C2(end, 2:end-1); solution.C3(end, 2:end-1);...
                     solution.C4(end, 2:end-1); solution.C5(end, 2:end-1)];
        dV = V_gas/N;
        dn = (P_sim_end./R./T_sim_end .* dV) .* y_sim_end;
        n_tot_final = sum(dn, 2);
        %
        %% Calculation of total adsorbed amount
        n_tot_in = trapz(t, ndot_tot_0.*y_in, 1);
        n_tot_out = trapz(t, ndot_tot_out.*y_out, 1);

        n_ads_tot = n_tot_init + n_tot_in - n_tot_out - n_tot_final';
        
        %% Adsorption loading (mol/kg)
        n_loading = n_ads_tot ./ mass_ads;

        % % Normalizing MF
        % C1 = solution.C1./y0(1);
        % C2 = solution.C2./y0(2);
        % C3 = solution.C3./y0(3);
        % C4 = solution.C4./y0(4);
        % C5 = solution.C5./y0(5);
        % 
        % [t_size, N_size] = size(C1);                
        % norm_conc = [C1(:, end), C2(:, end), C3(:, end), C4(:, end)];
        % norm_conc = norm_conc(:, 1:num_comp);
        % 
        % area_1 = zeros(1, num_comp);
        % area_2 = zeros(1, num_comp);
        % 
        % % Checking for roll up behaviour
        % for j=1:num_comp
        %     % Area before the roll-up curve
        %     idx_1 = find(norm_conc(:, j)>=1.0, 1, "first");
        %     area_1(j) = trapz(norm_conc(1:idx_1, j), t(1:idx_1));
        % 
        %     % Area under the roll-up curve
        %     idx_2 = find(norm_conc(idx_1+1:end, j)<=1.0, 1, "first");
        %     area_2(j) = trapz(t(idx_1:idx_2), norm_conc(idx_1:idx_2, j)-1.0);
        % end
        % amount_ads = area_1 - area_2;     
    end
    %% EVENT FUNCTION
    function [value, isterminal, direction] = event_function(~, ~)
        if progress_bar.CancelRequested
            value = 0;
            % error("Simulation stopped by user!!!")
        else
            value = 1;
        end
        cached_p0 = cached_p0_local;
        isterminal = 1;
        direction = 0;
    end
%
end

