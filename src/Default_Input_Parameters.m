function default_params = Default_Input_Parameters()   
    %% Column and Bed Parameters 
    L           = 1.0                 ;   % Length of the column [m]   
    P_0         = 101325              ;   % Adsorption pressure [Pa]
    dia_in      = 0.01                ;   % Internal Diameter of column [m]
    T_wall      = 298.15              ;   % Temperature of column wall [K]
    
    ro_total   = 500                  ;   % Solid Bulk Density of Bed [kg/m3]
    dia_p      = 2*2.5e-04            ;   % Radius of the pellets [m]
    epsilon_b  = 0.48                 ;   % Void fraction bed
    epsilon_p  = 0.24                 ;   % Void fraction particle
    C_psolid   = 850                  ;   % Specific heat capacity of the solid [J/kg/K]
    h_in       = 71                   ;   % Heat Transfer Coefficient for Wall to gas heat transfer [W/m2/K]
    
    K_wall     = 0.1                  ;   % Thermal conductivity of wall [W/m/K]
    ro_wall    = 1000                 ;   % Mass density of column wall [kg/m3]
    C_pwall    = 500                  ;   % Specific heat capacity of wall [J/kg/K]
    dia_out    = 0.015                ;   % Outer diamater of column [m]
    h_out      = 10                   ;   % Heat Transfer Coefficient for Wall to ambient heat transfer [W/m2/K] 
    T_amb      = 298.15               ;   % Ambient temperature [K]
    %
    %% Number and names of Components
    num_comp = 2;
    name_comp = {'Carrier Gas'};
    for i=2:num_comp
        name_comp(end+1) = {sprintf('Component %d', i-1)};
    end
    %
    %% Feed gas parameters and constants
    T_0        = 298.15                         ;   % Feed temperature [oC]
    C_pg       = 42                             ;   % Specific heat of gas [J/mol/k]
    K_z        = 0.29                           ;   % Thermal conduction in gas phase [W/m/k]
    v_0        = 0.1385                         ;   % Feed Superficial velocity;
    mu         = 2.87e-5                        ;   % Viscosity of gas [Pa*s]
    D_m        = 1.60e-5                        ;   % Molecular diffusivity [m^2/s]
    %
    %% Feed Mole fraction and molar masses
    y_feed = ones(1, num_comp)./num_comp;
    mm_feed = ones(1, num_comp).*0.02;
    %
    %% Numerical Parameters
    N = 30;                         % No. of Nodes
    time = 100;                     % Time of Adsorption Process
    AbsTol = 1e-06;                 % Absolute Tolerance value for ode15s
    RelTol = 1e-06;                 % Relative Tolerance value for ode15s
    y_init = zeros(1, num_comp);    % Initial column composition
    y_init(1) = 1;                  % Initializing with carrier gas only
    %
    %% Isotherm Parameter Initialization
    num_max_params = 9;
    num_max_comp = 6;
    iso_models_avai = {'SS-Langmuir', 'DS-Langmuir',...
                       'SS-Langmuir-Freundlich', 'DS-Langmuir-Freundlich',...
                       'Quadratic', 'Temkin', 'BET', 'Sips', 'Toth'};
    isotherm_params = NaN(num_max_params, num_max_comp);

    % By default all the component isotherm model will be
    % SS-Langmuir
    index_value = find(strcmpi(iso_models_avai, 'SS-Langmuir'));
    isotherm_params(1, :) = index_value;
    %
    %% Mass Transfer Coefficients [1/s]
    MTC = 0.1.* ones(1, num_comp-1);   % Does not account carrier gas                      
    %
    %% Mixture Prediction Method
    mixture_mode = 1;               %(EDSL:0, IAST:1)
    %
    %% Calculation Type
    calc_mode = 0;                  %(Isothermal:0, Non-isothermal:1)
    %
    %% PARAMETERS PACKING
    % default_params.FeedGas = [T_0, 0, C_pg, K_z, v_0, mu, D_m, y_0_1, y_0_2, y_0_3, y_0_4, y_0_5];
    default_params.FeedGas = [T_0, 0, C_pg, K_z, v_0, mu, D_m];
    default_params.ColBed = [P_0, L, dia_in, T_wall, 0, ro_total,...
                             dia_p, 0, epsilon_b, epsilon_p, C_psolid, h_in,...
                             K_wall, ro_wall, C_pwall, dia_out, h_out, T_amb];
    
    default_params.NumPar = [N, time, AbsTol, RelTol];
    
    default_params.CompNum = num_comp;
    default_params.CompNames = name_comp;
    default_params.FeedFlowType = 0;
    default_params.FeedMF = y_feed;
    default_params.FeedMM = mm_feed;
    default_params.InitMF = y_init;
    default_params.IsothermParams = isotherm_params;
    default_params.MTC = MTC; % Removing Carrier Gas
    
    default_params.MixtureMode = mixture_mode;
    default_params.CalcMode = calc_mode;
end