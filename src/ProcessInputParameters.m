function InputParams = ProcessInputParameters(parameter_set)  
%% State all input parameters for the simulation
    max_num_comp = 5;

    num_comp = parameter_set.CompNum;
    feed_gas_props = parameter_set.FeedGas;
    column_bed = parameter_set.ColBed;
    numerical_par = parameter_set.NumPar;
    FeedFlowType = parameter_set.FeedFlowType;
    feed_MF = parameter_set.FeedMF;
    init_MF = parameter_set.InitMF;
    MTC = parameter_set.MTC;
    molar_mass = parameter_set.FeedMM;
    isotherm_params = parameter_set.IsothermParams;

    y_feed = zeros(1, max_num_comp);
    mm_feed = zeros(1, max_num_comp);
    y_init = zeros(1, max_num_comp);
    mass_trans_coef = zeros(1, max_num_comp);
    
    % Column and Feed properties and parameters
    % P_0         = column_bed(1) * 1000;                 % Adsorption pressure [kPa --> Pa]
    P_0         = column_bed(1);                        % Adsorption pressure [Pa]
    L           = column_bed(2)       ;                 % Length of the column [m]   
    dia_in      = column_bed(3)       ;                 % Internal Diameter of bed [m]
    T_wall      = column_bed(4)       ;                 % Temperature of Wall [oC]
    
    ro_total                = column_bed(6)   ;         % Solid Bulk Density of Bed
    dia_p                   = column_bed(7)   ;         % Radius of the pellets [m]
    epsilon_1               = column_bed(9)   ;         % Void fraction bed
    epsilon_2               = column_bed(10)  ;         % Void fraction particle
    epsilon_3  = epsilon_1 ...
                + (1-epsilon_1)*epsilon_2     ;         % Void fraction total
    C_psolid                = column_bed(11)  ;         % Specific heat capacity of the solid [J/kg/K]
    h_wall_gas              = column_bed(12)  ;         % Heat Transfer COefficient [W/m2/K]

    % Feed gas parameters and constants
    feed_temperature        = feed_gas_props(1)             ;   % Feed temperature of gas [oC]
    C_pg                    = feed_gas_props(3)             ;   % Specific heat of gas [J/mol/k]
    K_z                     = feed_gas_props(4)             ;   % Thermal conduction in gas phase [W/m/k]
    % v_0                     = feed_gas_props(5)             ;   % Feed Superficial velocity;
    mu                      = feed_gas_props(6)             ;   % Viscosity of gas [Pa*s]
    D_m                     = feed_gas_props(7)             ;   % Molecular diffusivity [m^2/s]
    
    y_feed(1:num_comp-1) = feed_MF(2:end)  ;       % Feed composition of adsorbing components
    y_feed(end)  = feed_MF(1);                     % Feed composition of carrier gas
    
    mm_feed(1:num_comp-1) = molar_mass(2:end)  ;   % Feed molar mass of adsorbing components
    mm_feed(end)  = molar_mass(1);                 % Feed molar mass of carrier gas

    y_init(1:num_comp-1) = init_MF(2:end)  ;       % Initial composition of adsorbing components
    y_init(end)  = init_MF(1);                     % initial composition of carrier gas

    mass_trans_coef(1:num_comp-1) = MTC;           % No carrier gas mass transfer coefficeint

    y_0_1 = y_feed(1);
    y_0_2 = y_feed(2);
    y_0_3 = y_feed(3);
    y_0_4 = y_feed(4);

    MW_1  = mm_feed(1)                            ;   % Molecular weight [kg/mol]
    MW_2  = mm_feed(2)                            ;   % Molecular weight [kg/mol]
	MW_3  = mm_feed(3)                            ;   % Molecular weight [kg/mol]
    MW_4  = mm_feed(4)                            ;   % Molecular weight [kg/mol]
    MW_5  = mm_feed(5)                            ;   % Molecular weight [kg/mol]

    y_init_1 = y_init(1);
    y_init_2 = y_init(2);
    y_init_3 = y_init(3);
    y_init_4 = y_init(4);
    
    % Feed gas parameters and constants
    R          = 8.314                          ;   % Universal gas constant [J/mol/K : Pa*m^3/mol/K]
    % T_0        = feed_temperature+273.15        ;   % Feed temperature of flue gas [K]
    T_0        = feed_temperature;              % Feed temperature [K]             
    C_tot      = P_0/R/T_0                      ;
    
    if FeedFlowType == 0
        v_0 = feed_gas_props(5);
    elseif FeedFlowType == 1
        V_std = feed_gas_props(5);          % The sccm value from user
        V_std = V_std ./ 60e6;              % cm3/min to m3/s 
        P_std = 1.01325e5;
        T_std = 273.15;
        
        % The volume flow under process condition will be calculated using ideal gas law
        V_process = V_std * (P_std/P_0) * (T_0/T_std);
        v_0 = V_process / (pi/4 * dia_in^2);
    end

    % Numerical Parameters
    N = numerical_par(1);         % No. of Nodes
    t_ads = numerical_par(2);     % Time of process    
    
    dz          = L / N                 ;   % Differntial length of the bed [m]
    
    % Temperature is already in K now.
    % T_wall      = T_wall+273.15         ;   % Temperature of column wall [K]
    
    P_inlet     = 1.02;
    
    
    ndot_0          = v_0 * C_tot               ;
    
  
    
    
    C_pa       = C_pg                           ;   % Specific heat of adsorbed phase [J/mol/k]
    feed_gas   = 'Constant Velocity'            ;   % Whether gas during the feed step has a constant pressure or velocity

    % Adsorbent parameters
    ro_s       = ro_total                       ;   % Density of the adsorbent [kg/m^3]
    r_p        = dia_p/2                        ;   % Radius of the pellets [m]
    C_ps       = C_psolid                       ;   % Specific heat capacity of the adsorbent [J/kg/K]
    q_s0       = max(isotherm_params(2:end-2, :), ...
                    [], "all")                  ;   % Molar loading scaling factor [mol/kg]
    if isnan(q_s0) || q_s0 <= 0
        q_s0 = 1;
    end
    dH_array = isotherm_params(end-1, 1:5);
    dH_array(isnan(dH_array)) = 0;
    
%% Distribute the values to the necessary variables
    Params     = zeros(54, 1) ;
    Params(1)  = N			  ;
    Params(2) = L			  ;
    Params(3) = dz           ;
    Params(4)  = T_0		  ;
    Params(5)  = mu			  ;
    Params(6)  = R			  ;
    Params(7) = v_0		  ;
    Params(8) = q_s0		  ;
    Params(9) = C_pg		  ;
    Params(10) = C_pa		  ;
    Params(11) = C_ps		  ;
    Params(12) = D_m		  ;
    Params(13) = K_z		  ;
    Params(14) = P_0		  ;
    Params(15) = T_wall	      ;
    Params(16) = dia_in	      ;
    Params(18) = P_inlet	  ;
    Params(19) = MW_1  	      ;
    Params(20) = MW_2  	      ;
    Params(21) = MW_3  	      ;
    Params(22) = MW_4  	      ;
    Params(23) = MW_5  	      ;
    Params(24) = y_0_1        ;
    Params(25) = y_0_2        ;
    Params(26) = y_0_3        ;
    Params(27) = y_0_4        ;
    Params(28)  = epsilon_1	  ;
    Params(29) = epsilon_2    ;   % Void fraction particle
    Params(30) = epsilon_3    ;   % Void fraction total bed
    Params(31)  = r_p		  ;
    Params(32)  = ro_s		  ;
    Params(36) = t_ads * v_0/L   ;
    Params(38) = ndot_0       ;
    Params(39) = h_wall_gas   ;
    Params(40) = y_init_1     ;
    Params(41) = y_init_2     ;
    Params(42) = y_init_3     ;
    Params(43) = y_init_4     ;
    Params(44:48) = mass_trans_coef;
    Params(49:53) = dH_array;
    if strcmpi(feed_gas, 'Constant Pressure') == 1
        Params(end) = 1 ;
    elseif strcmpi(feed_gas, 'Constant Velocity') == 1
        Params(end) = 0 ;
    else
        error('Please specify whether inlet velocity or pressure is constant for the feed step')
    end


%% Combine all lists into one variable that can easily be passed
    InputParams = Params          ;  
%   
end 