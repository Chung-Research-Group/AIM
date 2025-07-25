function derivatives = Adsorption_step_EDSL_Iso(~, state_vars, Params, isotherm_params_array, comp_num)
%% Retrieve process parameters
    N       = Params(1)  			  ;
    L       = Params(2) 			  ;
    T_0     = Params(4)  		  ;
    mu      = Params(5) 			  ;
    R       = Params(6) 			  ;
    v_0     = Params(7) 		  ; 
    q_s0    = Params(8) 		  ;
    % C_pg    = Params(9) 		  ;
    % C_pa    = Params(10) 		  ;
    % C_ps    = Params(11) 		  ;
    D_m     = Params(12) 		  ;
    % K_z     = Params(13) 		  ;
    P_0     = Params(14) 		  ;
    % T_wall  = Params(15) 	      ;
    % r_in  = Params(16)       ;
    P_inlet = Params(18) 	  ;
    MW_1    = Params(19)  	      ;
    MW_2    = Params(20)  	      ;
    MW_3    = Params(21)  	      ;
    MW_4    = Params(22)  	      ;
    MW_5    = Params(23)  	      ;
    y_0_1   = Params(24)        ;
    y_0_2   = Params(25)        ;
    y_0_3   = Params(26)        ;
    y_0_4   = Params(27)        ;
    epsilon = Params(28)	  ;
%     epsilon_2 = Params(29)    ;   % Void fraction particle
    epsilon_3 = Params(30)    ;   % Void fraction total bed
    r_p = Params(31)		  ;
    ro_s =  Params(32)		  ;
    % h_in = Params(39)   ;
    mass_trans_coeff = Params(44:48);
       
    ndot_0          =   P_0/R/T_0*v_0;   
    %   
%% Initialize state variables
    P    = zeros(N+2, 1) ;
    y_1  = zeros(N+2, 1) ;  
    y_2  = zeros(N+2, 1) ;  
    y_3  = zeros(N+2, 1) ;  
    y_4  = zeros(N+2, 1) ;  
    x1   = zeros(N+2, 1) ;
    x2   = zeros(N+2, 1) ;
    x3   = zeros(N+2, 1) ;
    x4   = zeros(N+2, 1) ;
    x5   = zeros(N+2, 1) ;
    % T    = zeros(N+2, 1) ;
    T    = ones(N+2, 1);

    P(1:N+2)    = state_vars(1:N+2)                 ;
    y_1(1:N+2)  = min(max(state_vars(N+3:2*N+4), 0), 1)     ;
    y_2(1:N+2)  = min(max(state_vars(2*N+5:3*N+6), 0), 1)   ;
    y_3(1:N+2)  = min(max(state_vars(3*N+7:4*N+8), 0), 1)   ;
    y_4(1:N+2)  = min(max(state_vars(4*N+9:5*N+10), 0), 1)  ;
    x1(1:N+2) = max(state_vars(5*N+11:6*N+12), 0)   ;
    x2(1:N+2) = max(state_vars(6*N+13:7*N+14), 0)   ;
    x3(1:N+2) = max(state_vars(7*N+15:8*N+16), 0)   ;
    x4(1:N+2) = max(state_vars(8*N+17:9*N+18), 0)   ;
    x5(1:N+2) = max(state_vars(9*N+19:10*N+20), 0)  ;
    % T(1:N+2)  = state_vars(10*N+21:11*N+22)         ;
  
%   
%% Initialize all variables used in the function
    % Temporal derivatives
    derivatives = zeros(10*N+20, 1) ;
    dPdt        = zeros(N+2, 1)    ;
    dPdt1       = zeros(N+2, 1)    ;
    dPdt2       = zeros(N+2, 1)    ;
    dPdt3       = zeros(N+2, 1)    ;
    dPdt4       = zeros(N+2, 1)    ;
    dy_1dt      = zeros(N+2, 1)    ;    
    dy_2dt      = zeros(N+2, 1)    ;    
    dy_3dt      = zeros(N+2, 1)    ;    
    dy_4dt      = zeros(N+2, 1)    ;    
    dx1dt       = zeros(N+2, 1)    ;
    dx2dt       = zeros(N+2, 1)    ;
    dx3dt       = zeros(N+2, 1)    ;
    dx4dt       = zeros(N+2, 1)    ;
    dx5dt       = zeros(N+2, 1)    ;
    dTdt        = zeros(N+2, 1)    ;
    % dTdt1       = zeros(N+2, 1)    ;
    % dTdt2       = zeros(N+2, 1)    ;
    % dTdt3       = zeros(N+2, 1)    ;
    % dTdt4       = zeros(N+2, 1)    ;
    % dTdt5       = zeros(N+2, 1)    ;
    % dTdt6       = zeros(N+2, 1)    ;
    % dTdt7       = zeros(N+2, 1)    ;

    % Spatial derivatives
    dpdz        = zeros(N+2, 1)    ;
    dpdzh       = zeros(N+1, 1)    ;
    dy_1dz        = zeros(N+2, 1)    ;
    dy_2dz        = zeros(N+2, 1)    ;
    dy_3dz        = zeros(N+2, 1)    ;
    dy_4dz        = zeros(N+2, 1)    ;
    d2y_1dz2      = zeros(N+2, 1)    ;
    d2y_2dz2      = zeros(N+2, 1)    ;
    d2y_3dz2      = zeros(N+2, 1)    ;
    d2y_4dz2      = zeros(N+2, 1)    ;
    dTdz          = zeros(N+2, 1)    ;
    % d2Tdz2        = zeros(N+2, 1)    ;
%     
%% Calculate all parameters used
    dz   = 1/N                                ;
    D_l  = 0.7*D_m + v_0*r_p  ;                
    Pe   = v_0*L/D_l;                          
    phi  = R*T_0*q_s0*ro_s/epsilon_3/P_0   ;      
    % sigma_ads = ro_s*C_pa*q_s0;
    % ro_g = P(1:N+2).*P_0/R./T(1:N+2)/T_0      ;
    % Pe_h = epsilon_3*v_0*L*ro_g(1)*C_pg/K_z       ;  
    MW_t = MW_1.*y_1 + MW_2.*y_2 + MW_3.*y_3...
            + MW_4.*y_4 + MW_5.*(1-y_1-y_2-y_3-y_4)  ;
%   
%% Boundary Conditions
    % INLET Boundary
    % y_1(1) = y_0_1;          
    % y_2(1) = y_0_2;          
    % y_3(1) = y_0_3;    
    % y_4(1) = y_0_4;    
%     T(1) = T_0/T_0;

    % Danckwert's Boundary Conditions
    y_1(1) = (y_1(2) + v_0/v_0/epsilon*Pe*y_0_1*0.5*dz)/(1 + v_0/v_0/epsilon*Pe*0.5*dz);              
    y_2(1) = (y_2(2) + v_0/v_0/epsilon*Pe*y_0_2*0.5*dz)/(1 + v_0/v_0/epsilon*Pe*0.5*dz);               
    y_3(1) = (y_3(2) + v_0/v_0/epsilon*Pe*y_0_3*0.5*dz)/(1 + v_0/v_0/epsilon*Pe*0.5*dz);        
    y_4(1) = (y_4(2) + v_0/v_0/epsilon*Pe*y_0_4*0.5*dz)/(1 + v_0/v_0/epsilon*Pe*0.5*dz);         
    % T(1) = (T(2) + v_0/v_0/epsilon*Pe_h*T_0/T_0*0.5*dz)/(1 + v_0/v_0/epsilon*Pe_h*0.5*dz);       
    
    if Params(end) == 1
        % Inlet pressure is specified as a constant
        P(1) = P_inlet ;
    elseif Params(end) == 0
        % Inlet velocity is specified as a constant
        vh = zeros(N+1, 1) ;
        

        a_1   = 150*mu*(1-epsilon)^2*dz*L/2/4/r_p^2/epsilon^3/T(1)/T_0/R ;
        a_2_1 = 1.75*(1-epsilon)/2/r_p/epsilon/epsilon/epsilon*dz*L/2    ;
        a_2   = a_2_1/R/T(1)/T_0*ndot_0*MW_t(1)                               ;
    
        a =  a_1+a_2             ;
        b =  P(2)/T(1)*P_0/R/T_0 ;
        c = -ndot_0              ;
    
        vh(1) = (-b+sqrt(b^2-4*a*c))/2/a/v_0 ;

        a_p = a_1*T(1)*T_0*R      ;
        b_p = a_2_1*MW_t(1)/R/T(1)/T_0 ;
    
        P(1)  = ((a_p*vh(1)*v_0+P(2)*P_0)./(1-b_p*vh(1)*v_0*vh(1)*v_0))/P_0; 
        
    else
        error('Please specify whether inlet velocity or pressure is constant for the feed step')
    end
%   

%%  OUTLET Boundary
    y_1(N+2) = y_1(N+1) ;
    y_2(N+2) = y_2(N+1) ;
    y_3(N+2) = y_3(N+1) ;
    y_4(N+2) = y_4(N+1) ;
    % T(N+2) = T(N+1)     ;
    if P(N+1) >= 1
        P(N+2) = 1      ;
    else
        P(N+2) = P(N+1) ;
    end
%   
%% Spatial Derivative Calculations
%%   Pressure: at the center of the volumes and the walls
    
    Ph          = WENO(P, 'upwind')      ;
    
    dpdz(2:N+1) = (Ph(2:N+1)-Ph(1:N))/dz ;
    dpdzh(2:N)  = (P(3:N+1)-P(2:N))/dz   ;
    dpdzh(1)    =   2*(P(2)-P(1))/dz      ;
    dpdzh(N+1)  =   2*(P(N+2)-P(N+1))/dz  ;
%   
%%  Mole Fraction calculation for chemical species at the center of the volumes and walls  
%   1st Component
    der_out = spatial_der(y_1, dy_1dz, d2y_1dz2, N, dz) ;
    y_1h     = der_out{1}                            ;
    dy_1dz   = der_out{2}                               ;
    d2y_1dz2 = der_out{3}                               ;
%     
%   2nd Component   
    der_out  = spatial_der(y_2, dy_2dz, d2y_2dz2, N, dz) ;
    y_2h     = der_out{1}                                ;
    dy_2dz   = der_out{2}                                ;
    d2y_2dz2 = der_out{3}                                ;
    
%   3rd Component      
    der_out  = spatial_der(y_3, dy_3dz, d2y_3dz2, N, dz) ;
    y_3h     = der_out{1}   ;
    dy_3dz   = der_out{2}   ;
    d2y_3dz2 = der_out{3}   ;
    
%   4th Component          
    der_out  = spatial_der(y_4, dy_4dz, d2y_4dz2, N, dz) ;
    y_4h     = der_out{1}   ;
    dy_4dz   = der_out{2}   ;
    d2y_4dz2 = der_out{3}   ;
        
%   
%%  Temperature: at the center of the volumes
%   First Derivative
    
    % Th          = WENO(T, 'upwind')      ;
    Th = ones(size(Ph));
    % dTdz(2:N+1) = (Th(2:N+1)-Th(1:N))/dz ;

%   Second Derivative
    % d2Tdz2(3:N) = (T(4:N+1)+T(2:N-1)-2*T(3:N))/dz/dz ;
    % d2Tdz2(2)   =  4*(Th(2)+T(1)-2*T(2))/dz/dz       ;
    % d2Tdz2(N+1) =  4*(Th(N)+T(N+2)-2*T(N+1))/dz/dz   ;
%   
%% Velocity Calculations
    ro_gh          = (P_0/R/T_0)*Ph(1:N+1)./Th(1:N+1)               ;
    
    viscous_term   = 150*mu*(1-epsilon)^2/4/r_p^2/epsilon^2         ;
    kinetic_term_h = (ro_gh.*(MW_1.*y_1h + MW_2.*y_2h + MW_3.*y_3h...
                      + MW_4.*y_4h + MW_5.*(1-y_1h-y_2h-y_3h-y_4h)))...
                     .*(1.75*(1-epsilon)/2/r_p/epsilon)             ; 
                                                                       
    % Velocities at walls of volumes 
    vh = -sign(dpdzh).*(-viscous_term+(abs(viscous_term^2+4*kinetic_term_h... 
               .*abs(dpdzh)*P_0/L)).^(.5))/2./kinetic_term_h./v_0 ;   
    vh(1, 1);

%% Adsoption Term
    y_array = [y_1, y_2, y_3, y_4];
    eq_loading = zeros(size(y_1, 1),5);
    eq_loading(:, 1:comp_num-1) = EDSL_func(comp_num-1, P.*P_0, y_array(:, 1:comp_num-1), T.*T_0, isotherm_params_array, 0);
%% Temporal Derivatives
%   1.1) Calculate equilibrium molar loading for species
    q_eq_1   = eq_loading(:, 1);
    q_eq_2   = eq_loading(:, 2);
    q_eq_3   = eq_loading(:, 3);
    q_eq_4   = eq_loading(:, 4);
    q_eq_5   = eq_loading(:, 5);
%   
%%  
%   1.2) Calculate the LDF parameter
    k_1 = mass_trans_coeff(1) .* ones(length(y_1), 1);
    k_2 = mass_trans_coeff(2) .* ones(length(y_1), 1);
    k_3 = mass_trans_coeff(3) .* ones(length(y_1), 1);
    k_4 = mass_trans_coeff(4) .* ones(length(y_1), 1);
    k_5 = mass_trans_coeff(5) .* ones(length(y_1), 1);
%
%% 
    % Heat of Adsoprtion
%     deltaU_array = zeros(5, 1);
    % deltaU_array = isotherm_params_array(4, :);
    % 
    % deltaU_1 = deltaU_array(1);
    % deltaU_2 = deltaU_array(2);
    % deltaU_3 = deltaU_array(3);
    % deltaU_4 = deltaU_array(4);
    % deltaU_5 = deltaU_array(5);
    
    % Getting the delta H values from isotherm params array
    % dH_array = isotherm_params_array(end-1, 1:5);
    
%%  
%   1.3) Calculate the temporal derivative 
    dx1dt(2:N+1) = k_1(2:N+1) .* (q_eq_1(2:N+1) ./ q_s0 - x1(2:N+1)) .* L /v_0 ;
    dx2dt(2:N+1) = k_2(2:N+1) .* (q_eq_2(2:N+1) ./ q_s0 - x2(2:N+1)) .* L /v_0 ;
    dx3dt(2:N+1) = k_3(2:N+1) .* (q_eq_3(2:N+1) ./ q_s0 - x3(2:N+1)) .* L /v_0 ;
    dx4dt(2:N+1) = k_4(2:N+1) .* (q_eq_4(2:N+1) ./ q_s0 - x4(2:N+1)) .* L /v_0 ;   
    dx5dt(2:N+1) = k_5(2:N+1) .* (q_eq_5(2:N+1) ./ q_s0 - x5(2:N+1)) .* L /v_0 ;  
    
%%  Energy Balance

    % sink_term = (epsilon_3 .*ro_g(2:N+1) .* C_pg + ro_s * C_ps + ro_s * q_s0 * C_pa .* ...
    %                 (x1(2: N+1) + x2(2: N+1) + x3(2: N+1) + x4(2: N+1)...
    %                     + x5(2: N+1)));  
  
    % transfer_term = K_z ./ v_0 ./ L                             ;
    % dTdt1(2:N+1)  = transfer_term .* d2Tdz2(2:N+1) ./ sink_term ;
 
%   2.2) Calculate the temperature change due to advection
  
    PvT          = Ph(1:N+1).*vh(1:N+1)./Th(1:N+1) ;
    % Pv           = Ph(1:N+1).*vh(1:N+1)            ;
%     dTdt2(2:N+1) = -epsilon.*C_pg.*P_0./R./T_0.*((Pv(2:N+1)-Pv(1:N))- ... 
%                     T(2:N+1).*(PvT(2:N+1)-PvT(1:N)))./dz./sink_term      ;
%     dTdt2(2:N+1) = - epsilon .* C_pg .* P_0 .* (Pv(2:N+1)-Pv(1:N)) ./dz ...
%                     ./ (T_0 * R) ./ sink_term    ;
%    
%%   2.3) Calculate the temperature change due to adsorption/desoprtion enthalpy  
    % generation_term_1 = -(deltaU_1-R.* T .* T_0);
    % generation_term_2 = -(deltaU_2-R.* T .* T_0);
    % generation_term_3 = -(deltaU_3-R.* T .* T_0);
    % generation_term_4 = -(deltaU_4-R.* T .* T_0);
    % generation_term_5 = -(deltaU_5-R.* T .* T_0);

    % generation_term_1 = (dH_array(1));
    % generation_term_2 = (dH_array(2));
    % generation_term_3 = (dH_array(3));
    % generation_term_4 = (dH_array(4));
    % generation_term_5 = (dH_array(5));

    % dTdt3(2:N+1)      = (ro_s .* q_s0 ./ T_0) .* (generation_term_1(2:N+1).* dx1dt(2:N+1) +...
    %                      generation_term_2(2:N+1).* dx2dt(2:N+1) + ...
    %                      generation_term_3(2:N+1).* dx3dt(2:N+1) +...
    %                      generation_term_4(2:N+1).* dx4dt(2:N+1) + ...
    %                      generation_term_5(2:N+1).* dx5dt(2:N+1))./sink_term     ;
    % dTdt3(2:N+1)      = (ro_s .* q_s0 ./ T_0) .* (generation_term_1.* dx1dt(2:N+1) +...
    %                      generation_term_2.* dx2dt(2:N+1) + ...
    %                      generation_term_3.* dx3dt(2:N+1) +...
    %                      generation_term_4.* dx4dt(2:N+1) + ...
    %                      generation_term_5.* dx5dt(2:N+1))./sink_term     ; 
%    
%%  2.4) Temperature Change due to mass capacity change of gas

    % dTdt5(2:N+1) = -C_pg .* epsilon_3 .* P_0 .* dPdt(2:N+1)...
    %                ./ (R * T_0) ./ sink_term;
%%  2.5) Temperature change due to heat transfer to wall
    % T_w = (T_wall)/T_0  ;
    % dTdt6(2:N+1) = (2 * h_in * L / r_in /v_0) .* (T_w - T(2:N+1)) ./ sink_term;
    

%%  2.7) Temperature change due to mass transfer from gas to adsorbent
    % dTdt7(2:N+1) = - sigma_ads .* T(2:N+1) .* (dx1dt(2:N+1) + dx2dt(2:N+1) + ...
    %                        dx3dt(2:N+1) + dx4dt(2:N+1) + dx5dt(2:N+1)) ./ sink_term;
%%   2.8) Total sum of all temperature derivatives 
    
    % dTdt(2:N+1) = dTdt1(2:N+1) + dTdt2(2:N+1) + dTdt3(2:N+1) + dTdt4(2:N+1)...
    %                + dTdt5(2:N+1) + dTdt6(2:N+1) + dTdt7(2:N+1);
    
      dTdt(2:N+1) = 0;
     
%%  
%   3) Total mass balance
%%  
%   3.1) Calculate the change in pressure due to advection
    dPdt1(2:N+1) = -(epsilon / epsilon_3) .* T(2:N+1).*(PvT(2:N+1)-PvT(1:N))./dz  ;
%   
%%  
%   3.2) Calculate the change in pressure due to adsorption/desorption
    dPdt2(2:N+1) = -phi .* T(2:N+1).*(dx1dt(2:N+1) + dx2dt(2:N+1) + ...
                           dx3dt(2:N+1) + dx4dt(2:N+1) + dx5dt(2:N+1)) ;
%   
%%  
%   3.3) Calculate the change in pressure due to temperature changes
    dPdt3(2:N+1) = P(2:N+1).*dTdt(2:N+1)./T(2:N+1)            ;
%                       
%%
%   3.5) Total sum of all presure changes
    dPdt(2:N+1) = dPdt1(2:N+1) + dPdt2(2:N+1) + dPdt3(2:N+1) + dPdt4(2:N+1);
%   
%%  
 % *4) Component Mass Balance (Based on Mole Fraction)*
    dy_1dt(2:N+1) = molebalances(P, T, y_1, Ph, Th, y_1h, vh, dy_1dz...
             , d2y_1dz2, dpdz, dTdz, dPdt, dTdt, dx1dt, epsilon,...
             epsilon_3, Pe, phi, dz, N) ;

    dy_2dt(2:N+1) = molebalances(P, T, y_2, Ph, Th, y_2h, vh, dy_2dz...
             , d2y_2dz2, dpdz, dTdz, dPdt, dTdt, dx2dt, epsilon,...
             epsilon_3, Pe, phi, dz, N) ;

    dy_3dt(2:N+1) = molebalances(P, T, y_3, Ph, Th, y_3h, vh, dy_3dz...
             , d2y_3dz2, dpdz, dTdz, dPdt, dTdt, dx3dt, epsilon,...
             epsilon_3, Pe, phi, dz, N) ;

    dy_4dt(2:N+1) = molebalances(P, T, y_4, Ph, Th, y_4h, vh,...
                dy_4dz, d2y_4dz2, dpdz, dTdz, dPdt, dTdt, dx4dt,...
                epsilon, epsilon_3, Pe, phi, dz, N)  ;   

%   
%%  Boundary Derivatives
    dPdt(1)    = 0         ;
    dPdt(N+2)  = 0         ;
    dy_1dt(1)    = 0         ;
    dy_1dt(N+2)  = dy_1dt(N+1) ;
    dy_2dt(1)    = 0         ;
    dy_2dt(N+2)  = dy_2dt(N+1) ;
    dy_3dt(1)    = 0         ;
    dy_3dt(N+2)  = dy_3dt(N+1) ;
    dy_4dt(1)    = 0         ;
    dy_4dt(N+2)  = dy_4dt(N+1) ;
    dx1dt(1)   = 0         ;
    dx2dt(1)   = 0         ;
    dx3dt(1)   = 0         ;
    dx4dt(1)   = 0         ;
    dx5dt(1)   = 0         ;
    dx1dt(N+2) = 0         ;
    dx2dt(N+2) = 0         ;
    dx3dt(N+2) = 0         ;
    dx4dt(N+2) = 0         ;
    dx5dt(N+2) = 0         ;
    % dTdt(1)    = 0         ;
    % dTdt(N+2)  = dTdt(N+1) ;
%   
%%  Export derivatives to output
    derivatives(1:N+2)              = dPdt(1:N+2)  ;
    derivatives(N+3:2*N+4)          = dy_1dt(1:N+2);
    derivatives(2*N+5:3*N+6)        = dy_2dt(1:N+2);
    derivatives(3*N+7:4*N+8)        = dy_3dt(1:N+2);
    derivatives(4*N+9:5*N+10)       = dy_4dt(1:N+2);
    derivatives(5*N+11:6*N+12)      = dx1dt(1:N+2) ;
    derivatives(6*N+13:7*N+14)      = dx2dt(1:N+2) ;
    derivatives(7*N+15:8*N+16)      = dx3dt(1:N+2) ;
    derivatives(8*N+17:9*N+18)      = dx4dt(1:N+2) ;
    derivatives(9*N+19:10*N+20)     = dx5dt(1:N+2) ;
    % derivatives(10*N+21:11*N+22)    = dTdt(1:N+2)  ;
%
%% Supplementary Function
    function der_out = spatial_der(y, dydz, d2ydz2, N, dz)
        % First Derivative
        yh          = WENO(y, 'upwind')      ;
        dydz(2:N+1) = (yh(2:N+1)-yh(1:N))/dz ;
    
            % Second Derivative
        d2ydz2(3:N) = (y(4:N+1)+y(2:N-1)-2*y(3:N))/dz/dz ;
        d2ydz2(2)   = (y(3)-y(2))/dz/dz                  ;
        d2ydz2(N+1) = (y(N)-y(N+1))/dz/dz                ;
       
        der_out = {yh, dydz, d2ydz2};                                               
    end
end 