function Jac_S = Jacobian(N, CalcType)
%The function return sthe sparsity pattern for the ODE function being
%solved by ODE15s solver
%  Jac_s = [a_11 - - - , a_1n;
%                  -
%                  -
%          a_n1  - - -  a_nn]
%  The columns of the above matrix represents the the nodes. 
%  The rows of Jacobian matrix represents the ODE function vector.            
% 
% Input:
%       N = Number of Nodes
% Output:
%       Jac_S = Jacobian Pattern scheme for ODE function

%% Multi Node Jacobian Scheme
%  Because of upwind or downwind flow in the column, the rate of Pressure, 
%  Temperature, and mole fraction depends on at the node itself and adjacent nodes
%  Hence, Jacobian involves diagonal and diagonal adjacent elements.

    d = ones(N+2, 4);

    J_P = full(spdiags(d, -2:1, N+2, N+2)); % Pressure Jacobian Segment
    J_y_1 = J_P;                            % Comp 1 Mole fraction Jacobian Segment
    J_y_2 = J_P;                            % Comp 2 fraction Jacobian Segment
    J_y_3  = J_P;                           % Comp 3 fraction Jacobian Segment
    J_y_4 = J_P;                            % Comp 4 Mole fraction Jacobian Segment
    J_T = J_P;                              % Temperature Jacobian Segment
    J_Tw = J_P;                             % Wall Temperature Jacobian Segment
    
%% Single Node Jacobian Scheme    
    % Single Node Jacobian Scheme is used for molar adsorption loading as    
    % the molar loading depends only on the single local node
    J_x = full(spdiags(ones(N+2, 1), 0, N+2, N+2));
    J_x(1, 1) = 0;
    J_x(N+2, N+2) = 0;
%% Zero Node Jacobian Term
    J_0 = zeros(N+2);

%% Overall Jacobian Scheme
    if CalcType
        Jac_S = [J_P, J_y_1, J_y_2, J_y_3, J_y_4, J_x, J_x, J_x, J_x, J_x, J_T, J_0;
                 J_P, J_y_1, J_y_2, J_y_3, J_y_4, J_x, J_x, J_x, J_x, J_x, J_T, J_0;
                 J_P, J_y_1, J_y_2, J_y_3, J_y_4, J_x, J_x, J_x, J_x, J_x, J_T, J_0;
                 J_P, J_y_1, J_y_2, J_y_3, J_y_4, J_x, J_x, J_x, J_x, J_x, J_T, J_0;
                 J_P, J_y_1, J_y_2, J_y_3, J_y_4, J_x, J_x, J_x, J_x, J_x, J_T, J_0;
                 J_x, J_x, J_x, J_x, J_x, J_x, J_0, J_0, J_0, J_0, J_T, J_0;
                 J_x, J_x, J_x, J_x, J_x, J_0, J_x, J_0, J_0, J_0, J_T, J_0;
                 J_x, J_x, J_x, J_x, J_x, J_0, J_0, J_x, J_0, J_0, J_T, J_0;
                 J_x, J_x, J_x, J_x, J_x, J_0, J_0, J_0, J_x, J_0, J_T, J_0;
                 J_x, J_x, J_x, J_x, J_x, J_0, J_0, J_0, J_0, J_x, J_T, J_0;
                 J_P, J_y_1, J_y_2, J_y_3, J_y_4, J_x, J_x, J_x, J_x, J_x, J_T, J_Tw;
                 J_0, J_0, J_0, J_0, J_0, J_0, J_0, J_0, J_0, J_0, J_T, J_Tw];
    else
        Jac_S = [J_P, J_y_1, J_y_2, J_y_3, J_y_4, J_x, J_x, J_x, J_x, J_x;
                 J_P, J_y_1, J_y_2, J_y_3, J_y_4, J_x, J_x, J_x, J_x, J_x;
                 J_P, J_y_1, J_y_2, J_y_3, J_y_4, J_x, J_x, J_x, J_x, J_x;
                 J_P, J_y_1, J_y_2, J_y_3, J_y_4, J_x, J_x, J_x, J_x, J_x;
                 J_P, J_y_1, J_y_2, J_y_3, J_y_4, J_x, J_x, J_x, J_x, J_x;
                 J_x, J_x, J_x, J_x, J_x, J_x, J_0, J_0, J_0, J_0;
                 J_x, J_x, J_x, J_x, J_x, J_0, J_x, J_0, J_0, J_0;
                 J_x, J_x, J_x, J_x, J_x, J_0, J_0, J_x, J_0, J_0;
                 J_x, J_x, J_x, J_x, J_x, J_0, J_0, J_0, J_x, J_0;
                 J_x, J_x, J_x, J_x, J_x, J_0, J_0, J_0, J_0, J_x];
    end
%% Jacobian Scheme Modification based on BCs
%   The boundaty nodes are invariant with time.

%   Pressure Inlet 
    Jac_S(1,:) = 0;
    Jac_S(:,1) = 0;

%   Pressure Oulet 
    Jac_S(N+2,:) = 0;
    Jac_S(:, N+2) = 0;

%   Comp 1 Inlet
    Jac_S(N+3,:) = 0;
    Jac_S(:, N+3) = 0;

%   Comp 1 Oulet
    Jac_S(2*N+4,:) = Jac_S(2*N+3);
    Jac_S(:, 2*N+4) = 0;

%   Comp 2 Inlet
    Jac_S(2*N+5,:) = 0;
    Jac_S(:, 2*N+5) = 0;

%   Comp 2 Oulet
    Jac_S(3*N+6,:) = Jac_S(3*N+5,:);
    Jac_S(:, 3*N+6) = 0;

%   Comp 3 Inlet
    Jac_S(3*N+7,:) = 0;
    Jac_S(:, 3*N+7) = 0;

%   Comp 3 Oulet 
    Jac_S(4*N+8,:) = Jac_S(4*N+7,:);
    Jac_S(:,4*N+8) = 0;

%   Comp 4 Inlet 
    Jac_S(4*N+9,:) = 0;
    Jac_S(:, 4*N+9) = 0;

%   Comp 4 Oulet
    Jac_S(5*N+10,:) = Jac_S(5*N+9,:);
    Jac_S(:, 5*N+10) = 0;

    if CalcType
    %   Temperature Inlet
        Jac_S(10*N+21,:) = 0;
        Jac_S(:, 10*N+21) = 0;
    
    %   Temperature Outlet
        Jac_S(11*N+22,:) = Jac_S(11*N+21,:);
        Jac_S(:, 11*N+22) = 0;
    
    %   Wall Temperature Inlet
        Jac_S(11*N+23,:) = 0;
        Jac_S(:, 11*N+23) = 0;
    
    %   Wall Temperature Outlet
        Jac_S(12*N+24,:) = 0;
        Jac_S(:, 12*N+24) = 0;
    end

    Jac_S = sparse(Jac_S);
end

