function solution = Isotherm_functions(num_comp, iso_params, P, T, T_flag, calc_mode)
%% This functions calculates and returns:
% 1: Pure compoenent loading for the given pressure and temperature.
% Pure component reduced grand potential the given pressure and
% temperature.

%     pure_loadings = zeros(length(P), length(iso_params(1, :)));
%     reduced_grand_potential = zeros(length(P), length(iso_params(1, :)));
    iso_type_flag = iso_params(1, :);
    R = 8.314;
    if T_flag
        dH = iso_params(end-1, 1:num_comp);
        T_ref = iso_params(end, 1:num_comp);
        % Normalzing the partial pressures based on CC equation
        P = P.*exp(-1.0 * dH/R .* (1./T - 1./T_ref));             
    end
  
    %% Pure Loading
    if calc_mode == 1
        pure_loadings = zeros(size(P));
        
        % Loop for all the components
        for i=1:num_comp    
            % SS and DS Langmuir Case 
            if iso_type_flag(i) == 1 || iso_type_flag(i) == 2
               % For site d , the index is same as Langmuir-Freundlich
                pure_loadings(:, i) = iso_params(2, i).*iso_params(3, i).*P(:, i) ./ (1+iso_params(3, i).*P(:, i)) ...
                                    + iso_params(4, i).*iso_params(5, i).*P(:, i) ./ (1+iso_params(5, i).*P(:, i));
               
            % SS and DS Langmuir-Freunclich Case
            elseif iso_type_flag(i) == 3 || iso_type_flag(i) == 4  
                pure_loadings(:, i) = iso_params(2, i).*iso_params(3, i).*P(:, i).^(iso_params(4, i)) ./ (1+iso_params(3, i).*P(:, i).^(iso_params(4, i)))...
                                    + iso_params(5, i).*iso_params(6, i).*P(:, i).^(iso_params(7, i)) ./ (1+iso_params(6, i).*P(:, i).^(iso_params(7, i)));
           
            % Quadratic Case
            elseif iso_type_flag(i) == 5
                pure_loadings(:, i) = iso_params(2, i).*(iso_params(3, i).*P(:, i) + 2.*iso_params(4, i).*(P(:, i).^2))...
                                   ./ (1 + iso_params(3, i).*P(:, i) + iso_params(4, i).*(P(:, i).^2));
                
            % Temkin Approximation Case
            elseif iso_type_flag(i) == 6 
                m = iso_params(2, i);
                b = iso_params(3, i);
                thetha = iso_params(4, i);

                lang_term = b.*P(:, i) ./ (1+b.*P(:, i)); 
                first_term = m .* lang_term;
                first_multiplier = m.*thetha .* (lang_term).^2;
                second_multiplier = lang_term - 1;
                pure_loadings(:, i) = first_term + first_multiplier .* second_multiplier;
  
            % BET Case
            elseif iso_type_flag(i) == 7
                m = iso_params(2, i);
                b_surface = iso_params(3, i);
                b_layers = iso_params(4, i);
                
                pure_loadings(:, i) = m .* b_surface .* P(:, i) ./ (1 - b_layers.*P(:, i)) ./ (1 - b_layers.*P(:, i) + b_surface.*P(:, i));
            
            % Sips Case
            elseif iso_type_flag(i) == 8
                m = iso_params(2, i);
                b = iso_params(3, i);
                n = iso_params(4, i);
                
                pure_loadings(:, i) = m.*(b.*P(:, i)).^(1/n) ./ (1+(b.*P(:, i)).^(1/n));
            
            % Toth Case
            elseif iso_type_flag(i) == 9
                m = iso_params(2, i);
                b = iso_params(3, i);
                n = iso_params(4, i);

                pure_loadings(:, i) = m.*b.*P(:, i) ./ (1+(b.*P(:, i)).^n).^(1/n);
            end    
        end
        solution = pure_loadings;
    
    %% Reduced Potential Calculation    
    else
        reduced_grand_potential = zeros(size(P));
        % Loop for all the components
        for i=1:num_comp   
            % SS and DS-Langmuir Case 
            if iso_type_flag(i) == 1 || iso_type_flag(i) == 2 
                reduced_grand_potential(:, i) = iso_params(2, i) .* log(1 + iso_params(3, i).*P(:, i)) ...
                                              + iso_params(4, i) .* log(1 + iso_params(5, i).*P(:, i));
           
            % SS and DS-Langmuir-Freunclich Case
            elseif iso_type_flag(i) == 3 || iso_type_flag(i) == 4
                tol = 1e-10;     % To avoid undefined values in case the chosen model is SS-Langmuir Freundlich
                reduced_grand_potential(:, i) = iso_params(2, i)/iso_params(4, i) .* log(1 + iso_params(3, i).*P(:, i).^iso_params(4, i)) ...
                                              + iso_params(5, i)/(iso_params(7, i)+tol) .* log(1 + iso_params(6, i).*P(:, i).^iso_params(7, i));
            
            % Quadratic Case
            elseif iso_type_flag(i) == 5
                reduced_grand_potential(:, i) = iso_params(2, i) .* log(1 + iso_params(3, i).*P(:, i) + iso_params(4, i).*(P(:, i).^2));
            
            % Temkin Approximation Case
            elseif iso_type_flag(i) == 6 
                reduced_grand_potential(:, i) = iso_params(2, i) .* (log(1 + iso_params(3, i).*P(:, i)) - 0.5*iso_params(4, i) .* (iso_params(3, i).*P(:, i)/(1 + iso_params(3, i).*P(:, i))).^2);
            
            % BET Case
            elseif iso_type_flag(i) == 7
                reduced_grand_potential(:, i) = iso_params(2, i) .* log((1 + iso_params(3, i).*P(:, i) - iso_params(4, i).*P(:, i)) ./ (1 - iso_params(4, i).*P(:, i)));
            

            % Sips Case
            elseif iso_type_flag(i) == 8
                reduced_grand_potential(:, i) = iso_params(2, i)*iso_params(4, i) .* log(1 + (iso_params(3, i).*P(:, i)).^(1/iso_params(4, i)));
            
            % Toth Case
            elseif iso_type_flag(i) == 9
                % Implemented as available in RUPTURA
                temp = iso_params(3, i).*P(:, i);
                theta_1 = temp ./ (1 + temp.^iso_params(4, i)).^(1/iso_params(4, i));
                theta_pow = theta_1.^iso_params(4, i);

                temp_psi = iso_params(2, i) .* (theta_1-(theta_1./iso_params(4, i)) .* log(1-theta_pow));

                temp1 = iso_params(2, i).*theta_1;
                temp2 = 0;

                for k=1:100     % Sum of first 100 terms
                    temp1 = temp1.*theta_pow;
                    temp2 = temp2 + iso_params(4, i);
                    temp_psi = temp_psi - temp1./(temp2.*(temp2+1));
                end
                reduced_grand_potential(:, i) = temp_psi;
            end
        end
        solution = reduced_grand_potential; 
    end
end